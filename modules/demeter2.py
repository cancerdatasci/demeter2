import sys
import os
import time
from json_tricks.np import dump, load
from functools import reduce
import numpy as np
import tensorflow as tf
from sklearn.linear_model import LinearRegression
# from scipy.sparse import hstack, csr_matrix, csr
import pandas as pd

import edward as ed
from edward.models import Normal

if "../modules" not in sys.path:
    sys.path.append("../modules")
from preprocess import *

from ScipyOptimizerInterface import ScipyOptimizerInterface 


def run_sklearn_optim(optimizer, feed_dict, sess, loss, print_freq = 10):
    '''Run sklearn optimizer and keep track of loss.
    INPUTS:
        optimizer: optimizer op
        feed_dict: 
        sess: tf session
        loss: loss op
        print_freq: print loss per n iters
    OUTPUTS:
        dict of info on optimization results'''
    global_cnt = 0
    def callback(loss):
        nonlocal global_cnt
        if global_cnt % print_freq == 0:
            print(loss)
        sys.stdout.flush()
        global_cnt += 1
    results = optimizer.minimize(sess, feed_dict = feed_dict, fetches = [loss], loss_callback = callback)    
    return(results)


def make_sparse_tensor(csr_mat):
    '''Take a sparse matrix in csr format and makes a tf.SparseTensor'''
    coo_mat = csr_mat.tocoo()
    inds = np.concatenate([coo_mat.row[:,None], coo_mat.col[:,None]], axis = 1)
    vals = tf.to_float(coo_mat.data)
    sp_tens = tf.SparseTensor(indices=inds, values=vals, dense_shape=coo_mat.shape)
    return(sp_tens)


def update_param_dict(defaults, new_vals):
    '''Take a default dict and a dict of values to override and return the updated dict'''
    if new_vals is not None:
        assert(all(k in defaults.keys() for k in new_vals.keys()))
        defaults.update(new_vals)
    return(defaults)


def merge_dicts(orig_dict, add_dict):
    '''Update a dict with key-value pairs from a new dict'''
    new_dict = orig_dict.copy()
    new_dict.update(add_dict)
    return(new_dict)


def SSMD(pop1, pop2):
    '''Calculate SSMD between two samples'''
    beta = np.nanmean(pop1) - np.nanmean(pop2)
    beta = beta / np.sqrt(np.nanstd(pop1) ** 2 + np.nanstd(pop2) ** 2)
    return(beta)


def make_eval_masks(LFC_mats, test_ind_sets, inverse = False):
    '''Make boolean array to mask data used for training'''
    eval_masks = []
    if test_ind_sets is None:
        test_ind_sets = [np.array([], dtype = int) for _ in range(len(LFC_mats))]
    for LFC_mat, test_inds in zip(LFC_mats, test_ind_sets):
        cur_eval_mask = np.ones_like(LFC_mat, dtype=bool)
        cur_eval_mask[test_inds] = False
        if inverse: #if you want to evaluate on the test set 
            cur_eval_mask = ~cur_eval_mask
        cur_eval_mask[np.isnan(LFC_mat)] = False
        eval_masks.append(cur_eval_mask)
    return(eval_masks)


def compute_hairpins_per_gene_CL(LFC_mats, sparse_mat, unique_hp_seqs, unique_CLs, unique_genes):
    '''Estimate number of measured hairpin LFCs per gene/CL or seed/CL
    INPUTS:
        LFC_mats: list of hairpin LFC mats
        sparse_mat: sparse matrix mapping hairpins to genes/seeds
        unique_hp_seqs: ordered list of unique hairpins
        unique_CLs: ordered list of unique CLs
        unique_genes: ordered list of unique genes
    OUTPUTS:
        matrix with number used hairpins per gene/CL
    '''
    hp_ind_map = {name: ind for ind, name in enumerate(unique_hp_seqs)}
    CL_ind_map = {name: ind for ind, name in enumerate(unique_CLs)}
    n_hps_per_gene = np.zeros((len(unique_genes), len(unique_CLs)))
    for LFC_mat in LFC_mats:  
        cur_hp_set = [hp_ind_map[x] for x in LFC_mat.index.values]
        cur_CL_set = [CL_ind_map[x] for x in LFC_mat.columns.values]
        cur_hps_per_gene = sparse_mat[cur_hp_set,:].transpose().dot(~np.isnan(LFC_mat))
        n_hps_per_gene[:, cur_CL_set] = n_hps_per_gene[:, cur_CL_set] + cur_hps_per_gene
    return(n_hps_per_gene)


def map_effects(scores, CL_inds, sparse_mat):
    '''Apply hairpin mapping to predicted scores'''
    cur_scores = tf.gather(scores, CL_inds)
    return(tf.sparse_tensor_dense_matmul(sparse_mat, cur_scores, adjoint_a = False, adjoint_b = True))


#*******************  DEFINE DEMETER2 MODEl CLASS *************************#
class demeter:
    '''Class implementing a DEMETER2 model'''
    def default_reg_params(self):
        params = {
                'hairpin_l2_lambda': 0,
                'hp_unpred_l2_lambda': 0,
                'CL_l2_lambda': 0,
                'gene_l2_lambda': 0, #L2 penalty on across-CL avg
                'rel_gene_l2_lambda': 0, #L2 penalty on deviation from mean
                'seed_l2_lambda': 0,
                'rel_seed_l2_lambda': 0 #L2 penalty on deviation from mean
                }
        return(params)

    def default_optim_params(self):    
        params = {'precision': 'double',
                  'maxiter': 2000,
                  'print_freq': 50,
                  'ftol': 1e-7}
        return(params)


    def __init__(self, LFC_mats, gene_matrix, seed_matrix, gene_sets, data_names = None, 
        reg_params = None, optim_params = None, test_inds = None, log_file = None):
        '''
        Create a demeter model instance
        INPUTS:
            LFC_mats: List of matrices [hairpins x CLs] of observed LFC values, one per batch/dataset
            gene_matrix: [n_hairpins, n_genes] gene-target mapping, as a scipy csr sparse matrix. 
            seed_matrix: [n_hairpins, n_seeds] seed-target mapping, as a scipy csr sparse matrix. 
            gene_sets: dict with two entries 'pos' and 'neg'. Each are arrays of Gene names specifying positive and negative control gene sets respectively
            data_names: dict of names for different entities (genes, CLs, hps)
            reg_params: dict of regularization parameters. 
                Specify optional lambdas [hairpin_l2_lambda, CL_l2_lambda, gene_l2_lambda, seed_l2_lambda, rel_gene_l2_lambda, rel_seed_l2_lambda]
            optim_params: dict of optimization parameters
            test_inds: list of tuples specifying indices in the LFC data matrices to set aside for testing (set to None if not using xval)
            log_file: path of log file
        '''
        self.min_hairpins_per = 2 #minimum number of hairpins per gene/seed to use for estimation of gene/seed effects
        self.min_slope = 0.01 #minimum slope term (prevents them from getting set to 0 during optimization)
        reg_params = update_param_dict(self.default_reg_params(), reg_params)
        self.reg_params = reg_params
        optim_params = update_param_dict(self.default_optim_params(), optim_params)
        self.optim_params = optim_params
        if data_names is not None:
            self.data_names = data_names

        self.log_file = log_file
        if self.log_file is not None:
            self._log_file = open(log_file, 'w')
        else:
            self._log_file = None
 
        #init containers for storing stats across training iters
        self.R2_vals = {'train': [], 'test': [], 'train_ms': [], 'test_ms': []} #store R2 evals in a dict
        self.loss_evals = []
        self.SSMD = {'train': [], 'test': []}

        self.gene_sets = gene_sets

        if self.optim_params['precision'] == 'double':
            self.float = tf.float64
        elif self.optim_params['precision'] == 'single':
            self.float = tf.float32
        else:
            raise('invalid float type')

        self.test_inds = test_inds

        self.all_CL_names = get_CL_names(LFC_mats)
        self.all_CL_batches = get_CL_batches(LFC_mats)
        self.all_hp_seqs = get_hp_names(LFC_mats)
        self.all_hp_batches = get_hp_batches(LFC_mats)
        self.n_CLs = len(data_names['CLs'])
        self.n_CL_batches = len(self.all_CL_names)
        self.n_hp_batches = len(self.all_hp_seqs)


        #BUILD GRAPH
        self.g = tf.Graph()
        self.sess = tf.Session(graph = self.g)        
        with self.g.as_default():
            self.n_hairpins, self.n_genes = gene_matrix.shape
            _, self.n_seeds = seed_matrix.shape

            #calculate number of genes and seeds with data for each CL
            self.n_used_hairpins_per_gene = compute_hairpins_per_gene_CL(
                LFC_mats, gene_matrix, data_names['hps'], data_names['CLs'], data_names['genes'])
            self.n_targeted_genes = np.sum(self.n_used_hairpins_per_gene >= self.min_hairpins_per, axis = 0) 
            self.n_used_hairpins_per_seed = compute_hairpins_per_gene_CL(
                LFC_mats, seed_matrix, data_names['hps'], data_names['CLs'], data_names['seeds'])
            self.n_targeted_seeds = np.sum(self.n_used_hairpins_per_seed >= self.min_hairpins_per, axis = 0) 

            #define parameter inits
            init_params = {
                        'gene_score': tf.zeros([self.n_CLs, self.n_genes], self.float),
                        'seed_score': tf.zeros([self.n_CLs, self.n_seeds], self.float),
                        'gene_score_avgs': tf.zeros([1, self.n_genes], self.float),
                        'seed_score_avgs': tf.zeros([1, self.n_seeds], self.float),
                        'CL_offset': tf.zeros([self.n_CL_batches, 1], self.float),
                        'CL_slope': tf.ones([self.n_CL_batches, 1], self.float),
                        'gene_slope': tf.ones([self.n_CLs, 1], self.float),
                        'CL_noise_vars': tf.ones([self.n_CL_batches, 1], self.float),
                        'hairpin_offset': tf.zeros([self.n_hp_batches, 1], self.float),
                        'hairpin_unpred': tf.zeros([self.n_hairpins, 1], self.float),
                        'guide_Geff': tf.ones([self.n_hairpins, 1], self.float),
                        'guide_Seff': tf.ones([self.n_hairpins, 1], self.float)
                        }
 
            self.obs = [tf.placeholder(self.float, dset.shape, name = "obs_" + str(ii)) \
                        for ii, dset in enumerate(LFC_mats)]
            self.eval_mask = [tf.placeholder('bool', dset.shape, name = "eval_mask_" + str(ii)) \
                        for ii, dset in enumerate(LFC_mats)]

            #Define variables
            self.gene_score = tf.Variable(init_params['gene_score'], dtype = self.float, name = 'gene_score')
            self.seed_score = tf.Variable(init_params['seed_score'], dtype = self.float, name = 'seed_score')
            self.gene_score_avgs = tf.Variable(init_params['gene_score_avgs'], dtype = self.float, name = 'gene_score_avgs')
            self.seed_score_avgs = tf.Variable(init_params['seed_score_avgs'], dtype = self.float, name = 'seed_score_avgs')
            self.CL_offset = tf.Variable(init_params['CL_offset'], dtype = self.float, name = 'CL_offset')
            self.CL_slope = tf.Variable(init_params['CL_slope'], dtype = self.float, name = 'CL_slope')
            self.gene_slope = tf.Variable(init_params['gene_slope'], dtype = self.float, name = 'gene_slope')
            self.hairpin_offset = tf.Variable(init_params['hairpin_offset'], dtype = self.float, name = 'hairpin_offset')
            self.hairpin_unpred = tf.Variable(init_params['hairpin_unpred'], dtype = self.float, name = 'hairpin_offset')
            self.guide_Geff = tf.Variable(init_params['guide_Geff'], dtype = self.float, name = 'guide_Geff')
            self.guide_Seff = tf.Variable(init_params['guide_Seff'], dtype = self.float, name = 'guide_Seff')

            self.CL_noise_vars = tf.Variable(init_params['CL_noise_vars'], dtype = self.float, name = 'noise_vars')

            self.n_Geffs = self.n_hairpins
            self.n_Seffs = self.n_hairpins

            #maps from name to index value
            self.hp_ind_map = {name: ind for ind, name in enumerate(data_names['hps'])}
            self.CL_ind_map = {name: ind for ind, name in enumerate(data_names['CLs'])}

            #make list of sparse gene and seed maps for each LFC dataset
            gene_maps = [self.make_sparse_submap(gene_matrix, LFC_mat.index.values) for LFC_mat in LFC_mats]
            seed_maps = [self.make_sparse_submap(seed_matrix, LFC_mat.index.values) for LFC_mat in LFC_mats]

            #op that is the per-CL gene effect scaled by gene-KD slope (used for re-estimating gene slope)
            self.ind_gene_effects = tf.multiply(self.gene_score_avgs + self.gene_score, self.gene_slope)

            #package a dict of the model params
            mod_params = {
                'CL_noise_vars': self.CL_noise_vars,
                'CL_offset': self.CL_offset,
                'CL_slope': self.CL_slope,
                'guide_Seff': self.guide_Seff,
                'guide_Geff': self.guide_Geff,
                'hairpin_offset': self.hairpin_offset,
                'hairpin_unpred': self.hairpin_unpred,
                'gene_slope': self.gene_slope,
                'seed_score_avgs': self.seed_score_avgs,
                'seed_score': self.seed_score,
                'gene_score_avgs': self.gene_score_avgs,
                'gene_score': self.gene_score}

            #LOOP OVER DATASETS AND BUILD SUBGRAPH FOR EACH
            dataset_nLLs = []
            dataset_SS = []
            self.shRNA_R2 = []
            self.shRNA_nLL = []
            self.shRNA_oSS = []
            self.pred = []
            hp_offset = 0
            CL_offset = 0
            for ii in range(len(self.obs)):
                cur_pred = self.get_dataset_pred(
                    mod_params,
                    gene_maps[ii],
                    seed_maps[ii],
                    LFC_mats[ii].index.values,
                    LFC_mats[ii].columns.values,
                    hp_offset, 
                    CL_offset)

                cur_nLL, cur_SS = self.get_dataset_LL(
                    mod_params,
                    self.obs[ii],
                    cur_pred, 
                    self.eval_mask[ii],
                    CL_offset)
                cur_shRNA_R2, cur_shRNA_nLL, cur_shRNA_SS = self.get_shRNA_R2(
                    mod_params,
                    self.obs[ii],
                    cur_pred,
                    CL_offset)
                self.shRNA_R2.append(cur_shRNA_R2)
                self.shRNA_nLL.append(cur_shRNA_nLL)
                self.shRNA_oSS.append(cur_shRNA_SS)
                dataset_nLLs.append(cur_nLL)
                dataset_SS.append(cur_SS)
                self.pred.append(cur_pred)
                hp_offset += LFC_mats[ii].shape[0]
                CL_offset += LFC_mats[ii].shape[1]
            self.nLL = tf.add_n(dataset_nLLs) #sum negative log-like across datasets
            tot_SS = tf.add_n(dataset_SS) #sum squared error

            #LOOP OVER DATASETS AND BUILD GENE-AVG SUBGRAPHS
            dataset_avg_nLLs = []
            hp_offset = 0
            CL_offset = 0
            for ii in range(len(self.obs)):
                cur_pred = self.get_dataset_pred(
                    mod_params,
                    gene_maps[ii],
                    seed_maps[ii],
                    LFC_mats[ii].index.values,
                    LFC_mats[ii].columns.values,
                    hp_offset, 
                    CL_offset,
                    just_avg_scores = True)
                cur_nLL, cur_SS = self.get_dataset_LL(
                    mod_params,
                    self.obs[ii],
                    cur_pred, 
                    self.eval_mask[ii],
                    CL_offset)
                dataset_avg_nLLs.append(cur_nLL)
                hp_offset += LFC_mats[ii].shape[0]
                CL_offset += LFC_mats[ii].shape[1]
            self.avg_effect_loss = tf.add_n(dataset_avg_nLLs)

            #calc R2
            self.R2 = 1 - self.nLL / tot_SS

            #calc regularization penalty
            self.CL_l2_lambda = tf.Variable(reg_params['CL_l2_lambda'], dtype = self.float)
            self.hairpin_l2_lambda = tf.Variable(reg_params['hairpin_l2_lambda'], dtype = self.float)
            self.hp_unpred_l2_lambda = tf.Variable(reg_params['hp_unpred_l2_lambda'], dtype = self.float)
            self.rel_gene_l2_lambda = tf.Variable(reg_params['rel_gene_l2_lambda'], dtype = self.float)
            self.rel_seed_l2_lambda = tf.Variable(reg_params['rel_seed_l2_lambda'], dtype = self.float)
            self.gene_l2_lambda = tf.Variable(reg_params['gene_l2_lambda'], dtype = self.float)
            self.seed_l2_lambda = tf.Variable(reg_params['seed_l2_lambda'], dtype = self.float)

            self.pen = 0
            self.pen += 0.5 * tf.reduce_sum(tf.pow(self.CL_offset, 2)) * self.CL_l2_lambda
            self.pen += 0.5 * tf.reduce_sum(tf.pow(self.hairpin_offset, 2)) * self.hairpin_l2_lambda
            self.pen += 0.5 * tf.reduce_sum(tf.pow(self.hairpin_unpred, 2)) * self.hp_unpred_l2_lambda
            self.pen += 0.5 * tf.reduce_sum(tf.pow(self.gene_score, 2)) * self.rel_gene_l2_lambda
            self.pen += 0.5 * tf.reduce_sum(tf.pow(self.seed_score, 2)) * self.rel_seed_l2_lambda
            self.pen += 0.5 * tf.reduce_sum(tf.pow(self.gene_score_avgs, 2)) * self.gene_l2_lambda
            self.pen += 0.5 * tf.reduce_sum(tf.pow(self.seed_score_avgs, 2)) * self.seed_l2_lambda

            #get total loss as likelihood plus penalty
            self.loss = self.nLL + self.pen 
            self.avg_effect_loss += self.pen

            #make optimizer op for score estimation
            score_var_list = [self.gene_score, self.gene_score_avgs, self.seed_score, self.seed_score_avgs,
                        self.hairpin_offset, self.hairpin_unpred, self.CL_offset]
            self.score_optim = ScipyOptimizerInterface(self.loss, 
                                            options = self.optim_params, 
                                            var_list = score_var_list,
                                            method = 'L-BFGS-B')

            #make optimizer op for estimating guide efficacies
            guide_var_list = [self.hairpin_offset, self.hairpin_unpred, self.CL_offset, self.guide_Geff, self.guide_Seff]
            n_uncon = self.n_hp_batches + self.n_hairpins + self.n_CL_batches
            n_bcon = self.n_Geffs + self.n_Seffs
            bound_constraints = np.concatenate([
                    np.tile([None, None], [n_uncon, 1]),
                    np.tile([0, 1], [n_bcon, 1])], 
                axis = 0)
            self.guide_optim = ScipyOptimizerInterface(self.loss, 
                                            options = self.optim_params, 
                                            var_list = guide_var_list,
                                            method = 'L-BFGS-B',
                                            bounds = bound_constraints)


            #make optimizer ops for estimating gene slopes
            gene_slope_var_list = [self.hairpin_offset, self.hairpin_unpred, self.CL_offset, self.gene_slope]
            n_uncon = self.n_hp_batches + self.n_hairpins + self.n_CL_batches
            n_pcon = self.n_CLs
            bound_constraints = np.concatenate([
                    np.tile([None, None], [n_uncon, 1]),
                    np.tile([0, None], [n_pcon, 1])], 
                axis = 0)
            self.gene_slope_optim = ScipyOptimizerInterface(self.avg_effect_loss, 
                                            options = self.optim_params, 
                                            var_list = gene_slope_var_list,
                                            method = 'L-BFGS-B',
                                            bounds = bound_constraints)

            #make optimizer ops for estimating CL slopes
            ov_slope_var_list = [self.hairpin_offset, self.CL_offset, self.CL_slope]
            n_uncon = self.n_hp_batches + self.n_CL_batches
            n_pcon = self.n_CL_batches
            bound_constraints = np.concatenate([
                    np.tile([None, None], [n_uncon, 1]),
                    np.tile([0, None], [n_pcon, 1])], 
                axis = 0)
            self.ov_slope_optim = ScipyOptimizerInterface(self.avg_effect_loss, 
                                            options = self.optim_params, 
                                            var_list = ov_slope_var_list,
                                            method = 'L-BFGS-B',
                                            bounds = bound_constraints)

            init = tf.global_variables_initializer()  

        self.sess.run(init)
        if self._log_file is not None:
            self._log_file.close()


    def write(self, data, silent = False):
        '''Internal method to print to stdout and logfile
        INPUTS:
            data: string to print_freq
            silent: print to terminal?
        '''
        if not silent:
            print(data)
        if self._log_file is not None:
            if self._log_file.closed:
                self._log_file = open(self.log_file, 'a')
            self._log_file.write(data + '\n')


    def predict(self):
        '''Get model prediction'''
        return(self.sess.run(self.pred))


    def get_SSMD(self, use_test = False):
        '''Calculate SSMD for a set of gene avgs, given sets of positive and negative controls'''
        gene_score, gene_score_avgs, CL_noise_vars = \
            self.sess.run([self.gene_score, self.gene_score_avgs, self.CL_noise_vars])
        gene_scores = gene_score + gene_score_avgs
        noise_vars_per_CL = pd.DataFrame({'noise_vars': CL_noise_vars.flatten(), 'CL_name': self.all_CL_names}) \
            .groupby('CL_name').mean().ix[self.data_names['CLs'],:].values
        weights = noise_vars_per_CL.reshape(-1,1) / np.nanmean(noise_vars_per_CL)
        weight_avg = np.nanmean(gene_scores * weights, axis = 0)
        if use_test: #if using cross-val on set of control genes
            pop1 = weight_avg[np.in1d(self.data_names['genes'], self.gene_sets['neg_test'])]
            pop2 = weight_avg[np.in1d(self.data_names['genes'], self.gene_sets['pos_test'])]
        else:
            pop1 = weight_avg[np.in1d(self.data_names['genes'], self.gene_sets['neg'])]
            pop2 = weight_avg[np.in1d(self.data_names['genes'], self.gene_sets['pos'])]
        return(SSMD(pop1, pop2))



    def fit(self, LFC_mats, fit_params = 'scores', ignore_test = False):
        '''
        Train subset of model parameters
        INPUTS:
            LFC_mats: List of [n_hairpins, n_CLs] training data sets of measured hairpin-level LFCs
            fit_params: model parameter set to estimate  ['scores', 'guide_effs', 'gene_slopes', 'ov_slopes', 'noise_vars', 'gene_slopes_ML', 'ov_slopes_ML']
            ignore_test: optional fit to all data even if test_inds are defined 
        '''
        poss_fit_params = ['scores', 'guide_effs', 'gene_slopes', 'ov_slopes', 'noise_vars', 'gene_slopes_ML', 'ov_slopes_ML']
        assert(fit_params in poss_fit_params)
        
        if self.log_file is not None:
            self._log_file = open(self.log_file, 'a')

        if self.test_inds is not None and not ignore_test:
            train_eval_masks = make_eval_masks(LFC_mats, self.test_inds)
            test_eval_masks = make_eval_masks(LFC_mats, self.test_inds, inverse = True)
        else:
            train_eval_masks = make_eval_masks(LFC_mats, None)

        LFC_mats_no_na = []
        for LFC_mat in LFC_mats:
            cur = LFC_mat.copy()
            cur[np.isnan(cur)] = 0
            LFC_mats_no_na.append(cur)

        feed_dict = {i: d for i, d in zip(self.obs, LFC_mats_no_na)}
        train_mask_dict = {i: d for i, d in zip(self.eval_mask, train_eval_masks)}
        train_dict = merge_dicts(feed_dict, train_mask_dict)
        if self.test_inds is not None and not ignore_test:
            test_mask_dict = {i: d for i, d in zip(self.eval_mask, test_eval_masks)}
            test_dict = merge_dicts(feed_dict, test_mask_dict)

        if fit_params == 'scores':
            R2_evals = self.sess.run(self.R2, feed_dict = train_dict)
            self.write('Init R2: {}'.format(R2_evals))
            if self.test_inds and not ignore_test:
                R2_evals = self.sess.run(self.R2, feed_dict = test_dict)
                self.write('Init Test R2: {}'.format(R2_evals))

        t0 = time.time()
        if fit_params == 'scores':
            self.write('Fitting model scores')
            optim_res = self._fit_scores(train_dict)
        elif fit_params == 'guide_effs':    
            self.write('Fitting guide efficacies')
            optim_res = self._fit_guide_efficacies(train_dict)
        elif fit_params == 'gene_slopes':
            self._fit_gene_slopes()
        elif fit_params == 'ov_slopes':
            self._fit_ov_slopes(LFC_data, ignore_test = ignore_test)
        elif fit_params == 'gene_slopes_ML':
            optim_res = self._fit_gene_slopes_ML(train_dict)
        elif fit_params == 'ov_slopes_ML':
            optim_res = self._fit_ov_slopes_ML(train_dict)
        elif fit_params == 'noise_vars':
            self._fit_noise_vars(LFC_mats, ignore_test = ignore_test)
        elif fit_params == 'slopes':
            self._fit_slopes(train_dict)
        if fit_params in ['scores', 'guide_effs', 'gene_slopes_ML', 'ov_slopes_ML']:
            self.write(optim_res['message'].decode('utf-8'))
            self.write('Optimization finished after: {} sec, {} iter, {} fevals'.format(int(time.time() - t0), 
                                                                            optim_res['nit'],optim_res['nfev']))           
        if fit_params == 'scores':
            R2_evals = self.sess.run(self.R2, feed_dict = train_dict)
            self.R2_vals['train'].append(R2_evals)
            self.write('New R2: {}'.format(R2_evals))
            if self.test_inds and not ignore_test:
                R2_evals = self.sess.run(self.R2, feed_dict = test_dict)
                self.R2_vals['test'].append(R2_evals)
                self.write('New Test R2: {}'.format(R2_evals))

            self.SSMD['train'].append(self.get_SSMD(use_test = False))
            self.write('Train SSMD: {}'.format(self.SSMD['train'][-1]))
            if ('pos_test' in self.gene_sets) and (len(self.gene_sets['pos_test']) > 0):
                self.SSMD['test'].append(self.get_SSMD(use_test = True))
                self.write('Test SSMD: {}'.format(self.SSMD['test'][-1]))
     
        self.loss_evals.append(self.sess.run(self.loss, feed_dict = train_dict))
        if self._log_file is not None:
            self._log_file.close()

 
    def _fit_scores(self, feed_dict):
        '''
        Fit scores + intercepts using BFGS 
        INPUTS:
            feed_dict: input data
            optim_params: dict of optimization parameters
        '''        
        optim_res = run_sklearn_optim(self.score_optim, feed_dict, self.sess, self.loss, 
            print_freq = self.optim_params['print_freq'])
        
        return(optim_res)


    def _fit_gene_slopes_ML(self, feed_dict):
        '''
        Fit slopes + intercepts using BFGS 
        INPUTS:
            feed_dict: input data
            optim_params: dict of optimization parameters
        '''        
        init_gene_slopes = self.sess.run([self.gene_slope])
        optim_res = run_sklearn_optim(self.gene_slope_optim, feed_dict, self.sess, self.avg_effect_loss, 
            print_freq = self.optim_params['print_freq'])
        new_gene_slopes = self.sess.run(self.gene_slope)
        self.write('init gene slopes avg: {}, new gene slope avg: {}'.format(np.mean(init_gene_slopes), np.mean(new_gene_slopes)))               
        new_gene_slopes[new_gene_slopes < self.min_slope] = self.min_slope #constrain to be non negative
        # new_gene_slopes = euclidean_proj_simplex(new_gene_slopes.flatten(), s=self.n_CLs).reshape(1,-1)
        new_gene_slopes = new_gene_slopes / np.nanmean(new_gene_slopes)
        _=self.sess.run(self.gene_slope.assign(new_gene_slopes.reshape(-1,1)))
        return(optim_res)


    def _fit_ov_slopes_ML(self, feed_dict):
        '''
        Fit slopes + intercepts using BFGS 
        INPUTS:
            feed_dict: input data
            optim_params: dict of optimization parameters
        '''        
        init_CL_slopes = self.sess.run([self.CL_slope])
        optim_res = run_sklearn_optim(self.ov_slope_optim, feed_dict, self.sess, self.avg_effect_loss, 
            print_freq = self.optim_params['print_freq'])
        new_CL_slopes = self.sess.run(self.CL_slope)
        self.write('init ov slopes avg: {}, new ov slope avg: {}'.format(np.mean(init_CL_slopes), np.mean(new_CL_slopes)))               
        new_CL_slopes[new_CL_slopes < self.min_slope] = self.min_slope
        # new_CL_slopes = euclidean_proj_simplex(new_CL_slopes.flatten(), s=self.n_CLs).reshape(1,-1)
        new_CL_slopes = new_CL_slopes / np.nanmean(new_CL_slopes)
        _=self.sess.run(self.CL_slope.assign(new_CL_slopes))
        return(optim_res)


    def _fit_gene_slopes(self):
        '''Re-estimate gene score slope terms using pos/neg control gene set median separation'''        
        init_gene_slopes, init_gene_effects = self.sess.run([self.gene_slope, self.ind_gene_effects])

        # NA out gene scores for cell lines where we dont have targeting guides
        init_gene_effects[self.n_used_hairpins_per_gene.transpose() < self.min_hairpins_per] = np.nan

        #estimate centers of positive and negative gene set distributions
        pos_med = np.nanmedian(init_gene_effects[:, np.in1d(self.data_names['genes'], self.gene_sets['pos'])], axis = 1)
        neg_med = np.nanmedian(init_gene_effects[:, np.in1d(self.data_names['genes'], self.gene_sets['neg'])], axis = 1)
        new_gene_slopes = neg_med - pos_med
        self.write('negative gene slopes: {}/{}'.format(np.sum(new_gene_slopes < 0), self.n_CLs))
        self.write('init gene slopes avg: {}, new gene slope avg: {}'.format(np.mean(init_gene_slopes), np.mean(new_gene_slopes)))       
        new_gene_slopes = new_gene_slopes / np.nanmean(new_gene_slopes) #normalize to have mean 1
        _=self.sess.run(self.gene_slope.assign(new_gene_slopes.reshape(-1,1)))


    def _fit_guide_efficacies(self, feed_dict):
        '''
        Fit guide_efficacies + intercepts using BFGS 
        INPUTS:
            feed_dict: input data
            optim_params: dict of optimization parameters
        '''
        init_guide_Geffs, init_guide_Seffs = self.sess.run([self.guide_Geff, self.guide_Seff])
        optim_res = run_sklearn_optim(self.guide_optim, feed_dict, self.sess, self.loss, 
            print_freq = self.optim_params['print_freq']) 
        new_guide_Geffs, new_guide_Seffs = self.sess.run([self.guide_Geff, self.guide_Seff])
        self.write('init avg Geff: {} Seff: {}, new avg Geff: {} Seff: {}'.format(np.mean(init_guide_Geffs), 
            np.mean(init_guide_Seffs), np.mean(new_guide_Geffs), np.mean(new_guide_Seffs))) 
        return(optim_res)


    def _fit_noise_vars(self, LFC_mats, ignore_test = False):
        '''Estimate noise variance per CL'''

        tot_SSE = np.zeros(self.n_CL_batches)
        tot_used_hps = np.zeros(self.n_CL_batches)
        batch_offset = 0
        for batch_ii, (LFC_mat, pred_mat) in enumerate(zip(LFC_mats, self.predict())):
            cur_CL_inds = np.arange(LFC_mat.shape[1]) + batch_offset
            cur_d = LFC_mat.values.copy()
            if not ignore_test and self.test_inds is not None:
                cur_d[self.test_inds[batch_ii]] = np.nan
            tot_SSE[cur_CL_inds] += np.nansum((pred_mat - cur_d)**2, axis = 0)
            tot_used_hps[cur_CL_inds] += np.sum(~np.isnan(cur_d), axis = 0)
            batch_offset += LFC_mat.shape[1]

        # dof = tot_used_hps - self.n_targeted_genes - self.n_targeted_seeds - 1 #dof per CL (approximate)
        dof = tot_used_hps #dof per CL (gives biased estimate)
        per_CL_noise_var = tot_SSE / np.max(np.concatenate([dof.reshape(-1,1), np.ones((self.n_CL_batches, 1))], axis = 1), axis = 1)
        self.sess.run(self.CL_noise_vars.assign(per_CL_noise_var.reshape(-1,1).astype(np.float32)))
 

    def compute_R2_stats(self, LFC_mats):
        '''
        Computes R2 values per CL, and per hairpin 
        '''
        self.CL_R2_df = pd.DataFrame()
        self.hp_R2_df = pd.DataFrame()
        for batch_id, (LFC_data, pred) in enumerate(zip(LFC_mats, self.predict())):
            resids = LFC_data - pred 
            CL_noise_var = np.nanvar(resids, axis = 0)
            CL_tot_var = np.nanvar(LFC_data, axis = 0)
            CL_R2 = 1 - CL_noise_var / CL_tot_var
            self.CL_R2_df = pd.concat([self.CL_R2_df, 
                                  pd.DataFrame({'CCLE_ID': LFC_data.columns.values, 
                                                'batch_id': np.ones_like(CL_R2)*batch_id,
                                                'R2': CL_R2})])

            hp_noise_var = np.nanvar(resids, axis = 1)
            hp_tot_var = np.nanvar(LFC_data, axis = 1)
            hp_R2 = 1 - hp_noise_var / hp_tot_var
            self.hp_R2_df = pd.concat([self.hp_R2_df, 
                                  pd.DataFrame({'hp_seq': LFC_data.index.values, 
                                                'batch_id': np.ones_like(hp_R2)*batch_id,
                                                'R2': hp_R2})])


 

    def init_slopes(self, LFC_mats):
        '''Get initial estimates of CL slopes by regressing each CL's data on within-batch avg'''
        lm = LinearRegression(fit_intercept=True)

        #first get overall slope adjustment per batch
        if len(LFC_mats) > 0:
            common_hps = reduce(np.intersect1d, [LFC_mat.index.values for LFC_mat in LFC_mats])
        else:
            common_hps = LFC_mats[0].index.values
        if len(common_hps) > 100:
            per_batch_avgs = np.ones((len(LFC_mats), len(common_hps)))
            for ii, LFC_mat in enumerate(LFC_mats):
                cur_d = LFC_mat.ix[np.in1d(LFC_mat.index.values, common_hps),:]
                cur_d = cur_d.reset_index().drop_duplicates(subset='index', keep='last').set_index('index')
                per_batch_avgs[ii,:] = np.nanmean(cur_d.ix[common_hps,:].values, axis = 1)
 
            ov_avg = np.nanmean(per_batch_avgs, axis = 0)
            batch_slopes = np.ones(per_batch_avgs.shape[0])
            for ii in range(per_batch_avgs.shape[0]):
                uset = np.where(~np.isnan(per_batch_avgs[ii,:]))[0]
                lm.fit(ov_avg.reshape(-1,1)[uset,:], per_batch_avgs[ii,uset].transpose())
                batch_slopes[ii] = lm.coef_
        else:
            batch_slopes = np.array([np.nanstd(LFC_mat.values) for LFC_mat in LFC_mats])
            batch_slopes = batch_slopes / np.nanmean(batch_slopes)

        CL_slopes = np.ones(self.n_CL_batches)
        CL_offset = 0
        for batch_ii, LFC_mat in enumerate(LFC_mats):
            avg_hp = np.nanmean(LFC_mat, axis = 1)
            for ii in np.arange(LFC_mat.shape[1]):
                uvals = np.where(~np.isnan(LFC_mat.values[:,ii]))[0]
                lm.fit(avg_hp.reshape(-1,1)[uvals,:], LFC_mat.values[uvals,ii])
                CL_slopes[ii + CL_offset] = lm.coef_ * batch_slopes[batch_ii]
            CL_offset += LFC_mat.shape[1]
                
        _=self.sess.run(self.CL_slope.assign(CL_slopes.reshape(-1,1)))


    def save(self, results_dir, save_perf_only = False, edward = False):
        '''
        Write parameter matrices to text files. Also serialize other model params to json file at specified path
        '''
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)

        if not edward:
            other_df = {}
            other_df['reg_params'] = self.reg_params
            other_df['R2_vals'] = self.R2_vals
            other_df['optim_params'] = self.optim_params
            other_df['loss_evals'] = self.loss_evals
            other_df['SSMD'] = self.SSMD
            with open(os.path.join(results_dir, 'other_info.json'), 'w') as f:
                dump(other_df, f, primitives = True, allow_nan = True)
 
 
        if not save_perf_only: #if not just saving performance params
            CL_df = pd.DataFrame({
                'CCLE_ID': self.data_names['CLs'],
                'gene_slope': self.gene_slope.eval(session = self.sess).flatten()
                # 'noise_vars': self.CL_noise_vars.eval(session = self.sess).flatten()
            })
            CL_df.to_csv(os.path.join(results_dir, 'CL_data.csv'), index = False)

            CL_batch_df = pd.DataFrame({'CCLE_ID': self.all_CL_names,
                           'CL_slope': self.CL_slope.eval(session = self.sess).flatten(),
                           'CL_offset': self.CL_offset.eval(session = self.sess).flatten(),
                           'CL_batch': self.all_CL_batches,
                           'noise_vars': self.CL_noise_vars.eval(session = self.sess).flatten()})
            if hasattr(self, 'CL_R2'):
                CL_batch_df['R2'] = self.CL_R2_df['R2']
            if edward:
                CL_batch_df['offset_mean'] = self.q_CL_offset.loc.eval().flatten()
                CL_batch_df['offset_sd'] = self.q_CL_offset.scale.eval().flatten()
            CL_batch_df.to_csv(os.path.join(results_dir, 'CL_batch_data.csv'), index = False)

            hp_df = pd.DataFrame({
                'hp': self.data_names['hps'],
                'unpred_offset': self.hairpin_unpred.eval(session = self.sess).flatten(),
                'Geff': self.guide_Geff.eval(session = self.sess).flatten(),
                'Seff': self.guide_Seff.eval(session = self.sess).flatten()
            })
            if edward:
                hp_df['unpred_offset_mean'] = self.q_hairpin_unpred.loc.eval()
                hp_df['unpred_offset_sd'] = self.q_hairpin_unpred.scale.eval()
            hp_df.to_csv(os.path.join(results_dir, 'hp_data.csv'), index = False)

            hp_batch_df = pd.DataFrame({
                'hp': self.all_hp_seqs,
                'hp_batch': self.all_hp_batches,
                'hairpin_offset': self.hairpin_offset.eval(session = self.sess).flatten()
                })                       
            if hasattr(self, 'hp_R2'):
                hp_batch_df['R2'] = self.hp_R2_df['R2']            
            if edward:
                hp_batch_df['hairpin_offset_mean'] = self.q_hairpin_offset.loc.eval()
                hp_batch_df['hairpin_offset_sd'] = self.q_hairpin_offset.scale.eval()
            hp_batch_df.to_csv(os.path.join(results_dir, 'hp_batch_data.csv'), index = False)

          
            per_gene_df = pd.DataFrame({'avg': self.gene_score_avgs.eval(session = self.sess).flatten()},
                index = self.data_names['genes'])
            if edward:
                gene_mean_df = pd.DataFrame((self.q_gene_score.loc.eval() + \
                    self.q_gene_score_avgs.loc.eval()).transpose(),
                    index = self.data_names['genes'], columns = self.data_names['CLs'])
                gene_sd_df = pd.DataFrame(np.sqrt(self.q_gene_score.scale.eval()**2 + \
                    self.q_gene_score_avgs.scale.eval()**2).transpose(),
                    index = self.data_names['genes'], columns = self.data_names['CLs'])
                gene_mean_df = gene_mean_df.where(self.n_used_hairpins_per_gene >= self.min_hairpins_per, other = np.nan)
                gene_sd_df = gene_sd_df.where(self.n_used_hairpins_per_gene >= self.min_hairpins_per, other = np.nan)
                gene_mean_df.to_csv(os.path.join(results_dir, 'gene_means.csv'))
                gene_sd_df.to_csv(os.path.join(results_dir, 'gene_SDs.csv'))
                per_gene_df['SD'] = self.q_gene_score_avgs.scale.eval().flatten()
            else:
                gene_df = pd.DataFrame((self.gene_score.eval(session = self.sess) + \
                                self.gene_score_avgs.eval(session = self.sess)).transpose(),
                              index = self.data_names['genes'], columns = self.data_names['CLs'])
                gene_df = gene_df.where(self.n_used_hairpins_per_gene >= self.min_hairpins_per, other = np.nan)
                gene_df.to_csv(os.path.join(results_dir, 'gene_data.csv'))
            per_gene_df.to_csv(os.path.join(results_dir, 'per_gene_data.csv'))

            if edward:
                seed_mean_df = pd.DataFrame((self.q_seed_score.loc.eval() + \
                    self.q_seed_score_avgs.loc.eval()).transpose(),
                    index = self.data_names['seeds'], columns = self.data_names['CLs'])
                seed_sd_df = pd.DataFrame(np.sqrt(self.q_seed_score.scale.eval()**2 + \
                    self.q_seed_score_avgs.scale.eval()**2).transpose(),
                    index = self.data_names['seeds'], columns = self.data_names['CLs'])
                seed_mean_df = seed_mean_df.where(self.n_used_hairpins_per_seed >= self.min_hairpins_per, other = np.nan)
                seed_sd_df = seed_sd_df.where(self.n_used_hairpins_per_seed >= self.min_hairpins_per, other = np.nan)
                seed_mean_df.to_csv(os.path.join(results_dir, 'seed_means.csv'))
                seed_sd_df.to_csv(os.path.join(results_dir, 'seed_SDs.csv'))
            else:
                seed_df = pd.DataFrame((self.seed_score.eval(session = self.sess) + \
                                self.seed_score_avgs.eval(session = self.sess)).transpose(),
                              index = self.data_names['seeds'], columns = self.data_names['CLs'])
                seed_df = seed_df.where(self.n_used_hairpins_per_seed >= self.min_hairpins_per, other = np.nan)
                seed_df.to_csv(os.path.join(results_dir, 'seed_data.csv'))

 
    def make_edward_model(self, LFC_mats, gene_matrix, seed_matrix, data_names):
        '''Create a Bayesian model in edward, using current parameter estimates to initialize'''

        #define priors on parameters
        gene_score = Normal(loc=tf.zeros([self.n_CLs, self.n_genes], dtype = self.float), 
            scale = np.sqrt(1.0/self.reg_params['rel_gene_l2_lambda']) * tf.ones([self.n_CLs, self.n_genes], dtype = self.float))
        seed_score = Normal(loc=tf.zeros([self.n_CLs, self.n_seeds], dtype = self.float), 
            scale = np.sqrt(1.0/self.reg_params['rel_seed_l2_lambda']) * tf.ones([self.n_CLs, self.n_seeds], dtype = self.float))
        gene_score_avgs = Normal(loc=tf.zeros([1, self.n_genes], dtype = self.float), 
            scale = np.sqrt(1.0/self.reg_params['gene_l2_lambda']) * tf.ones([1, self.n_genes], dtype = self.float))
        seed_score_avgs = Normal(loc=tf.zeros([1, self.n_seeds], dtype = self.float), 
            scale = np.sqrt(1.0/self.reg_params['seed_l2_lambda']) * tf.ones([1, self.n_seeds], dtype = self.float))
        hairpin_offset = Normal(loc=tf.zeros([self.n_hp_batches, 1], dtype = self.float), 
            scale = np.sqrt(1.0/self.reg_params['hairpin_l2_lambda']) * tf.ones([self.n_hp_batches, 1], dtype = self.float))
        hairpin_unpred = Normal(loc=tf.zeros([self.n_hairpins, 1], dtype = self.float), 
            scale = np.sqrt(1.0/self.reg_params['hp_unpred_l2_lambda']) * tf.ones([self.n_hairpins, 1], dtype = self.float))
        CL_offset = Normal(loc=tf.zeros([self.n_CL_batches, 1], dtype = self.float),
            scale = np.sqrt(1.0/self.reg_params['CL_l2_lambda']) * tf.ones([self.n_CL_batches, 1], dtype = self.float))

        #parameters we dont try to fit here
        CL_slope = tf.constant(self.CL_slope.eval(session = self.sess))
        gene_slope = tf.constant(self.gene_slope.eval(session = self.sess))
        # region_weights = tf.constant(self.region_weights.eval(session = self.sess))
        guide_Geff = tf.constant(self.guide_Geff.eval(session = self.sess))
        guide_Seff = tf.constant(self.guide_Seff.eval(session = self.sess))

        gene_maps = [self.make_sparse_submap(gene_matrix, LFC_mat.index.values) for LFC_mat in LFC_mats]
        seed_maps = [self.make_sparse_submap(seed_matrix, LFC_mat.index.values) for LFC_mat in LFC_mats]

        self.noise_sigma = tf.sqrt(tf.exp(tf.Variable(np.log(self.CL_noise_vars.eval(session = self.sess)), dtype = self.float)))

        mod_params = {
            'CL_offset': CL_offset,
            'CL_slope': CL_slope,
            'guide_Seff': guide_Seff,
            'guide_Geff': guide_Geff,
            'hairpin_offset': hairpin_offset,
            'hairpin_unpred': hairpin_unpred,
            'gene_slope': gene_slope,
            'seed_score_avgs': seed_score_avgs,
            'seed_score': seed_score,
            'gene_score_avgs': gene_score_avgs,
            'gene_score': gene_score
            }

        bool_masks = [tf.logical_not(tf.is_nan(LFC_mat.values)) for LFC_mat in LFC_mats]
        y_list = []
        CL_cnt = 0
        hp_cnt = 0
        for ii in range(len(self.obs)):
            cur_pred = self.get_dataset_pred(
                mod_params,
                gene_maps[ii],
                seed_maps[ii],
                LFC_mats[ii].index.values,
                LFC_mats[ii].columns.values,
                hp_cnt, 
                CL_cnt)
            cur_CL_inds = np.arange(LFC_mats[ii].shape[1]) + CL_cnt
            hp_cnt += LFC_mats[ii].shape[0]
            CL_cnt += LFC_mats[ii].shape[1]
            cur_sigma = tf.transpose(tf.gather(self.noise_sigma, cur_CL_inds)) * tf.ones_like(cur_pred)
            y_list.append(Normal(loc=tf.boolean_mask(cur_pred, bool_masks[ii]), 
                scale = tf.boolean_mask(cur_sigma, bool_masks[ii])))

        LFC_mats_no_na = []
        for LFC_mat in LFC_mats:
            cur = LFC_mat.values.copy()
            cur[np.isnan(cur)] = 0
            LFC_mats_no_na.append(cur)

        # obs_list = [tf.constant(LFC_mat, dtype = 'float') for LFC_mat in LFC_mats_no_na]
        obs_list = [tf.placeholder(self.float, dset.shape) for dset in LFC_mats_no_na]

        #posterior approximating distributions (fully factorized gaussian)
        self.q_gene_score = Normal(loc=tf.Variable(self.gene_score.eval(session = self.sess), dtype = self.float), 
                        scale=tf.nn.softplus(tf.Variable(tf.random_normal([self.n_CLs, self.n_genes], dtype = self.float))))
        self.q_seed_score = Normal(loc=tf.Variable(self.seed_score.eval(session = self.sess), dtype = self.float), 
                        scale=tf.nn.softplus(tf.Variable(tf.random_normal([self.n_CLs, self.n_seeds], dtype = self.float))))
        self.q_gene_score_avgs = Normal(loc=tf.Variable(self.gene_score_avgs.eval(session = self.sess), dtype = self.float), 
                        scale=tf.nn.softplus(tf.Variable(tf.random_normal([1, self.n_genes], dtype = self.float))))
        self.q_seed_score_avgs = Normal(loc=tf.Variable(self.seed_score_avgs.eval(session = self.sess), dtype = self.float), 
                        scale=tf.nn.softplus(tf.Variable(tf.random_normal([1, self.n_seeds], dtype = self.float))))
        self.q_hairpin_offset = Normal(loc=tf.Variable(self.hairpin_offset.eval(session = self.sess), dtype = self.float), 
                        scale=tf.nn.softplus(tf.Variable(tf.random_normal([self.n_hp_batches, 1], dtype = self.float))))
        self.q_hairpin_unpred = Normal(loc=tf.Variable(self.hairpin_unpred.eval(session = self.sess), dtype = self.float), 
                        scale=tf.nn.softplus(tf.Variable(tf.random_normal([self.n_hairpins, 1], dtype = self.float))))
        self.q_CL_offset = Normal(loc=tf.Variable(self.CL_offset.eval(session = self.sess), dtype = self.float), 
                        scale=tf.nn.softplus(tf.Variable(tf.random_normal([self.n_CL_batches, 1], dtype = self.float))))

        data_dict = {i: tf.boolean_mask(d, m) for i, d, m in zip(y_list, obs_list, bool_masks)}
        for i, d in zip(obs_list, LFC_mats):
            data_dict.update({i: d})
            
        self.inference = ed.KLqp({gene_score: self.q_gene_score, 
                                seed_score: self.q_seed_score,
                                gene_score_avgs: self.q_gene_score_avgs,
                                seed_score_avgs: self.q_seed_score_avgs,
                                hairpin_offset: self.q_hairpin_offset,
                                hairpin_unpred: self.q_hairpin_unpred,
                                CL_offset: self.q_CL_offset}, 
                            data=data_dict)
        self.inference.initialize()
        tf.global_variables_initializer().run()


    def run_edward_inference(self, n_iter = 1000, print_freq = 100):      
        loss_evals = np.zeros(n_iter)
        orig_GS = self.gene_score.eval(session = self.sess).flatten()
        tot_GS_var = np.var(orig_GS)
        delta_G_R2 = 1 - np.var(orig_GS - self.q_gene_score.mean().eval().flatten()) / tot_GS_var
        self.write('Init DeltaG_R2: {}'.format(delta_G_R2))  
        for ii in range(n_iter):
            info_dict = self.inference.update()
            loss_evals[ii] = info_dict['loss']
            if ii % print_freq == 0:
                delta_G_R2 = 1 - np.var(orig_GS - self.q_gene_score.mean().eval().flatten()) / tot_GS_var
                self.write('It: {}, DeltaG_R2: {}, Loss: {}'.format(ii, delta_G_R2, loss_evals[ii]))  
        return(loss_evals)  


    ## prediction
    def get_dataset_pred(self, mod_params, gene_map, seed_map, cur_hp_seqs, cur_CL_names, hp_offset, CL_offset, just_avg_scores = False):
        cur_hp_inds = np.array([self.hp_ind_map[x] for x in cur_hp_seqs])
        cur_CL_inds = np.array([self.CL_ind_map[x] for x in cur_CL_names])
        CL_batch_range = CL_offset + np.arange(len(cur_CL_names))
        hp_batch_range = hp_offset + np.arange(len(cur_hp_seqs))

        if just_avg_scores:
            cur_gene_effect = map_effects(mod_params['gene_score_avgs'] + tf.zeros_like(mod_params['gene_score'], dtype = self.float), 
                cur_CL_inds, gene_map)
        else:
            cur_gene_effect = map_effects(mod_params['gene_score'] + mod_params['gene_score_avgs'], 
                cur_CL_inds, gene_map)
        cur_gene_effect = tf.multiply(cur_gene_effect, 
            tf.gather(mod_params['guide_Geff'],cur_hp_inds))
        cur_gene_effect = tf.multiply(cur_gene_effect, 
            tf.reshape(tf.gather(mod_params['gene_slope'], cur_CL_inds), [1, -1])) 

        if just_avg_scores:
            cur_seed_effect = map_effects(mod_params['seed_score_avgs'] + tf.zeros_like(mod_params['seed_score'], dtype = self.float), 
                cur_CL_inds, seed_map)
        else:
            cur_seed_effect = map_effects(mod_params['seed_score'] + mod_params['seed_score_avgs'], 
                cur_CL_inds, seed_map)
        cur_seed_effect = tf.multiply(cur_seed_effect, 
            tf.gather(mod_params['guide_Seff'], cur_hp_inds))
        
        #total KD effect of each hp
        cur_KD_effect = tf.gather(mod_params['hairpin_unpred'], cur_hp_inds) + cur_gene_effect + cur_seed_effect 
        cur_KD_effect = tf.multiply(cur_KD_effect, 
            tf.reshape(tf.gather(mod_params['CL_slope'], CL_batch_range), [1, -1]))
       
        cur_pred = cur_KD_effect + tf.reshape(tf.gather(mod_params['CL_offset'], CL_batch_range), [1, -1]) + \
                tf.gather(mod_params['hairpin_offset'], hp_batch_range) #add offset terms

        return(cur_pred)


    def get_dataset_LL(self, mod_params, LFC_mat, preds, cur_eval_mask, CL_offset):
        CL_batch_range = CL_offset + np.arange(LFC_mat.get_shape().as_list()[1])
        cur_nLL = 0.5 * tf.reduce_sum(
            tf.boolean_mask(
                tf.multiply(tf.pow(preds - LFC_mat, 2), 
                    1/tf.reshape(tf.gather(mod_params['CL_noise_vars'], CL_batch_range), [1, -1])),
                cur_eval_mask
            ))
        
        cur_SS = 0.5 * tf.reduce_sum(
            tf.boolean_mask(
                tf.multiply(tf.pow(LFC_mat, 2), 
                    1/tf.reshape(tf.gather(mod_params['CL_noise_vars'], CL_batch_range), [1, -1])),
                cur_eval_mask
            ))

        return(cur_nLL, cur_SS)


    def get_shRNA_R2(self, mod_params, LFC_mat, preds, CL_offset):
        CL_batch_range = CL_offset + np.arange(LFC_mat.get_shape().as_list()[1])
        preds_ms = preds - tf.reduce_mean(preds, axis = 1, keep_dims = True)
        LFC_mat_ms = LFC_mat - tf.reduce_mean(LFC_mat, axis = 1, keep_dims = True)
        cur_nLL = 0.5 * tf.reduce_sum(
                tf.multiply(tf.pow(preds_ms - LFC_mat_ms, 2), 
                    1/tf.reshape(tf.gather(mod_params['CL_noise_vars'], CL_batch_range), [1, -1])),
                axis = 1)
        
        cur_SS = 0.5 * tf.reduce_sum(
                tf.multiply(tf.pow(LFC_mat_ms, 2), 
                    1/tf.reshape(tf.gather(mod_params['CL_noise_vars'], CL_batch_range), [1, -1])),
                axis = 1)
        cur_R2 = 1 - tf.div(cur_nLL , cur_SS)
        return(cur_R2, cur_nLL, cur_SS)



    def make_sparse_submap(self, sparse_hp_mat, cur_hp_seqs):
        '''Extract a set of rows for specific hairpins from a sparse matrix'''
        cur_hp_inds = np.array([self.hp_ind_map[x] for x in cur_hp_seqs])
        map_coo = sparse_hp_mat[cur_hp_inds,:].tocoo()
        map_inds = np.concatenate([map_coo.row[:,None], map_coo.col[:,None]], axis = 1)
        if self.float == tf.float32:
            vals = tf.to_float(map_coo.data)
        else:
            vals = tf.to_double(map_coo.data)
        return(tf.SparseTensor(indices=map_inds, values=vals, dense_shape=map_coo.shape))
