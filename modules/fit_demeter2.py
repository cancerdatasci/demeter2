import sys
import os
import argparse
import time
import pickle
import re
import numpy as np
import pandas as pd

from sklearn.linear_model import LinearRegression
import tensorflow as tf

if "../modules" not in sys.path:
    sys.path.append("../modules")
import preprocess
from taigapy import TaigaClient

import demeter2

import warnings
warnings.filterwarnings('ignore')

tf.logging.set_verbosity(tf.logging.INFO)

alt_tol = 1e-5 #convergence relative tolerance (on change in cost) for blockwise coord descent procedure
fit_CL_slopes = True 
regions = [(10,17),(11,18)] #seed sequence regions
data_dir = './data/'

parser = argparse.ArgumentParser()
parser.add_argument("--model_dir", help="directory name to save model result", default = None)
parser.add_argument("--data_files", help="file containing paths for each data file", default = None)
parser.add_argument("--sh_targets", help="file containing shRNA-to-gene mapping", default = None)
parser.add_argument("--pos_cons", help="file containing positive control genes", default = None)
parser.add_argument("--neg_cons", help="file containing negative control genes", default = None)
parser.add_argument("--zscore", help="zscore data for each CL", action = 'store_true')
parser.add_argument("--mean_sub", help="mean-subtract LFC data per hairpins", action = 'store_true')
parser.add_argument("--max_genes", help="max number of genes to use", type = int, default = np.Inf)
parser.add_argument("--max_CLs", help="max number of CLs to use", type = int, default = np.Inf)
parser.add_argument("--xv_frac", help="fraction of data to use for test R2 eval", type = float, default = 0.0)
parser.add_argument("--gene_xv_frac", help="fraction of pos and neg control genes to use for test SSMD eval", type = float, default = 0.0)
parser.add_argument("--no_edward", dest='fit_edward', help="also fit bayesian model", action = 'store_false')
parser.add_argument("--fit_noise", dest='fit_noise', help="alternating estimation of noise variances", action = 'store_true')
parser.add_argument("--save_perf_only", dest='save_perf_only', help="only save model performance measures", action = 'store_true')
parser.add_argument("--max_alt_iters", help="Maximum number of alternating optimization iterations", type = int, default = 20)
parser.add_argument("--print_freq", help="Print progress per n iterations", type = int, default = 50)
parser.add_argument("--random_seed", help="random seed", type = int, default = 0)
parser.add_argument("--ML_gene_slope", dest = 'ML_gene_slope', help="use ML estimate for gene slopes", action = 'store_true')
parser.add_argument("--float_precision", dest = 'float_precision', help="set float precision ('single' or 'double')", default = 'double')
parser.add_argument("--no_shRNA_effs", dest = 'fit_shRNA_effs', help="dont fit shRNA efficacy terms", action = 'store_false')

parser.set_defaults(zscore=False)
parser.set_defaults(fit_edward=True)
parser.set_defaults(fit_noise=False)
parser.set_defaults(mean_sub=False)
parser.set_defaults(save_perf_only=False)
parser.set_defaults(ML_gene_slope=False)

#regularization parameters
parser.add_argument("--gene_l2_lambda", type = float, 
	help="l2 penalty on (avg) gene scores", default = 1.0) 
parser.add_argument("--seed_l2_lambda", type = float, 
	help="l2 penalty on (avg) seed scores", default = 1.0) 
parser.add_argument("--rel_gene_l2_lambda", type = float, 
    help="l2 penalty on relative gene scores", default = 1.0) 
parser.add_argument("--rel_seed_l2_lambda", type = float, 
	help="l2 penalty on relative seed scores", default = 2.0) 
parser.add_argument("--hairpin_l2_lambda", type = float, 
    help="l2 penalty on per-hp offset", default = 10.0) 
parser.add_argument("--hp_unpred_l2_lambda", type = float, 
    help="l2 penalty on per-hp offset (unpredicted seed effects)", default = 10.0) 
parser.add_argument("--CL_l2_lambda", type = float, 
    help="l2 penalty on per-CL offset", default = 0.001)

#optimization parameters
parser.add_argument("--maxiter", type = int, 
	help="Maximum number of iterations for optimization", default = 2000)
parser.add_argument("--edward_iters", type = int, 
    help="Number of edward iterations", default = 10000)

args = parser.parse_args()

t0 = time.time()
def total_elapsed():
    print('Time elapsed: {}'.format(int(time.time() - t0)))

if args.model_dir is None:
    model_path = './'
else:
    model_path = args.model_dir
if not os.path.exists(model_path):
    os.makedirs(model_path)
log_file = os.path.join(model_path, 'logfile')
print('Logging at ' + log_file)

np.random.seed(args.random_seed) #set RNG seed

if os.path.exists('/taiga/token'):
    tc = TaigaClient(token_path = '/taiga/token')
elif os.path.exists('taiga/token'):
    tc = TaigaClient(token_path = 'taiga/token')

input_data = pd.read_csv(args.data_files, sep = '\t').dropna(axis = 0, how = 'all')
def get_dset(input_info):
    if len(input_info) == 1:
        print('Loading dataset from: {}'.format(input_info['data_file_paths']))
        cur_data = pd.read_csv(input_info['data_file_paths'], index_col = 0)
    else:
        print('Fetching taiga dataset: {} file: {}, version: {}'.format(input_info['data_set'], 
            input_info['data_file'], 
            int(input_info['version'])))
        if pd.isnull(input_info['data_file']):
            cur_data = tc.get(
                name = input_info['data_set'],
                version = int(input_info['version']))
        else:
            cur_data = tc.get(
                name = input_info['data_set'],
                file = input_info['data_file'],
                version = int(input_info['version']))
    return(cur_data)
dsets = [get_dset(x) for ii, x in input_data.iterrows()]
print('Read {} data files'.format(len(dsets)))

if (args.sh_targets is not None) and os.path.exists(args.sh_targets):
    sh_targets = pd.read_csv(args.sh_targets)
else:
    sh_targets = tc.get(name='gpp-shrna-mapping-8759', version=6, file='shmap_19mer_noXLOC_Entrezonly')
sh_targets.rename(columns = {'Barcode Sequence': 'SEQ', 'Gene ID': 'Gene_ID'}, inplace=True)

sh_targets.dropna(subset = ['Gene_ID'], inplace=True)

#load curated pos and negative control gene sets
print('Loading positive/negative control sets')
#load Entrez IDs for pos and neg con genes
if (args.pos_cons is not None) and os.path.exists(args.pos_cons):
    pos_con_genes = pd.read_csv(args.pos_cons)['Gene_ID'].values
else:
    pos_con_genes = tc.get(name='demeter2-pos-neg-controls-a5c6', version=1, file='hart_pos_controls')['Gene_ID'].values
if (args.neg_cons is not None) and os.path.exists(args.neg_cons):
    neg_con_genes = pd.read_csv(args.neg_cons)['Gene_ID'].values
else:
    neg_con_genes = tc.get(name='demeter2-pos-neg-controls-a5c6', version=1, file='hart_neg_controls')['Gene_ID'].values

if (args.gene_xv_frac > 0):
    train_neg_con_genes = np.random.choice(neg_con_genes, np.round((1-args.gene_xv_frac) * len(neg_con_genes)), replace = False)
    train_pos_con_genes = np.random.choice(pos_con_genes, np.round((1-args.gene_xv_frac) * len(pos_con_genes)), replace = False)
else:
    train_neg_con_genes = neg_con_genes
    train_pos_con_genes = pos_con_genes
test_neg_con_genes = np.setdiff1d(neg_con_genes, train_neg_con_genes)
test_pos_con_genes = np.setdiff1d(pos_con_genes, train_pos_con_genes)
print('Using {} positive and {} negative control genes for training'.format(len(train_pos_con_genes), 
    len(train_neg_con_genes)))
if args.gene_xv_frac > 0:
    print('Using {} positive and {} negative control genes for testing'.format(len(test_pos_con_genes), 
        len(test_neg_con_genes)))

#parse data
print('Making processed data')
data = preprocess.make_demeter2_data(dsets, sh_targets)

#subsample data if desired
if not (np.isinf(args.max_genes) and np.isinf(args.max_CLs)):
    if not np.isinf(args.max_genes):
        n_gene_sample = np.min([args.max_genes, len(data['unique_genes'])])
        used_genes = np.random.choice(data['unique_genes'], size = n_gene_sample, replace = False)
        used_hps = np.where(np.sum(data['gene_matrix'][:, np.in1d(data['unique_genes'], used_genes)], axis = 1) >  0)
        used_hp_seqs = data['unique_hp_seqs'][used_hps[0]]
    else:
        used_hp_seqs = data['unique_hp_seqs']
    if not np.isinf(args.max_CLs):
        n_CL_sample = np.min([args.max_CLs, len(data['unique_CLs'])])
        used_CLs = np.random.choice(data['unique_CLs'], size = n_CL_sample, replace = False)
    else:
        used_CLs = data['unique_CLs']
    new_dsets = [LFC_mat.ix[np.in1d(LFC_mat.index.values, used_hp_seqs),np.in1d(LFC_mat.columns.values, used_CLs)] for LFC_mat in data['LFC_mats'] if np.any(np.in1d(LFC_mat.columns.values, used_CLs))]
    data = preprocess.make_demeter2_data(new_dsets, sh_targets)

data_names = {'genes': data['unique_genes'],
             'CLs': data['unique_CLs'],
             'hps': data['unique_hp_seqs'],
             'seeds': data['unique_seed_seqs']}


#get train and test indices (random subset of datapoints)
if args.xv_frac == 0:
    test_inds = None
else:
    test_inds = []
    for LFC_mat in data['LFC_mats']:
        usable = np.where(~np.isnan(LFC_mat))
        n_usable = len(usable[0])
        test_set = np.random.choice(np.arange(n_usable), np.round(n_usable * args.xv_frac), replace = False)
        test_inds.append((usable[0][test_set], usable[1][test_set]))

#regularization params
reg_params = {
    'CL_l2_lambda': args.CL_l2_lambda,
    'rel_gene_l2_lambda': args.rel_gene_l2_lambda,
    'gene_l2_lambda': args.gene_l2_lambda,
    'seed_l2_lambda': args.seed_l2_lambda,
    'rel_seed_l2_lambda': args.rel_seed_l2_lambda,
    'hairpin_l2_lambda': args.hairpin_l2_lambda,
    'hp_unpred_l2_lambda': args.hp_unpred_l2_lambda
}

optim_params = {
    'precision': args.float_precision,
    'maxiter': args.maxiter,
    'print_freq': args.print_freq
}

gene_controls = {'pos': train_pos_con_genes, 'neg': train_neg_con_genes,
                'pos_test': test_pos_con_genes, 'neg_test': test_neg_con_genes}

print('Initializing model')
#init model
mod = demeter2.demeter(data['LFC_mats'], data['gene_matrix'], data['seed_matrix'], 
                   gene_sets = gene_controls, data_names = data_names, reg_params = reg_params, 
                   optim_params = optim_params, test_inds = test_inds, log_file = log_file)
print('Model initialized')

    
if fit_CL_slopes:
    #initialize overall slope estimates
    mod.init_slopes(data['LFC_mats'])

def fit_iter(ignore_test = False):
    if args.fit_shRNA_effs:
        mod.fit(data['LFC_mats'], fit_params = 'guide_effs', ignore_test = ignore_test) #estimate hairpin efficacies
    if args.ML_gene_slope:
        mod.fit(data['LFC_mats'], fit_params = 'gene_slopes_ML', ignore_test = ignore_test) #estimate gene slopes
    else:
        mod.fit(data['LFC_mats'], fit_params = 'gene_slopes', ignore_test = ignore_test) #estimate gene slopes
    if fit_CL_slopes:
        mod.fit(data['LFC_mats'], fit_params = 'ov_slopes_ML', ignore_test = ignore_test) #estimate overall slopes
    if args.fit_noise:
        mod.fit(data['LFC_mats'], fit_params = 'noise_vars', ignore_test = ignore_test) #estimate noise variances
    mod.fit(data['LFC_mats'], fit_params = 'scores', ignore_test = ignore_test) #re-estimate gene and seed scores

#get initial estimates of scores (gene and seed) and intercepts given slopes and noise
mod.fit(data['LFC_mats'], fit_params = 'scores')

#run blockwise optimization until convergence
cur_loss = mod.loss_evals[-1]
for ii in range(args.max_alt_iters):
    fit_iter(ignore_test = False)
    new_loss = mod.loss_evals[-1]
    frac_change = (cur_loss - new_loss) / np.max(np.array(np.abs([cur_loss, new_loss, 1])))
    mod.write('\nIt {}, prev loss: {}, new loss: {}, change: {}\n'.format(ii, cur_loss, new_loss, frac_change))
    cur_loss = new_loss
    if frac_change <= alt_tol:
        mod.write('Alternating optimization converged')
        break

if args.xv_frac > 0:
    mod.write('Fitting final iteration to all data')
    fit_iter(ignore_test = True)

#save model
mod.compute_R2_stats(data['LFC_mats'])
mod.write('saving MAP model to: {}'.format(model_path))
mod.save(model_path, save_perf_only = args.save_perf_only, edward = False)

if args.fit_edward:
    
    mod.write('Making edward model')
    total_elapsed()
    mod.fit(data['LFC_mats'], fit_params = 'noise_vars') #estimate noise variances
    mod.make_edward_model(data['LFC_mats'], data['gene_matrix'], data['seed_matrix'], data_names)
    mod.write('Training edward model')
    mod.edward_loss_evals = mod.run_edward_inference(n_iter = args.edward_iters)

    mod.write('saving edward model')
    mod.save(model_path, save_perf_only = args.save_perf_only, edward = True)
    total_elapsed()


