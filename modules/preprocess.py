import numpy as np
from scipy.sparse import csr_matrix, csr
import pandas as pd
import re
import taigapy


def get_hp_names(LFC_list):
    return(np.concatenate([df.index.values for df in LFC_list]))

def get_CL_names(LFC_list):
    return(np.concatenate([df.columns.values for df in LFC_list]))

def get_hp_batches(LFC_list):
    return(np.concatenate([ii * np.ones(df.shape[0]) for ii, df in enumerate(LFC_list)]))

def get_CL_batches(LFC_list):
    return(np.concatenate([ii * np.ones(df.shape[1]) for ii, df in enumerate(LFC_list)]))
 
def remove_promiscuous_hairpins(LFC_mats, sh_map, max_genes_targeted):
    '''Remove hairpins that target more than predefined number of genes from data matrices'''
    genes_targeted = sh_map.groupby('SEQ')['Gene_ID'].count()
    promisc_hps = genes_targeted[genes_targeted > max_genes_targeted].index.values
    print('Eliminating {} promiscuous hairpins'.format(len(promisc_hps)))
    LFC_data = [dset.ix[np.setdiff1d(dset.index.values, promisc_hps),:] for dset in LFC_mats]
    return(LFC_data)


def make_gene_matrix(sh_map, full_hp_set):
    '''Make sparse sh-to-gene mapping matrix'''
    #get unique IDs for each hp and gene
    un_hp_seqs, hp_inverse = np.unique(sh_map['SEQ'].values, return_inverse=True)
    unique_genes, gene_inverse = np.unique(sh_map['Gene_ID'].values.astype(str), return_inverse = True)

    #add in hairpins that dont appear in gene matrix
    add_hps = np.setdiff1d(full_hp_set, un_hp_seqs)
    un_hp_seqs = np.concatenate([un_hp_seqs, add_hps])

    gene_matrix = csr_matrix((np.ones_like(hp_inverse), (hp_inverse, gene_inverse)), shape = (len(un_hp_seqs), len(unique_genes)))
    return(gene_matrix, un_hp_seqs, unique_genes)


def make_seed_matrix(hp_seqs, seed_regions):
    '''Make sparse mapping from hairpins to unique seed sequences, using specified seed regions'''
    seed_seqs = np.empty([len(hp_seqs) * len(seed_regions)],dtype = 'object')
    seed_hp = np.empty([len(hp_seqs) * len(seed_regions)], dtype = 'int')
    cnt = 0
    for (start, stop) in seed_regions:
        for hps_i, hps in enumerate(hp_seqs):
            seed_seqs[cnt] = hps[start:stop]
            seed_hp[cnt] = hps_i
            cnt += 1
    unique_seeds, seed_inds = np.unique(seed_seqs, return_inverse = True)
    seed_matrix = csr_matrix((np.ones(len(seed_seqs)), (seed_hp, seed_inds)), 
                    shape=(len(hp_seqs), len(unique_seeds)))
    return(seed_matrix, unique_seeds)


def make_demeter2_data(LFC_list, target_df, max_genes_targeted = 10, seed_regions=((10,17),(11,18))):
    '''
    Make data inputs for DEMETER2.
    INPUTS:
        LFC_list: list of data matrices. Each with hairpin sequences as rownames, and cell lines as column names, entries are LFCs
        target_df: dataframe specifying hairpin-to-gene mapping as pairs of hairpin-gene relationships
        max_genes_targeted: Maximum number of genes a hairpin can map to before we remove it
        seed_regions: sequence regions defining the seed sequences used
    OUTPUTS:
        data dict with:
            LFC_mats: Same as inputs matrices but with bad hps removed
            gene_matrix: sparse matrix mapping hairpins to genes
            seed_matrix: sparse matrix mapping hairpins to seed sequences
            unique_hp_seqs: hp sequence of unique hps (rows of gene and seed matrices)
            unique_genes: genes making columns of gene matrix
            unique_seed_seqs: unique seed sequences making columns of seed_matrix
    '''
    LFC_mats = remove_promiscuous_hairpins(LFC_list, target_df, max_genes_targeted = max_genes_targeted)

    all_hps = get_hp_names(LFC_mats)

    target_df = target_df.ix[np.in1d(target_df['SEQ'].values, all_hps), ['SEQ', 'Gene_ID']]

    init = target_df.shape[0]
    target_df = target_df[~target_df['Gene_ID'].str.contains('NO_CURRENT')]
    print('Eliminated {}/{} non-targeting hairpins from map'.format(init-target_df.shape[0], init))

    target_df = collapse_identical_genes(target_df)

    gene_matrix, unique_hp_seqs, unique_genes = make_gene_matrix(target_df, all_hps)

    seed_matrix, unique_seed_seqs = make_seed_matrix(unique_hp_seqs, seed_regions)

    unique_CLs = np.unique(get_CL_names(LFC_mats))
    print('Creating dataset with:\n{} hairpins\n{} CLs\n{} genes\n{} seeds'.format(
        len(unique_hp_seqs),
        len(unique_CLs),
        len(unique_genes),
        len(unique_seed_seqs)))

    return({'LFC_mats': LFC_mats, 
            'gene_matrix': gene_matrix, 
            'seed_matrix': seed_matrix, 
            'unique_hp_seqs': unique_hp_seqs, 
            'unique_seed_seqs': unique_seed_seqs, 
            'unique_genes': unique_genes,
            'unique_CLs': unique_CLs})


def read_gct(filename):
    # given a GCT, return HairpinsAndGenes
    with open(filename, "rt") as fd:
        df = pd.read_csv(fd, delimiter="\t", skiprows=2, index_col=0)
        return df.drop("Description", axis=1)


def collapse_identical_genes(sh_map):
    '''Find genes targeted by identical groups of hairpins and collapse them to single 'gene-family' representation'''
    
    new_sh_map = sh_map.copy()

    #get set of hps targeting each gene as a concatenated, ordered string
    a = sh_map.sort_values(['SEQ'])
    a = sh_map.groupby('Gene_ID')['SEQ'].unique() 
    a = a.apply(lambda x: '_'.join(x.astype(str)))

    #find groups of genes corresponding to recurrent hp-sets
    avc = a.value_counts()
    gf_hp_groups = avc[avc > 1].index.values
    a_gf = a[np.in1d(a.values, gf_hp_groups)].to_frame().reset_index()
    gene_families = a_gf.groupby('SEQ')['Gene_ID'].unique().values

    print('Identified {} gene families'.format(len(gene_families)))

    #make map from genes to gene-families
    gf_ID_map = {}
    for gf in gene_families:
        for gene in gf:
            gf_ID_map[gene] = '&'.join(gf.astype(str))

    #add in singlets to map
    other_genes = [x for x in sh_map['Gene_ID'].values if x not in gf_ID_map]
    gf_ID_map.update({gene: gene for gene in other_genes})

    #apply map
    new_sh_map['Gene_ID'] = new_sh_map['Gene_ID'].map(gf_ID_map)
    new_sh_map.drop_duplicates(inplace = True)
    return(new_sh_map)
