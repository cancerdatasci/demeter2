# This file's purpose is two-fold:
# 1. It creates the directory structure necessary to run the pipeline and stores variable names globally so they can be accessed from any script
# that calls source('~/Documents/demeter2/scripts/pull_data.R')
# 2. It provides functions to pulls various files from Taiga to run the pipeline
library(taigr)

# For CCLE annotations
load('~/Documents/DRIVE ATARiS Validation/Annotations.RData')

# These locations must be changed
base <- '/Users/zho/Desktop/' # Base for storing RNAi pipeline--this should be an existing location

# Location of demeter files for kubeque
kubeque_d1_fname <- paste0(base, 'demeter/kubeque_d1_files')

# find sequences that target too many genes (>10) and filter them out
remove_promiscuous_seqs <- function(map, datatype) {
  if (datatype == 'Achilles' || datatype == 'DRIVE') {
    promiscuous_seqs <- map %>% 
      dplyr::select(seq, Gene = Entrez_ID)
  } else if (datatype == 'Marcotte') {
    promiscuous_seqs <- map %>% 
      dplyr::select(seq, Gene = Gene_Entrez_ID)
  }
  promiscuous_seqs %<>% 
    unique() %>% 
    group_by(seq) %>% 
    dplyr::summarise(num_genes_targeted = n()) %>% 
    filter(num_genes_targeted > 10) %>% 
    .[['seq']]
  
  clean_map <- map[!(map$seq %in% promiscuous_seqs),]
  print(paste0('Removed ', nrow(map) - nrow(clean_map), ' promiscuous sequences from ', datatype, ' hairpin map'))
  return(clean_map)
}

# Temporary function to deal with mapping discrepancies -- returns first Gene/Entrez ID pair found
get_first_match <- function(map) {
  pairs <- map %>% select(Gene, Entrez_ID) %>% unique()
  good_pairs <- pairs[match(unique(pairs$Entrez_ID), pairs$Entrez_ID),]
  print(paste0('get_first_match removed ', nrow(pairs) - nrow(good_pairs), ' entries'))
  return(map %>% filter(Entrez_ID %in% good_pairs$Entrez_ID))
}

# Deal with hairpin map. To be modified
get_hairpin_map <- function() {
  map <- load.from.taiga(data.name='gpp-shrna-mapping-8759', data.version=1, data.file='CP1175_20171102_19mer')
  map %<>% select(Gene = `Gene Symbol`, Entrez_ID = `Gene ID`, seq = `Barcode Sequence`)
  map <- map[-(grep('NO_CUR', map$Gene)),] # filter out any genes containing "NO_CURRENT..."
  hairpin_map <<- remove_promiscuous_seqs(map, 'Achilles')
}

# Directory that contains everything needed to run models
models_directory <- paste0(base, 'rnai_pipeline_models/')

# Directory structure for the pipeline
root <- paste0(base, 'rnai_pipeline/') # Root to the entire pipeline

# Directories--everything is based on the root
DRIVE_directory <- paste0(root, 'DRIVE/') # DRIVE directory
Achilles_directory <- paste0(root, 'Achilles/') # Achilles directory
Marcotte_directory <- paste0(root, 'Marcotte/') # Marcotte directory
combined_directory <- paste0(root, 'Combined/') # Combined directory

DRIVE_LFC_directory <- paste0(DRIVE_directory, 'LFC/') # DRIVE LFC directory
# Achilles_LFC_directory <- paste0(Achilles_directory, 'LFC/') # Achilles LFC directory
# Marcotte_LFC_directory <- paste0(Marcotte_directory, 'LFC/') # Marcotte LFC directory

DRIVE_GA_directory <- paste0(DRIVE_directory, 'GA/') # DRIVE Gene Averaging directory
Achilles_GA_directory <- paste0(Achilles_directory, 'GA/') # Achilles Gene Averaging directory
Marcotte_GA_directory <- paste0(Marcotte_directory, 'GA/') # Marcotte Gene Averaging directory
combined_GA_directory <- paste0(combined_directory, 'GA/') # Combined Gene Averaging directory

DRIVE_RSA_directory <- paste0(DRIVE_directory, 'RSA/') # DRIVE RSA Results directory
Achilles_RSA_directory <- paste0(Achilles_directory, 'RSA/') # Achilles RSA Results directory
Marcotte_RSA_directory <- paste0(Marcotte_directory, 'RSA/') # Marcotte RSA Results directory

DRIVE_RSA_inputs_directory <- paste0(DRIVE_RSA_directory, 'inputs/') # DRIVE RSA Inputs directory
Achilles_RSA_inputs_directory <- paste0(Achilles_RSA_directory, 'inputs/') # Achilles RSA Inputs directory
Marcotte_RSA_inputs_directory <- paste0(Marcotte_RSA_directory, 'inputs/') # Marcotte RSA Inputs directory

DRIVE_RSA_outputs_directory <- paste0(DRIVE_RSA_directory, 'outputs/') # DRIVE RSA outputs directory
Achilles_RSA_outputs_directory <- paste0(Achilles_RSA_directory, 'outputs/') # Achilles RSA outputs directory
Marcotte_RSA_outputs_directory <- paste0(Marcotte_RSA_directory, 'outputs/') # Marcotte RSA outputs directory

DRIVE_ATARiS_directory <- paste0(DRIVE_directory, 'ATARiS/') # DRIVE ATARiS Results directory
Achilles_ATARiS_directory <- paste0(Achilles_directory, 'ATARiS/') # Achilles ATARiS Results directory
# Marcotte_ATARiS_directory <- paste0(Marcotte_directory, 'ATARiS/') # Marcotte ATARiS Results directory

DRIVE_DEMETER1_directory <- paste0(DRIVE_directory, 'DEMETER1/') # DRIVE DEMETER1 Results directory
Achilles_DEMETER1_directory <- paste0(Achilles_directory, 'DEMETER1/') # Achilles DEMETER1 Results directory
# Marcotte_DEMETER1_directory <- paste0(Marcotte_directory, 'DEMETER1/') # Marcotte DEMETER1 Results directory

# DRIVE_DEMETER2_directory <- paste0(DRIVE_directory, 'DEMETER2/') # DRIVE DEMETER2 Results directory
# Achilles_DEMETER2_directory <- paste0(Achilles_directory, 'DEMETER2/') # Achilles DEMETER2 Results directory
# Marcotte_DEMETER2_directory <- paste0(Marcotte_directory, 'DEMETER2/') # Marcotte DEMETER2 Results directory

# Batch names
Achilles_batch_names <- c('98k', '55k_batch1', '55k_batch2')
DRIVE_batch_names <- c('BGPD', 'poolA', 'poolB')

# Raw DRIVE file
DRIVE_counts_fname <- paste0(DRIVE_directory, 'DriveCountData.RDS')

# Intermediate replicate median-collapsed DRIVE counts
rep_medians_fname <- paste0(DRIVE_directory, 'replicate_counts.Rds')

# Gene Averaging files

Achilles_GA_fname <- paste0(Achilles_GA_directory, 'Achilles_GA_scores.csv')
DRIVE_GA_fname <- paste0(DRIVE_GA_directory, 'DRIVE_GA_scores.csv')
Marcotte_GA_fname <- paste0(Marcotte_GA_directory, 'Marcotte_GA_scores.csv')
Achilles_DRIVE_Marcotte_RDS_fname <- paste0(combined_GA_directory, 'Achilles_DRIVE_Marcotte_GA_mat.RDS')

# RSA files
Achilles_RSA_input_fname <- paste0(Achilles_RSA_directory, 'Achilles_RSA_input.csv')
DRIVE_RSA_input_fname <- paste0(DRIVE_RSA_directory, 'DRIVE_RSA_input.csv')
Achilles_RSA_scores_fname <- paste0(Achilles_RSA_directory, 'Achilles_RSA_scores.csv')

DRIVE_RSA_scores_fname <- paste0(DRIVE_RSA_directory, 'DRIVE_RSA_scores.csv')

# ATARiS files

Achilles_55k_preATARiS_fname <- paste0(Achilles_ATARiS_directory, 'Achilles_55k_preATARiS.gct')
Achilles_98k_preATARiS_fname <- paste0(Achilles_ATARiS_directory, 'Achilles_98k_preATARiS.gct')

DRIVE_preATARiS_fname <- paste0(DRIVE_ATARiS_directory, 'DRIVE_preATARiS.gct')

Achilles_55k_intermediate_fname <- paste0(Achilles_directory, 'Achilles_55k_intermediate.RDS')
Achilles_98k_intermediate_fname <- paste0(Achilles_directory, 'Achilles_98k_intermediate.RDS')
DRIVE_intermediate_fname <- paste0(DRIVE_directory, 'DRIVE_intermediate.RDS')

Achilles_ATARiS_scores_fname <- paste0(Achilles_ATARiS_directory, 'Achilles_ATARiS_scores.csv')

DRIVE_ATARiS_scores_fname <- paste0(DRIVE_ATARiS_directory, 'DRIVE_ATARiS_scores.csv')

# DEMETER1 files

Achilles_DEMETER1_batches_fname <- paste0(Achilles_DEMETER1_directory, 'Achilles_batches.csv')
DRIVE_DEMETER1_batches_fname <- paste0(DRIVE_DEMETER1_directory, 'DRIVE_batches.csv')

Achilles_DEMETER1_mat_fnames <- paste0(Achilles_DEMETER1_directory, 'Achilles_', Achilles_batch_names, '.gct')
Achilles_DEMETER1_scores_fname <- paste0(Achilles_DEMETER1_directory, 'Achilles_DEMETER1_scores.csv')
DRIVE_DEMETER1_mat_fname <- paste0(DRIVE_DEMETER1_directory, 'DRIVE_mat.gct')
DRIVE_DEMETER1_scores_fname <- paste0(DRIVE_DEMETER1_directory, 'DRIVE_DEMETER1_scores.csv')

# Create the following directories if they do not already exist

dirnames <- c(root, DRIVE_directory, Achilles_directory, Marcotte_directory, combined_directory,
              DRIVE_LFC_directory,
              DRIVE_GA_directory, Achilles_GA_directory, Marcotte_GA_directory, combined_GA_directory,
              DRIVE_RSA_directory, Achilles_RSA_directory, Marcotte_RSA_directory,
              DRIVE_RSA_inputs_directory, Achilles_RSA_inputs_directory, Marcotte_RSA_inputs_directory,
              DRIVE_RSA_outputs_directory, Achilles_RSA_outputs_directory, Marcotte_RSA_outputs_directory,
              DRIVE_ATARiS_directory, Achilles_ATARiS_directory,
              DRIVE_DEMETER1_directory, Achilles_DEMETER1_directory)

dirnames %>% l_ply(function(dirname) {
  ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
})

# Load LFC Matrices, stored on Taiga

get_LFC_matrices <- function() {
  # Achilles
  a98k_LFC_mat <<- load.from.taiga(data.name='achilles-98k-repcollapsed-lfc-19ce', data.version=1)
  a55k1_LFC_mat <<- load.from.taiga(data.name='achilles-55k-batch1-repcollapsed-lfc-d708', data.version=1)
  a55k2_LFC_mat <<- load.from.taiga(data.name='achilles-55k-batch2-repcollapsed-lfc-bd7f', data.version=1)
  # DRIVE
  BGPD_LFC_mat <<- load.from.taiga(data.name='drive-lfc-matrices-3867', data.version=4, data.file='BGPD_LFC_mat')
  poolA_LFC_mat <<- load.from.taiga(data.name='drive-lfc-matrices-3867', data.version=4, data.file='poolA_LFC_mat')
  poolB_LFC_mat <<- load.from.taiga(data.name='drive-lfc-matrices-3867', data.version=4, data.file='poolB_LFC_mat')
  # Marcotte
  # Preprocess Marcotte LFC matrix
  Marcotte_LFC_mat <- load.from.taiga(data.name='marcotte-demeter-z-scored-expanded-gene-sols-e905', data.version=7, data.file='log.effect.removed.bad')
  # Convert cleannames to CCLE IDs
  colnames(Marcotte_LFC_mat) <- ifelse(is.na(CleanCellLineName(colnames(Marcotte_LFC_mat))), colnames(Marcotte_LFC_mat), CleanCellLineName(colnames(Marcotte_LFC_mat)))
  Marcotte_LFC_mat <<- Marcotte_LFC_mat
}