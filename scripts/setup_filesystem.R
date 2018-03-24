#################################
# Set up Directory Structure
#################################

# This setup file creates the directory structure necessary to run the pipeline and stores variable names globally
# so they can be accessed from any script that calls source('~/Documents/demeter2/scripts/pull_data.R')

# This location must be changed
base <- '/Users/zho/Desktop/' # Base for storing RNAi pipeline

# Downloaded raw DRIVE file from https://data.mendeley.com/datasets/y3ds55n88r/4
path_to_DRIVE_counts <- '/Users/zho/Desktop/rnai_pipeline/DRIVE/DriveCountData.RDS'

# Batch names
Achilles_batch_names <- c('98k', '55k_batch1', '55k_batch2')
DRIVE_batch_names <- c('BGPD', 'poolA', 'poolB')

# Directory structure for the pipeline
root <- paste0(base, 'other_models/')

# Make directories--everything is based on the root
datatypes <- c('Achilles', 'DRIVE', 'Marcotte', 'Combined')
paths <- paste0(root, datatypes, '/') %>% as.list() %>% setNames(paste0(datatypes, '_directory'))

# Set up GA directories
for (datatype in datatypes) {
  paths <- c(paths, paste0(root, datatype, '/GA/') %>% as.list() %>% setNames(paste0(datatype, '_GA_directory')))
  paths <- c(paths, paste0(root, datatype, '/GA/', datatype, '_GA_scores.csv') %>% as.list() %>% setNames(paste0(datatype, '_GA_fname')))
}

for (datatype in c('Achilles', 'DRIVE')) {
  # Set up basic RSA directories
  paths <- c(paths, paste0(root, datatype, '/RSA/') %>% as.list() %>% setNames(paste0(datatype, '_RSA_directory')))
  paths <- c(paths, paste0(root, datatype, '/RSA/inputs/') %>% as.list() %>% setNames(paste0(datatype, '_RSA_inputs_directory')))
  paths <- c(paths, paste0(root, datatype, '/RSA/outputs/') %>% as.list() %>% setNames(paste0(datatype, '_RSA_outputs_directory')))
  paths <- c(paths, paste0(root, datatype, '/RSA/', datatype, '_RSA_input.csv') %>% as.list() %>% setNames(paste0(datatype, '_RSA_input_fname')))
  paths <- c(paths, paste0(root, datatype, '/RSA/', datatype, '_RSA_scores.csv') %>% as.list() %>% setNames(paste0(datatype, '_RSA_scores_fname')))
  # Set up basic ATARiS directories
  paths <- c(paths, paste0(root, datatype, '/ATARiS/') %>% as.list() %>% setNames(paste0(datatype, '_ATARiS_directory')))
  paths <- c(paths, paste0(root, datatype, '/ATARiS/', datatype, '_ATARiS_scores.csv') %>% as.list() %>% setNames(paste0(datatype, '_ATARiS_scores_fname')))
  paths <- c(paths, paste0(root, datatype, '/ATARiS/output/') %>% as.list() %>% setNames(paste0(datatype, '_ATARiS_output_directory')))
  # Set up basic DEMETER1 directories
  paths <- c(paths, paste0(root, datatype, '/DEMETER1/') %>% as.list() %>% setNames(paste0(datatype, '_DEMETER1_directory')))
  paths <- c(paths, paste0(root, datatype, '/DEMETER1/', datatype, '_DEMETER1_scores.csv') %>% as.list() %>% setNames(paste0(datatype, '_DEMETER1_scores_fname')))
  paths <- c(paths, paste0(root, datatype, '/DEMETER1/', datatype, '_DEMETER1_batches.csv') %>% as.list() %>% setNames(paste0(datatype, '_DEMETER1_batches_fname')))
  paths <- c(paths, paste0(root, datatype, '/DEMETER1/output/') %>% as.list() %>% setNames(paste0(datatype, '_DEMETER1_output_directory')))
}

# Miscellaneous

# LFC directory for DRIVE data
paths$DRIVE_LFC_directory <- paste0(paths$DRIVE_directory, 'LFC/') # DRIVE LFC directory

# Raw DRIVE file, downloaded from https://data.mendeley.com/datasets/y3ds55n88r/4
paths$DRIVE_counts_fname <- path_to_DRIVE_counts

# Intermediate replicate median-collapsed DRIVE counts
paths$rep_medians_fname <- paste0(paths$DRIVE_directory, 'replicate_counts.Rds')

# Combined GA file for running on the cloud
paths$Combined_RDS_fname <- paste0(paths$Combined_GA_directory, 'Combined_GA_mat.RDS')

# ATARiS files, split by 55k and 98k batches for Achilles
paths$Achilles_55k_preATARiS_fname <- paste0(paths$Achilles_ATARiS_directory, 'Achilles_55k_preATARiS.gct')
paths$Achilles_98k_preATARiS_fname <- paste0(paths$Achilles_ATARiS_directory, 'Achilles_98k_preATARiS.gct')
paths$DRIVE_preATARiS_fname <- paste0(paths$DRIVE_ATARiS_directory, 'DRIVE_preATARiS.gct')

# DEMETER1 files, split by batch for Achilles
paths$Achilles_DEMETER1_mat_fnames <- paste0(paths$Achilles_DEMETER1_directory, 'Achilles_', Achilles_batch_names, '.gct')
paths$DRIVE_DEMETER1_mat_fname <- paste0(paths$DRIVE_DEMETER1_directory, 'DRIVE_mat.gct')

# Create the following directories if they do not already exist
dirnames <- c(root, paths[grep('directory', names(paths), value = T)]) %>% 
  l_ply(function(dirname) {
    ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
  })

# Create global variable

paths <<- paths
