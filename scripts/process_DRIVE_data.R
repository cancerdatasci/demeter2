library(tidyverse)
library(plyr)
library(magrittr) 
library(reshape2)
library(tibble)
library(dplyr)
library(data.table)
library(ggplot2)
library(readr)
library(stringr)
library(edgeR)

source('~/Documents/github/demeter2/scripts/setup_filesystem.R')

# Process DRIVE Data for the pipeline. From raw counts to LFC estimate per hairpin/CL/batch.
DRIVE_df <- read_rds(paths$DRIVE_counts_fname)

# Quality Control: remove all CLEANNAME f36p entries because it has an abnormally large number
# of sample counts that equal zero.
DRIVE_df %<>% filter(CLEANNAME != 'f36p')

exp_info <- DRIVE_df %>% 
  dplyr::select(EXPERIMENT_ID, POOL, CLEANNAME) %>% 
  unique()

exp_info$CCLE_ID <- ifelse(is.na(CleanCellLineName(exp_info$CLEANNAME)), 
                           exp_info$CLEANNAME, 
                           CleanCellLineName(exp_info$CLEANNAME))

hp_info <- DRIVE_df %>% 
  dplyr::select(SHRNA_ID, SEQ) %>% 
  unique()

# Another update for cell line names in DRIVE
new_name_table <- read.csv('~/Documents/demeter2/name_change_map.csv', check.names = FALSE, stringsAsFactors = FALSE)
new_name_map <- new_name_table$new_name %>% set_names(new_name_table$old_name)

prior_count <- 10
plasmid_threshold <- 20000000

LFCs <- DRIVE_df %>% dlply(.(POOL), function(df) {
  
  PLASMID_counts <- acast(df, SHRNA_ID ~ EXPERIMENT_ID, value.var = 'PLASMID_COUNT')
  SAMPLE_counts <- acast(df, SHRNA_ID ~ EXPERIMENT_ID, value.var = 'SAMPLE_COUNT')
  
  stopifnot(all(rownames(PLASMID_counts) == rownames(SAMPLE_counts)))
  stopifnot(all(colnames(PLASMID_counts) == colnames(SAMPLE_counts)))
  
  good_samples <- which(colSums(PLASMID_counts, na.rm = T) >= plasmid_threshold) # samples with >= 20 million total plasmid counts
  virtualLibrary <- edgeR::equalizeLibSizes(PLASMID_counts[,good_samples])
  virtualLibrary <- rowMeans(virtualLibrary$pseudo.counts, na.rm = T) # average pseudocounts across samples
  virtualLibrary <- round(virtualLibrary) # round to nearest integer
  bad_samples <- setdiff(seq(ncol(PLASMID_counts)), good_samples)
  PLASMID_counts[,bad_samples] <- t(matrix(rep(virtualLibrary, length(bad_samples)), nrow = length(bad_samples), ncol = length(virtualLibrary), byrow = T))
  
  PLASMID_lcpm <- cpm(PLASMID_counts, log = T, prior.count = prior_count)
  SAMPLE_lcpm <- cpm(SAMPLE_counts, log = T, prior.count = prior_count)
  LFC <- SAMPLE_lcpm - PLASMID_lcpm

  LFC[PLASMID_lcpm < 1] <- NA # remove shRNAs with < 1 log plasmid counts per million
  
  long <- melt(LFC) %>% 
    set_colnames(c('SHRNA_ID', 'EXPERIMENT_ID', 'LFC')) %>% 
    left_join(exp_info, by = 'EXPERIMENT_ID') %>% 
    left_join(hp_info, by = 'SHRNA_ID')
  
  # Median collapse replicates in two dimensions
  # 1. EXPERIMENT_IDs with the same POOL/CL combination
  # 2. SHRNA_IDs with the same SEQ
  rep_medians <- long %>% 
    group_by(POOL, CCLE_ID, SEQ) %>% 
    dplyr::summarise(median_LFC = median(LFC, na.rm = T))
  median_LFC_mat <- rep_medians %>% acast(SEQ ~ CCLE_ID)
  
  # update colnames in drive
  colnames(median_LFC_mat) %<>% plyr::revalue(new_name_map)
  
  write.csv(as.data.frame(median_LFC_mat), paste0(paths$DRIVE_LFC_directory, unique(long$POOL), '_LFC_mat.csv'))
  return(median_LFC_mat)
})
