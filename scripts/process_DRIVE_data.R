library(tidyverse)
library(plyr)
library(magrittr) 
library(reshape2)
library(tibble)
library(dplyr)
library(data.table)
library(taigr)
library(ggplot2)
library(ggrepel)
library(readr)
library(stringr)
library(cowplot)
library(psych)
library(knitr)
library(limma)
library(edgeR)

source('~/Documents/demeter2/scripts/pull_data.R')

# Process DRIVE Data for the pipeline. From raw counts to LFC estimate per hairpin/CL/batch.

# Collapse data by taking the median of replicate POOL/CLEANNAME/SHRNA_ID combinations

DRIVE_df <- read_rds(DRIVE_counts_fname)

rep_medians <- DRIVE_df %>%
  group_by(POOL, CLEANNAME, SEQ) %>%
  dplyr::summarise(PLASMID_COUNT = median(PLASMID_COUNT, na.rm=T),
                   SAMPLE_COUNT = median(SAMPLE_COUNT, na.rm=T)) %>%
  ungroup()

# Store intermediate file for backup
write_rds(rep_medians, rep_medians_fname)

# Read from intermediate file if necessary
# rep_medians <- read_rds(rep_medians_fname)

# Quality Control: remove all CLEANNAME f36p entries because it has an abnormally large number
# of sample counts that equal zero.

rep_medians %<>% filter(CLEANNAME != 'f36p')

# Get hairpin map
get_hairpin_map()

# Split rep_medians by pool and get normalized LFCS per hairpin/cl

prior_count <- 10
plasmid_threshold <- 20000000

rep_medians %>% d_ply(.(POOL), function(df) {
  
  PLASMID_counts <- acast(df, SEQ ~ CLEANNAME, value.var = 'PLASMID_COUNT')
  SAMPLE_counts <- acast(df, SEQ ~ CLEANNAME, value.var = 'SAMPLE_COUNT')
  
  stopifnot(all(rownames(PLASMID_counts) == rownames(SAMPLE_counts)))
  stopifnot(all(colnames(PLASMID_counts) == colnames(SAMPLE_counts)))
  
  good_samples <- which(colSums(PLASMID_counts, na.rm = T) >= plasmid_threshold) # samples with >= 20 million total plasmid counts
  virtualLibrary <- edgeR::equalizeLibSizes(PLASMID_counts[,good_samples])
  virtualLibrary <- rowMeans(virtualLibrary$pseudo.counts, na.rm = T) # average pseudocounts across samples
  virtualLibrary <- round(virtualLibrary) # round to nearest integer
  bad_samples <- setdiff(seq(ncol(PLASMID_counts)), good_samples)
  PLASMID_counts[,bad_samples] <- matrix(rep(virtualLibrary, length(bad_samples)), nrow = length(bad_samples), ncol = length(virtualLibrary), byrow = T)
  
  PLASMID_lcpm <- cpm(PLASMID_counts, log = T, prior.count = prior_count)
  SAMPLE_lcpm <- cpm(SAMPLE_counts, log = T, prior.count = prior_count)
  LFC <- SAMPLE_lcpm - PLASMID_lcpm

  LFC[PLASMID_lcpm < 1] <- NA # remove shRNAs with < 1 log plasmid counts per million
  
  # Convert cleannames to CCLE IDs
  colnames(LFC) <- ifelse(is.na(CleanCellLineName(colnames(LFC))), colnames(LFC), CleanCellLineName(colnames(LFC)))
  write.csv(LFC, paste0(DRIVE_LFC_directory, unique(df$POOL), '_LFC_mat.csv'))
  
})
