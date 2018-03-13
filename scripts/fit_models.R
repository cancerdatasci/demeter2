library(plyr)
library(magrittr) 
library(reshape2)
library(tibble)
library(dplyr)
library(data.table)
library(taigr)
library(readr)
library(stringr)
library(dependr)

source('~/Documents/demeter2/scripts/pull_data.R')

# Retrieve data for all models
get_hairpin_map()
get_LFC_matrices()
Achilles_LFC_mats <- list(a98k_LFC_mat, a55k1_LFC_mat, a55k2_LFC_mat)
DRIVE_LFC_mats <- list(BGPD_LFC_mat, poolA_LFC_mat, poolB_LFC_mat)

# Fit each model on each datatype

############################################## Gene Averaging ##############################################

preprocess_for_GA <- function(LFC) {
  LFC %<>% scale(center = TRUE, scale = TRUE) # zscore normalize by CL
  # Convert to longform
  LFC %<>% as.data.frame() %>%
    rownames_to_column(var = 'seq') %>%
    melt(id.vars = 'seq', variable.name = 'CL', value.name = 'LFC')
  LFC %<>% left_join(hairpin_map, by = 'seq')
  LFC$Gene_Entrez_ID <- paste0(LFC$Gene, ' (', LFC$Entrez_ID, ')')
  return(LFC %>% select(CL, LFC, Gene_Entrez_ID))
}

fit_GA <- function(LFC, fname) {
  GA <- LFC %>%
    group_by(CL, Gene_Entrez_ID) %>%
    dplyr::summarise(avg_score = mean(LFC, na.rm=T), num_shRNAs = n())
  # ensure >1 hairpin/gene
  GA <- GA[GA$num_shRNAs > 1,]
  GA_mat <- GA %>% acast(Gene_Entrez_ID ~ CL, value.var = 'avg_score')
  write.csv(GA_mat, fname)
}

GA_wrapper <- function(datatype, fname, LFC_mats) {
  LFC <- mapply(preprocess_for_GA, LFC = LFC_mats, SIMPLIFY = FALSE)
  LFC <- do.call(rbind, LFC)
  fit_GA(LFC, fname)
}

GA_wrapper('Achilles', Achilles_GA_fname, Achilles_LFC_mats)
GA_wrapper('DRIVE', DRIVE_GA_fname, DRIVE_LFC_mats)
GA_wrapper('Marcotte', Marcotte_GA_fname, list(Marcotte_LFC_mat))

#### Achilles + DRIVE + Marcotte GA takes lots of memory to run. We ran this remotely.
#### Our pipeline involves uploading an intermeditate RDS file to Google Cloud, as shown below.

#### Our approach to running the GA model for the combined data

# GA_combination <- function() {
#   LFC_mats <- list(a98k_LFC_mat, a55k1_LFC_mat, a55k2_LFC_mat, BGPD_LFC_mat, poolA_LFC_mat, poolB_LFC_mat, Marcotte_LFC_mat)
#   LFC <- mapply(preprocess_for_GA, LFC = LFC_mats, SIMPLIFY = FALSE)
#   LFC <- do.call(rbind, LFC)
#   write_rds(LFC, 'Achilles_DRIVE_Marcotte_GA_mat.RDS')
# }
# 
# fit_GA <- function(LFC, fname) {
#   LFC %<>% mutate(CL = as.character(CL))
#   # Convert to CCLE names
#   LFC$CL <- ifelse(is.na(CleanCellLineName(LFC$CL)), LFC$CL, CleanCellLineName(LFC$CL))
#   GA <- LFC %>%
#     group_by(CL, Gene_Entrez_ID) %>%
#     dplyr::summarise(avg_score = mean(LFC, na.rm=T), num_shRNAs = n())
#   # ensure >1 hairpin/gene
#   GA <- GA[GA$num_shRNAs > 1,]
#   GA_mat <- GA %>% acast(Gene_Entrez_ID ~ CL, value.var = 'avg_score')
#   
#   write.csv(GA_mat, fname)
# }
# 
# GA_combination()
# LFC <- read_rds('Achilles_DRIVE_Marcotte_GA_mat.RDS')
# fit_GA(LFC, 'Achilles_DRIVE_Marcotte_GA_scores.csv')

############################################## Shared Functions ##############################################

# Make sure X doesn't get appended to beginning of column names
check_CCLE_colnames <- function(mat) {
  stopifnot(length(grep('^X[0-9]', colnames(mat))) == 0)
}

# Normalizes LFC matrix by CL and changes rownames to reflect the sequence, entrez id, and batch name
process_LFC <- function(LFC, datatype, batch_name, distinct_hp) {
  LFC %<>% scale(center = TRUE, scale = TRUE) # zscore normalize by CL
  # hairpin_map %<>% mutate(Gene_Entrez_ID = paste0(Gene, '_', Entrez_ID))
  cur <- hairpin_map %>% filter(seq %in% rownames(LFC)) %>% select(seq, Gene, Entrez_ID) %>% unique()
  new_LFC <- LFC[cur$seq, ] # Split into separate rows for hairpins targeting each gene
  if (distinct_hp) {
    # hairpins from different batches should be treated as distinct
    rownames(new_LFC) <- paste0(cur$seq, '@', cur$Gene, '@', cur$Entrez_ID, '@', batch_name)
  } else {
    # hairpins from different batches should be treated as the same
    rownames(new_LFC) <- paste0(cur$seq, '@', cur$Gene, '@', cur$Entrez_ID)
  }
  return(new_LFC)
}

############################################## RSA ##############################################

# Create a file for each cell line, run RSA for each cell line, and stitch outputs together
preprocess_for_RSA <- function(LFC_mats, dir) {
  LFC <- LFC_mats %>% llply(function(mat) {
    mat %<>% scale(center = TRUE, scale = TRUE) # zscore normalize by CL
    long <- melt(mat) %>% set_colnames(c('seq', 'CL', 'Score'))
    return(long)
  })
  long <- do.call(rbind, LFC)
  # Use "@" as a delimiter
  long %<>% left_join(hairpin_map, by = 'seq') %>%
    mutate(Entrez_ID_CL = paste0(Entrez_ID, '@', CL)) %>% 
    select(CL, Gene_ID = Entrez_ID_CL, Well_ID = seq, Score)
  long %>% d_ply(.(CL), function(df) {
    print(df$CL %>% unique())
    # Filter out genes targeted by only one hairpin
    gene_count <- table(df$Gene_ID)
    good_genes <- gene_count[gene_count > 1] %>% names()
    temp <- df %>% filter(Gene_ID %in% good_genes) %>% select(Gene_ID, Well_ID, Score)
    write.csv(temp, paste0(dir, unique(df$CL), '.csv'), row.names = F)
  })
}

postprocess_RSA_output <- function(dir, fname) {
  files <- list.files(dir, full.names = T) %>% as.list()
  RSA <- files %>% llply(function(f) {
    df <- read.csv(f, check.names = F)
    new <- df %>% cbind(colsplit(as.character(df$Gene_ID),"@",c("Entrez_ID","CL"))) %>%
      dplyr::select(CL, Entrez_ID, Score = LogP) %>%
      left_join(hairpin_map, by = 'Entrez_ID') %>%
      mutate(Gene_Entrez_ID = paste0(Gene, ' (', Entrez_ID, ')')) %>%
      select(CL, Gene_Entrez_ID, Score) %>% unique()
    return(new)
  })
  RSA <- do.call(rbind, RSA)
  mat <- RSA %>% acast(Gene_Entrez_ID ~ CL)
  write.csv(mat, fname)
}

preprocess_for_RSA(DRIVE_LFC_mats, DRIVE_RSA_inputs_directory)
preprocess_for_RSA(Achilles_LFC_mats, Achilles_RSA_inputs_directory)

#### INTERMEDIATE STEP: Use scripts in run_RSA to generate RSA scores for each cell line ####

postprocess_RSA_output(DRIVE_RSA_outputs_directory, DRIVE_RSA_scores_fname)
postprocess_RSA_output(Achilles_RSA_outputs_directory, Achilles_RSA_scores_fname)


############################################## ATARiS ##############################################

write_ATARiS_gct <- function(df, gct_fname) {
  gene_entrez_ids <- rownames(df) %>% llply(function(str) {
    chunks <- strsplit(str, '@')[[1]]
    return(paste0(chunks[2], ' (', chunks[3], ')'))
  }) %>% unlist()
  rownames(df) <- gsub('@', '_', rownames(df)) # Convert back to "_" for ATARiS
  check_CCLE_colnames(df)
  dependr::write.gct(df, gct_fname, gene_entrez_ids)
}

preprocess_for_ATARiS <- function(LFC_mats, datatype, batch_names) {
  distinct_hp <- ifelse(datatype == 'Achilles', FALSE, TRUE)
  mats <- mapply(process_LFC, LFC = LFC_mats, datatype = datatype, batch_name = batch_names, distinct_hp)
  if (datatype == 'DRIVE') {
    LFC <- melt(mats[[1]], value.name = 'LFC') %>% rbind(melt(mats[[2]], value.name = 'LFC')) %>% rbind(melt(mats[[3]], value.name = 'LFC')) %>% mutate(Var1 = as.character(Var1))
    drive_df <- acast(LFC, Var1 ~ Var2, value.var = 'LFC')
    write_ATARiS_gct(drive_df, DRIVE_preATARiS_fname) 
  } else if (datatype == 'Achilles') {
    # Achilles 98k batch
    a98k_LFC <- melt(mats[[1]], value.name = 'LFC')
    a98k_df <- acast(a98k_LFC, Var1 ~ Var2, value.var = 'LFC')
    write_ATARiS_gct(a98k_df, Achilles_98k_preATARiS_fname)
    # Achilles 55k batches
    a55k_LFC <- melt(mats[[2]], value.name = 'LFC') %>% rbind(melt(mats[[3]], value.name = 'LFC'))
    a55k_df <- acast(a55k_LFC, Var1 ~ Var2, value.var = 'LFC')
    write_ATARiS_gct(a55k_df, Achilles_55k_preATARiS_fname)
  }
}

# Save output in corresponding ATARiS directory

helper_filter <- function(fname) {
  df <- read.gct(fname)
  genes <- str_match(rownames(df), '.+\\)') %>% as.vector()
  rownames(df) <- genes
  # Deal with unwanted "X" preceding cell line names
  colnames(df) <- gsub("^X", "", colnames(df))
  # Keep first ATARiS solution for each gene
  filtered_df <- df[match(unique(genes), genes),]
  return(filtered_df)
}

# Stitches output from 55k and 98k batches to form a single matrix
postprocess_ATARiS_output <- function(gct_files, ataris_score_fname, stitch) {
  if (stitch) {
    mats <- mapply(helper_filter, fname = gct_files)
    scores <- melt(mats[[1]], value.name = 'score') %>% rbind(melt(mats[[2]], value.name = 'score'))
    df <- acast(scores, Var1 ~ Var2, value.var = 'score')
  } else {
    df <- helper_filter(gct_files)
  }
  check_CCLE_colnames(df)
  write.csv(df, ataris_score_fname)
}

preprocess_for_ATARiS(Achilles_LFC_mats, 'Achilles', Achilles_batch_names)
preprocess_for_ATARiS(DRIVE_LFC_mats, 'DRIVE', DRIVE_batch_names)

#### INTERMEDIATE STEP: Run ATARiS on Gene Pattern ####

postprocess_ATARiS_output(c(paste0(Achilles_ATARiS_directory, 'Achilles_55k.Gs.gct'), paste0(Achilles_ATARiS_directory, 'Achilles_98k.Gs.gct')), Achilles_ATARiS_scores_fname, TRUE)
postprocess_ATARiS_output(paste0(DRIVE_ATARiS_directory, 'DRIVE.Gs.gct'), DRIVE_ATARiS_scores_fname, FALSE)

############################################## DEMETER 1 ##############################################

# Create batch files for input to DEMETER1
make_batches_file <- function(LFC_mats, datatype) {
  CL <- c(colnames(LFC_mats[[1]]), colnames(LFC_mats[[2]]), colnames(LFC_mats[[3]]))
  if (datatype == 'Achilles') {
    batches <- data.frame(Name = CL, `DEMETER batch` = c(rep(1, ncol(LFC_mats[[1]])), rep(2, ncol(LFC_mats[[2]])), rep(3, ncol(LFC_mats[[3]]))), check.names = FALSE)
    fname <- Achilles_DEMETER1_batches_fname
  } else if (datatype == 'DRIVE') {
    # All CLs in DRIVE data will be considered to be in the same batch
    batches <- data.frame(Name = CL, `DEMETER batch` = rep(1, length(CL)), check.names = FALSE)
    fname <- DRIVE_DEMETER1_batches_fname
  }
  write.csv(batches, fname, row.names = FALSE)
}

write_DEMETER1_gct <- function(LFC, fname, datatype) {
  if (datatype == 'Achilles') {
    # Format: SEQ_ENTREZID
    rownames(LFC) <- gsub('@(.*)@', "_", rownames(LFC)) %>% make.names(unique = TRUE)
  } else if (datatype == 'DRIVE') {
    names <- rownames(LFC) %>% llply(function(str) {
      chunks <- strsplit(str, '@')[[1]]
      # Format: SEQ.POOL_ENTREZID
      return(paste0(chunks[1], '.', chunks[4], '_', chunks[3]))
    }) %>% unlist()
    rownames(LFC) <- names %>% make.names(unique = T)
  }
  check_CCLE_colnames(LFC)
  dependr::write.gct(LFC, fname, rownames(LFC))
}

# Write Achilles and DRIVE gct files for preprocessing. Achilles will have 3 files, DRIVE will have 1
preprocess_for_DEMETER1 <- function(LFC_mats, datatype, batch_names, mat_fnames) {
  
  make_batches_file(LFC_mats, datatype)    
  distinct_hp <- ifelse(datatype == 'Achilles', FALSE, TRUE)
  
  if (datatype == 'Achilles') {
    mats <- mapply(process_LFC, LFC = LFC_mats, datatype = datatype, batch_name = batch_names, distinct_hp)
    mapply(write_DEMETER1_gct, LFC = mats, fname = mat_fnames, datatype = datatype)
  } else if (datatype == 'DRIVE') {
    mats <- mapply(process_LFC, LFC = LFC_mats, datatype = datatype, batch_name = batch_names, distinct_hp)
    LFC <- melt(mats[[1]], value.name = 'LFC') %>% rbind(melt(mats[[2]], value.name = 'LFC')) %>% rbind(melt(mats[[3]], value.name = 'LFC')) %>% mutate(Var1 = as.character(Var1))
    df <- acast(LFC, Var1 ~ Var2, value.var = 'LFC')
    write_DEMETER1_gct(df, mat_fnames, datatype = datatype)
  }
}

# Convert rownames to Gene Symbol (Entrez ID)
postprocess_DEMETER1_output <- function(input_fname, output_fname) {
  df <- read.csv(input_fname, check.names = F, row.names = 1)
  # Deal with DEMETER1 parsing issue. Entrez IDs are after the last "&" in some rows.
  entrez_ids <- ifelse(grepl('&', rownames(df)), gsub(')', '', gsub('.+&', '', rownames(df))), rownames(df))
  map <- hairpin_map %>% select(Gene, Entrez_ID) %>% unique()
  rownames(map) <- map$Entrez_ID
  final_genes <- map[entrez_ids,]
  rownames(df) <- paste0(final_genes$Gene, ' (', final_genes$Entrez_ID, ')')
  write.csv(df, output_fname)
}

preprocess_for_DEMETER1(Achilles_LFC_mats, 'Achilles', Achilles_batch_names, Achilles_DEMETER1_mat_fnames)
preprocess_for_DEMETER1(DRIVE_LFC_mats, 'DRIVE', DRIVE_batch_names, DRIVE_DEMETER1_mat_fname)

#### INTERMEDIATE STEP: Run DEMETER on the cloud with output of preprocessing step (see run_DEMETER1 for more details) ####

# Create "output" directory in corresponding DEMETER1 directory and save DEMETER1 output
DRIVE_DEMETER1_output_fname <- '~/Desktop/rnai_pipeline/DRIVE/DEMETER1/output/GeneSols.csv'
Achilles_DEMETER1_output_fname <- '~/Desktop/rnai_pipeline/Achilles/DEMETER1/output/GeneSolsCleaned.csv'

postprocess_DEMETER1_output(DRIVE_DEMETER1_output_fname, DRIVE_DEMETER1_scores_fname)
postprocess_DEMETER1_output(Achilles_DEMETER1_output_fname, Achilles_DEMETER1_scores_fname)


