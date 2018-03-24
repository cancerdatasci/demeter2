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

source('~/Documents/github/demeter2/scripts/setup_filesystem.R')


######################
# Shared Functions
######################

# Stores functions for annotating CCLE cancer cell lines
load('~/Desktop/other_models/Annotations.RData')

# Load LFC Matrices from Taiga, store as global variables
get_LFC_matrices <- function() {
  # Achilles
  a98k_LFC_mat <<- load.from.taiga(data.name='achilles-98k-repcollapsed-lfc-19ce', data.version=1)
  a55k1_LFC_mat <<- load.from.taiga(data.name='achilles-55k-batch1-repcollapsed-lfc-d708', data.version=1)
  a55k2_LFC_mat <<- load.from.taiga(data.name='achilles-55k-batch2-repcollapsed-lfc-bd7f', data.version=1)
  Achilles_LFC_mats <<- list(a98k_LFC_mat, a55k1_LFC_mat, a55k2_LFC_mat)
  # DRIVE
  BGPD_LFC_mat <<- load.from.taiga(data.name='drive-lfc-matrices-3867', data.version=10, data.file='BGPD_LFC_mat')
  poolA_LFC_mat <<- load.from.taiga(data.name='drive-lfc-matrices-3867', data.version=10, data.file='poolA_LFC_mat')
  poolB_LFC_mat <<- load.from.taiga(data.name='drive-lfc-matrices-3867', data.version=10, data.file='poolB_LFC_mat')
  DRIVE_LFC_mats <<- list(BGPD_LFC_mat, poolA_LFC_mat, poolB_LFC_mat)
  # Marcotte
  Marcotte_LFC_mat <- load.from.taiga(data.name='marcotte-demeter-z-scored-expanded-gene-sols-e905', data.version=10, data.file='log.effect.removed.bad')
  # Convert cleannames to CCLE IDs
  colnames(Marcotte_LFC_mat) <- ifelse(is.na(CleanCellLineName(colnames(Marcotte_LFC_mat))), 
                                       colnames(Marcotte_LFC_mat), 
                                       CleanCellLineName(colnames(Marcotte_LFC_mat)))
  Marcotte_LFC_mat <<- Marcotte_LFC_mat
}

# Find sequences that target too many genes (>10) and filter them out
remove_promiscuous_seqs <- function(map) {
  promiscuous_seqs <- map %>% 
    dplyr::select(seq, Gene = Entrez_ID) %>% 
    unique() %>% 
    dplyr::group_by(seq) %>% 
    dplyr::summarise(num_genes_targeted = n()) %>% 
    dplyr::filter(num_genes_targeted > 10) %$%
    seq
  
  clean_map <- map %>% dplyr::filter(!(seq %in% promiscuous_seqs))
  print(paste0('Removed ', nrow(map) - nrow(clean_map), ' promiscuous sequences from hairpin map'))
  return(clean_map)
}

# Remove outdated genes to create hairpin map
get_hairpin_map <- function() {
  map <- load.from.taiga(data.name='gpp-shrna-mapping-8759', data.version=2, data.file='shmap_19mer_noXLOC')
  map %<>% dplyr::select(
    Gene = `Gene Symbol`, 
    Entrez_ID = `Gene ID`, 
    seq = `Barcode Sequence`) %>% 
    dplyr::filter(!grepl('NO_CURRENT', Gene)) %>% # filter out any genes containing "NO_CURRENT..."
    remove_promiscuous_seqs()
  hairpin_map <<- map
}

# Make sure X doesn't get appended to beginning of column names
check_CCLE_colnames <- function(mat) {
  stopifnot(length(grep('^X[0-9]', colnames(mat))) == 0)
}

# Normalizes LFC matrix by CL and changes rownames to reflect the sequence, entrez id, and batch name
process_LFC <- function(LFC, datatype, batch_name, distinct_hp) {
  LFC %<>% scale(center = TRUE, scale = TRUE) # zscore normalize by CL
  cur <- hairpin_map %>% 
    dplyr::filter(seq %in% rownames(LFC)) %>% 
    dplyr::select(seq, Gene, Entrez_ID) %>% 
    unique()
  
  new_LFC <- LFC[cur$seq, ] # Split into separate rows for hairpins targeting each gene
  rownames(new_LFC) <- if (distinct_hp) {
    # hairpins from different batches should be treated as distinct
    paste0(cur$seq, '@', cur$Gene, '@', cur$Entrez_ID, '@', batch_name)
  } else {
    # hairpins from different batches should be treated as the same
    paste0(cur$seq, '@', cur$Gene, '@', cur$Entrez_ID)
  }
  
  long <- new_LFC %>% 
    as.matrix() %>% 
    melt(varnames = c('Gene', 'CL'), value.name = 'LFC')
  return(long)
}


##############
# Fit Models
##############

# Retrieve data for all models and hairpin map globally
get_LFC_matrices()
get_hairpin_map()


###################
# Gene Averaging
###################

preprocess_for_GA <- function(LFC) {
  long <- LFC %>% scale(center = TRUE, scale = TRUE) %>% # zscore normalize by CL
    reshape2::melt() %>% 
    set_colnames(c('seq', 'CL', 'LFC')) %>% 
    dplyr::left_join(hairpin_map, by = 'seq') %>% 
    dplyr::mutate(Gene_Entrez_ID = paste0(Gene, ' (', Entrez_ID, ')')) %>% 
    dplyr::select(CL, LFC, Gene_Entrez_ID)
  return(long)
}

fit_GA <- function(LFC, fname) {
  GA <- LFC %>%
    dplyr::group_by(CL, Gene_Entrez_ID) %>%
    dplyr::summarise(avg_score = mean(LFC, na.rm=T), num_shRNAs = n()) %>% 
    dplyr::filter(num_shRNAs > 1) %>% # ensure >1 hairpin/gene
    reshape2::acast(Gene_Entrez_ID ~ CL, value.var = 'avg_score') %>% 
    write.csv(fname)
}

GA_wrapper <- function(datatype, fname, LFC_mats) {
  LFC <- mapply(preprocess_for_GA, LFC = LFC_mats, SIMPLIFY = FALSE)
  LFC <- do.call(rbind, LFC)
  fit_GA(LFC, fname)
}

GA_wrapper('Achilles', paths$Achilles_GA_fname, Achilles_LFC_mats)
GA_wrapper('DRIVE', paths$DRIVE_GA_fname, DRIVE_LFC_mats)
GA_wrapper('Marcotte', paths$Marcotte_GA_fname, list(Marcotte_LFC_mat))

#### Achilles + DRIVE + Marcotte GA takes lots of memory to run. We ran this remotely.
#### Our pipeline involves uploading an intermeditate RDS file to Google Cloud, as shown below.

#### Our approach to running the GA model for the combined data

# GA_combination <- function(LFC_mats) {
#   LFC <- mapply(preprocess_for_GA, LFC = LFC_mats, SIMPLIFY = FALSE)
#   LFC <- do.call(rbind, LFC)
#   write_rds(LFC, paths$Combined_RDS_fname)
# }
# 
# GA_combination(list(a98k_LFC_mat, a55k1_LFC_mat, a55k2_LFC_mat, 
#                     BGPD_LFC_mat, poolA_LFC_mat, poolB_LFC_mat, 
#                     Marcotte_LFC_mat))
# 
# LFC <- read_rds(paths$Combined_RDS_fname)
# fit_GA(LFC, paths$Combined_GA_fname)


##########
# RSA
##########

# Create a file for each cell line, run RSA for each cell line, and stitch outputs together
preprocess_for_RSA <- function(LFC_mats, dir) {
  LFC <- LFC_mats %>% llply(function(mat) {
    long <- mat %>% 
      scale(center = TRUE, scale = TRUE) %>% # zscore normalize by CL
      reshape2::melt() %>% 
      set_colnames(c('seq', 'CL', 'Score'))
    return(long)
  })
  
  long <- do.call(rbind, LFC) %>% 
    dplyr::left_join(hairpin_map, by = 'seq') %>%
    dplyr::mutate(Entrez_ID_CL = paste0(Entrez_ID, '@', CL)) %>% # Use "@" as a delimiter
    dplyr::select(CL, Gene_ID = Entrez_ID_CL, Well_ID = seq, Score)
  
  long %>% d_ply(.(CL), function(df) {
    print(df$CL %>% unique())
    # Filter out genes targeted by only one hairpin
    good_genes <- df %>% 
      dplyr::count(Gene_ID) %>% 
      dplyr::filter(n > 1) %$% Gene_ID
    
    df %>% 
      dplyr::filter(Gene_ID %in% good_genes) %>% 
      dplyr::select(Gene_ID, Well_ID, Score) %>% 
      write.csv(paste0(dir, unique(df$CL), '.csv'), row.names = F)
  })
}

postprocess_RSA_output <- function(dir, fname) {
  files <- list.files(dir, full.names = T) %>% as.list()
  RSA <- files %>% llply(function(f) {
    df <- read.csv(f, check.names = F)
    new <- df %>% cbind(colsplit(as.character(df$Gene_ID),"@",c("Entrez_ID","CL"))) %>%
      dplyr::select(CL, Entrez_ID, Score = LogP) %>%
      dplyr::left_join(hairpin_map, by = 'Entrez_ID') %>%
      dplyr::mutate(Gene_Entrez_ID = paste0(Gene, ' (', Entrez_ID, ')')) %>%
      dplyr::select(CL, Gene_Entrez_ID, Score) %>% unique()
    return(new)
  })
  do.call(rbind, RSA) %>% 
    reshape2::acast(Gene_Entrez_ID ~ CL, value.var = 'Score') %>% 
    write.csv(fname)
}

preprocess_for_RSA(Achilles_LFC_mats, paths$Achilles_RSA_inputs_directory)
preprocess_for_RSA(DRIVE_LFC_mats, paths$DRIVE_RSA_inputs_directory)

#### INTERMEDIATE STEP: Use scripts in run_RSA to generate RSA scores for each cell line ####

postprocess_RSA_output(paths$Achilles_RSA_outputs_directory, paths$Achilles_RSA_scores_fname)
postprocess_RSA_output(paths$DRIVE_RSA_outputs_directory, paths$DRIVE_RSA_scores_fname)


###########
# ATARiS
###########

write_ATARiS_gct <- function(df, gct_fname) {
  gene_entrez_ids <- rownames(df) %>% 
    llply(function(str) {
      chunks <- strsplit(str, '@')[[1]]
      return(paste0(chunks[2], ' (', chunks[3], ')'))
    }) %>% unlist()
  rownames(df) <- gsub('@', '_', rownames(df)) # Convert back to "_" for ATARiS
  check_CCLE_colnames(df)
  dependr::write.gct(df, gct_fname, gene_entrez_ids)
}

preprocess_for_ATARiS <- function(LFC_mats, datatype, batch_names) {
  distinct_hp <- ifelse(datatype == 'Achilles', FALSE, TRUE)
  LFCs <- mapply(process_LFC, LFC = LFC_mats, datatype = datatype, batch_name = batch_names, distinct_hp, SIMPLIFY = FALSE)
  if (datatype == 'DRIVE') {
    Reduce(rbind, LFCs) %>% 
      reshape2::acast(Gene ~ CL, value.var = 'LFC') %>% 
      write_ATARiS_gct(paths$DRIVE_preATARiS_fname)
  } else if (datatype == 'Achilles') {
    # Achilles 98k batch
    LFCs[[1]] %>% 
      reshape2::acast(Gene ~ CL, value.var = 'LFC') %>% 
      write_ATARiS_gct(paths$Achilles_98k_preATARiS_fname)
    # Achilles 55k batches
    Reduce(rbind, LFCs[2:3]) %>% 
      reshape2::acast(Gene ~ CL, value.var = 'LFC') %>% 
      write_ATARiS_gct(paths$Achilles_55k_preATARiS_fname)
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
  long <- filtered_df %>% 
    as.matrix() %>% 
    melt(varnames = c('Gene', 'CL'), value.name = 'Score')
  return(long)
}

# Stitches output from 55k and 98k batches to form a single matrix
postprocess_ATARiS_output <- function(gct_files, fname, stitch) {
  if (stitch) {
    LFCs <- mapply(helper_filter, fname = gct_files, SIMPLIFY = FALSE)
    df <- Reduce(rbind, LFCs) %>% 
      reshape2::acast(Gene ~ CL, value.var = 'Score')
  } else {
    df <- helper_filter(gct_files) %>% 
      reshape2::acast(Gene ~ CL, value.var = 'Score')
  }
  check_CCLE_colnames(df)
  write.csv(df, fname)
}

preprocess_for_ATARiS(Achilles_LFC_mats, 'Achilles', Achilles_batch_names)
preprocess_for_ATARiS(DRIVE_LFC_mats, 'DRIVE', DRIVE_batch_names)

#### INTERMEDIATE STEP: Run ATARiS on Gene Pattern ####

# Save ATARiS output in corresponding ATARiS "output" directory
# Use identifiers Achilles_55k, Achilles_98k, and DRIVE
postprocess_ATARiS_output(
  paste0(paths$Achilles_ATARiS_output_directory, c('Achilles_55k.Gs.gct', 'Achilles_98k.Gs.gct')), 
  paths$Achilles_ATARiS_scores_fname, TRUE)

postprocess_ATARiS_output(
  paste0(paths$DRIVE_ATARiS_output_directory, 'DRIVE.Gs.gct'), 
  paths$DRIVE_ATARiS_scores_fname, FALSE)

############
# DEMETER 1 
############

# Create batch files for input to DEMETER1
make_batches_file <- function(LFC_mats, datatype) {
  batches <- data.frame()
  for (i in 1:length(LFC_mats)) {
    cur <- LFC_mats[[i]]
    batches %<>% rbind(data.frame(
      Name = colnames(cur),
      `DEMETER batch` = rep(i, ncol(cur)),
      check.names = FALSE))
  }
  if (datatype == 'Achilles') {
    batches %>% write.csv(paths$Achilles_DEMETER1_batches_fname, row.names = FALSE)
  } else if (datatype == 'DRIVE') {
    # All CLs in DRIVE data will be considered to be in the same batch
    batches %<>% 
      dplyr::mutate(`DEMETER batch` = 1) %>% 
      unique() %>% 
      write.csv(paths$DRIVE_DEMETER1_batches_fname, row.names = FALSE)
  }
}

write_DEMETER1_gct <- function(LFC, fname, datatype) {
  LFC %<>% reshape2::acast(Gene ~ CL, value.var = 'LFC')
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
    LFCs <- mapply(process_LFC, LFC = LFC_mats, datatype = datatype, batch_name = batch_names, distinct_hp, SIMPLIFY = FALSE)
    mapply(write_DEMETER1_gct, LFC = LFCs, fname = mat_fnames, datatype = datatype)
  } else if (datatype == 'DRIVE') {
    LFCs <- mapply(process_LFC, LFC = LFC_mats, datatype = datatype, batch_name = batch_names, distinct_hp, SIMPLIFY = FALSE)
    Reduce(rbind, LFCs) %>% 
      write_DEMETER1_gct(mat_fnames, datatype = datatype)
  }
}

# Convert rownames to Gene Symbol (Entrez ID)
postprocess_DEMETER1_output <- function(input_fname, output_fname) {
  df <- read.csv(input_fname, check.names = F, row.names = 1)
  # Parse DEMETER1 output. Multiple Entrez IDs that map to the same gene are separated by "&" in some rows.
  entrez_ids <- ifelse(grepl('&', rownames(df)), gsub(')', '', gsub('.+&', '', rownames(df))), rownames(df))
  map <- hairpin_map %>% 
    dplyr::select(Gene, Entrez_ID) %>% 
    unique() %>% 
    dplyr::filter(Entrez_ID %in% entrez_ids) %>% 
    dplyr::mutate(Gene_Entrez_ID = paste0(Gene, ' (', Entrez_ID, ')')) %>% 
    dplyr::slice(match(entrez_ids, Entrez_ID))
  rownames(df) <- map$Gene_Entrez_ID
  write.csv(df, output_fname)
}

preprocess_for_DEMETER1(Achilles_LFC_mats, 'Achilles', Achilles_batch_names, paths$Achilles_DEMETER1_mat_fnames)
preprocess_for_DEMETER1(DRIVE_LFC_mats, 'DRIVE', DRIVE_batch_names, paths$DRIVE_DEMETER1_mat_fname)

#### INTERMEDIATE STEP: Run DEMETER on the cloud with output of preprocessing step (see run_DEMETER1 for more details) ####

# Save DEMETER1 output in corresponding DEMETER1 "output" directory
Achilles_DEMETER1_output_fname <- paste0(paths$Achilles_DEMETER1_output_directory, 'GeneSols.csv')
DRIVE_DEMETER1_output_fname <- paste0(paths$DRIVE_DEMETER1_output_directory, 'GeneSols.csv')

postprocess_DEMETER1_output(Achilles_DEMETER1_output_fname, paths$Achilles_DEMETER1_scores_fname)
postprocess_DEMETER1_output(DRIVE_DEMETER1_output_fname, paths$DRIVE_DEMETER1_scores_fname)


