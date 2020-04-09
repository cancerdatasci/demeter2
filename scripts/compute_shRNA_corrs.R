library(tidyverse)
library(plyr)
library(magrittr)
library(taigr)
doMC::registerDoMC(cores = 4)

##PARAMS
min_pairs <- 10 #minimum number of overlapping cell lines to compute LFC correlations between hairpin pairs


#shRNA to gene mapping file
sh_mapping <- load.from.taiga(data.name='gpp-shrna-mapping-8759', data.version=6, data.file='shmap_19mer_noXLOC_Entrezonly')


compute_hairpin_stats <- function(shRNA_LFCs, sh_mapping, D2_gene_means) {
  all_hps <- Reduce(union, lapply(shRNA_LFCs, rownames))
  all_CLs <- Reduce(union, lapply(shRNA_LFCs, colnames))
  
  #helper function for merging LFC matrices into a single aligned average
  add_LFC_data <- function(data_obj, new_mat) {
    uu <- !is.na(new_mat)
    data_obj$sum[rownames(new_mat), colnames(new_mat)][uu] %<>% 
      magrittr::add(new_mat[uu])
    data_obj$counts[rownames(new_mat), colnames(new_mat)][uu] %<>% 
      magrittr::add(1)
    return(data_obj)
  }
  
  data_obj <- list()
  data_obj$sum <- matrix(0, nrow = length(all_hps), ncol = length(all_CLs)) %>% 
    set_rownames(all_hps) %>% 
    set_colnames(all_CLs)
  data_obj$counts <- matrix(0, nrow = length(all_hps), ncol = length(all_CLs)) %>% 
    set_rownames(all_hps) %>% 
    set_colnames(all_CLs)
  for (mat in shRNA_LFCs) {
    data_obj <- add_LFC_data(data_obj, mat)
  }
  LFC_merged <- data_obj$sum %>% magrittr::divide_by(data_obj$counts)
  LFC_merged[data_obj$counts == 0] <- NA
  LFC_merged %<>% t()
  
  #mapping for relevant targeting shRNAs
  sh_targets <- sh_mapping %>% 
    filter(`Barcode Sequence` %in% all_hps, 
           !grepl('NO_CURRENT', `Gene ID`), 
           grepl('[0-9]', `Gene ID`))
  
  all_genes <- sh_targets %>% 
    dplyr::select(`Gene Symbol`, `Gene ID`) %>% 
    dplyr::distinct(`Gene ID`, .keep_all=T)
  
  res <- ddply(all_genes, 1, function(grow) {
    target_gene <- grow[['Gene ID']]
    cur_seqs <- sh_targets %>% 
      filter(`Gene ID` == target_gene) %>% 
      .[['Barcode Sequence']]
    
    n_pairs <- t(!is.na(LFC_merged[,cur_seqs])) %*% (!is.na(LFC_merged[,cur_seqs])) # same as count.pairwise(x,y) from psych package
    
    if (length(cur_seqs) >= 2) { #if there were at least 2 hairpins targeting the gene
      temp <- cor(LFC_merged[,cur_seqs], use = 'pairwise.complete.obs') #compute correlation matrix
      diag(temp) <- NA
      temp[n_pairs < min_pairs] <- NA #ignore any pairs with too few overlapping lines
      cor_stats <- adply(temp, 1, function(row) {
        data.frame(
          best_cor = max(row, na.rm=T),
          up_quant_cor = quantile(row, 0.75, na.rm=T)
        )
      }, .id = 'hp') %>% 
        mutate(best_cor = ifelse(is.infinite(best_cor), NA, best_cor))
    } else {
      cor_stats <- data.frame(hp = cur_seqs, best_cor = NA, up_quant_cor = NA)
    }
    
    #compute correlation between shRNA LFC profile and D2 gene score for the target gene
    cur_title <- paste0(grow[['Gene Symbol']], ' (', target_gene, ')')
    if (cur_title %in% rownames(D2_gene_means)) {
      cur_GS <- D2_gene_means[cur_title, rownames(LFC_merged)]
      cor_stats$D2_GE_cor <- cor(LFC_merged[,cur_seqs, drop=FALSE], cur_GS, use = 'pairwise.complete.obs')[,1]
    } else {
      cor_stats$D2_GE_cor <- NA
    }
    
    cor_stats %<>% mutate(`Gene ID` = grow[['Gene ID']],
                          `Gene Symbol` = grow[['Gene Symbol']])
    return(cor_stats)
  }, .parallel = TRUE)
  return(res)
}

## -----------DRIVE-------------------
DRIVE_D2 <- load.from.taiga(data.name='demeter2-drive-0591', data.version=10, data.file='gene_means_proc')
#load LFC data per hp, center and scale by screen
DRIVE_mat_BGPD <- load.from.taiga(data.name = 'drive-lfc-matrices-3867',
                                  data.version = 11,
                                  data.file = 'BGPD_LFC_mat') %>% scale(center = T, scale = T)
DRIVE_mat_poolA <- load.from.taiga(data.name = 'drive-lfc-matrices-3867',
                                   data.version = 11,
                                   data.file = 'poolA_LFC_mat') %>% scale(center = T, scale = T)
DRIVE_mat_poolB <- load.from.taiga(data.name = 'drive-lfc-matrices-3867',
                                   data.version = 11,
                                   data.file = 'poolB_LFC_mat') %>% scale(center = T, scale = T)
DRIVE_D2_hp <- load.from.taiga(data.name='demeter2-drive-0591', data.version=10, data.file='hp_data_comb') %>% as.data.frame() %>% rownames_to_column(var = 'hp')

DRIVE_hp_stats <- compute_hairpin_stats(list(DRIVE_mat_BGPD, DRIVE_mat_poolA, DRIVE_mat_poolB), sh_mapping, DRIVE_D2)

all_hps <- Reduce(union, lapply(list(DRIVE_mat_BGPD, DRIVE_mat_poolA, DRIVE_mat_poolB), rownames))
DRIVE_hp_stats %<>% left_join(sh_mapping %>% 
                                dplyr::select(hp = `Barcode Sequence`, `Gene ID`) %>% 
                                filter(hp %in% all_hps) %>% 
                                dplyr::group_by(hp) %>% 
                                dplyr::summarise(n_targs = n()), by = 'hp')
DRIVE_hp_stats %<>% left_join(DRIVE_D2_hp %>% dplyr::select(hp, Geff, Seff), by = 'hp')
write_csv(DRIVE_hp_stats, '~/CPDS/data/D2_figshare/DRIVE_hp_quality_metrics.csv')



## -----------Achilles-------------------
Ach_D2 <- load.from.taiga(data.name='demeter2-achilles-5386', data.version=11, data.file='gene_means_proc')

#load 3 batches of shRNA LFC data, center and scale per screen
Ach_mat_98k <- load.from.taiga(data.name = 'achilles-98k-repcollapsed-lfc-19ce', 
                               data.version = 1) %>% scale(center = T, scale = T)
Ach_mat_55k1 <- load.from.taiga(data.name = 'achilles-55k-batch1-repcollapsed-lfc-d708', 
                                data.version = 1) %>% scale(center = T, scale = T)
Ach_mat_55k2 <- load.from.taiga(data.name = 'achilles-55k-batch2-repcollapsed-lfc-bd7f', 
                                data.version = 1) %>% scale(center = T, scale = T)
Ach_D2_hp <- load.from.taiga(data.name='demeter2-achilles-5386', data.version=11, data.file='hp_data_comb') %>% as.data.frame() %>% rownames_to_column(var = 'hp')

Ach_hp_stats <- compute_hairpin_stats(list(Ach_mat_98k, Ach_mat_55k1, Ach_mat_55k2), sh_mapping, Ach_D2)

all_hps <- Reduce(union, lapply(list(Ach_mat_98k, Ach_mat_55k1, Ach_mat_55k2), rownames))
Ach_hp_stats %<>% left_join(sh_mapping %>% 
                              dplyr::select(hp = `Barcode Sequence`, `Gene ID`) %>% 
                              filter(hp %in% all_hps) %>% 
                              dplyr::group_by(hp) %>% 
                              dplyr::summarise(n_targs = n()), by = 'hp')
Ach_hp_stats %<>% left_join(Ach_D2_hp %>% dplyr::select(hp, Geff, Seff), by = 'hp')
write_csv(Ach_hp_stats, '~/CPDS/data/D2_figshare/Ach_hp_quality_metrics.csv')