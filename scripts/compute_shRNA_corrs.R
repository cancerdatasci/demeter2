library(knitr); 
library(plyr); 
library(magrittr); 
library(reshape2)
library(stringr)
library(tibble)
library(dependr)
library(dplyr)
library(taigr)
library(jmmBase)
library(ggplot2)
library(doMC)
doMC::registerDoMC(cores = 4)

sh_targets <- load.from.taiga(data.name = 'gpp-shrna-mapping-8759', 
                              data.version = 2, 
                              data.file = 'shmap_19mer_noXLOC') 
DRIVE_D2 <- load.from.taiga(data.name='demeter2-drive-0591', data.version=8, data.file='gene_means_proc')
DRIVE_used_genes <- str_match(rownames(DRIVE_D2), '\\((.+)\\)')[,2]

DRIVE_mat_BGPD <- load.from.taiga(data.name = 'drive-lfc-matrices-3867',
                               data.version = 9,
                               data.file = 'BGPD_LFC_mat') %>% scale(center = T, scale = T)
DRIVE_mat_poolA <- load.from.taiga(data.name = 'drive-lfc-matrices-3867',
                                  data.version = 9,
                                  data.file = 'poolA_LFC_mat') %>% scale(center = T, scale = T)
DRIVE_mat_poolB <- load.from.taiga(data.name = 'drive-lfc-matrices-3867',
                                   data.version = 9,
                                   data.file = 'poolB_LFC_mat') %>% scale(center = T, scale = T)

all_hps <- rownames(DRIVE_mat_BGPD) %>% 
  union(rownames(DRIVE_mat_poolA)) %>% 
  union(rownames(DRIVE_mat_poolB))
all_CLs <- colnames(DRIVE_mat_BGPD) %>% 
  union(colnames(DRIVE_mat_poolA)) %>% 
  union(colnames(DRIVE_mat_poolB))

#create merged dataframe (average LFC for same cell lines hairpins where repeating)
DRIVE_merged <- matrix(0, nrow = length(all_hps), ncol = length(all_CLs)) %>% 
  set_rownames(all_hps) %>% 
  set_colnames(all_CLs)
DRIVE_counts <- matrix(0, nrow = length(all_hps), ncol = length(all_CLs)) %>% 
  set_rownames(all_hps) %>% 
  set_colnames(all_CLs)
uu <- !is.na(DRIVE_mat_BGPD)
DRIVE_merged[rownames(DRIVE_mat_BGPD), colnames(DRIVE_mat_BGPD)][uu] %<>% 
  magrittr::add(DRIVE_mat_BGPD[uu])
DRIVE_counts[rownames(DRIVE_mat_BGPD), colnames(DRIVE_mat_BGPD)][uu] %<>% 
  magrittr::add(1)

uu <- !is.na(DRIVE_mat_poolA)
DRIVE_merged[rownames(DRIVE_mat_poolA), colnames(DRIVE_mat_poolA)][uu] %<>% 
  magrittr::add(DRIVE_mat_poolA[uu])
DRIVE_counts[rownames(DRIVE_mat_poolA), colnames(DRIVE_mat_poolA)][uu] %<>% 
  magrittr::add(1)

uu <- !is.na(DRIVE_mat_poolB)
DRIVE_merged[rownames(DRIVE_mat_poolB), colnames(DRIVE_mat_poolB)][uu] %<>% 
  magrittr::add(DRIVE_mat_poolB[uu])
DRIVE_counts[rownames(DRIVE_mat_poolB), colnames(DRIVE_mat_poolB)][uu] %<>% 
  magrittr::add(1)

DRIVE_merged %<>% magrittr::divide_by(DRIVE_counts)
DRIVE_merged[DRIVE_counts == 0] <- NA
DRIVE_merged %<>% t()

sh_targets %<>% 
  filter(`Barcode Sequence` %in% all_hps, 
         !grepl('NO_CURRENT', `Gene ID`), 
         `Gene ID` %in% DRIVE_used_genes,
          grepl('[0-9]', `Gene ID`))

all_genes <- sh_targets %>% 
  dplyr::select(`Gene Symbol`, `Gene ID`) %>% 
  dplyr::distinct(`Gene ID`, .keep_all=T)

res <- ddply(all_genes, 1, function(grow) {
  target_gene <- grow[['Gene ID']]
  cur_seqs <- sh_targets %>% 
    filter(`Gene ID` == target_gene) %>% 
    .[['Barcode Sequence']]
  
  temp <- cor(DRIVE_merged[,cur_seqs], use = 'pairwise.complete.obs')
  diag(temp) <- NA
  cor_stats <- adply(temp, 1, function(row) {
    data.frame(
      best_cor = max(row, na.rm=T),
      up_quant_cor = quantile(row, 0.75, na.rm=T)
    )
  }, .id = 'hp')
  cor_stats %<>% mutate(`Gene ID` = grow[['Gene ID']],
                        `Gene Symbol` = grow[['Gene Symbol']])
  return(cor_stats)
}, .parallel = TRUE)

other_data <- read.csv('~/CPDS/demeter2/kube_results/DRIVE_final/1/hp_R2_post.csv', stringsAsFactors = F, check.names = F)
comb <- full_join(res, other_data, by = 'hp')


######## Now for Achilles data
sh_targets <- load.from.taiga(data.name = 'gpp-shrna-mapping-8759', 
                              data.version = 2, 
                              data.file = 'shmap_19mer_noXLOC') 
Ach_D2 <- load.from.taiga(data.name='demeter2-achilles-5386', data.version=10, data.file='gene_means_proc')
Ach_used_genes <- str_match(rownames(Ach_D2), '\\((.+)\\)')[,2]

Ach_mat_98k <- load.from.taiga(data.name = 'achilles-98k-repcollapsed-lfc-19ce', 
                            data.version = 1) %>% scale(center = T, scale = T)
Ach_mat_55k1 <- load.from.taiga(data.name = 'achilles-55k-batch1-repcollapsed-lfc-d708', 
                                data.version = 1) %>% scale(center = T, scale = T)
Ach_mat_55k2 <- load.from.taiga(data.name = 'achilles-55k-batch2-repcollapsed-lfc-bd7f', 
                                data.version = 1) %>% scale(center = T, scale = T)

all_hps <- rownames(Ach_mat_98k) %>% 
  union(rownames(Ach_mat_55k1)) %>% 
  union(rownames(Ach_mat_55k2))
all_CLs <- colnames(Ach_mat_98k) %>% 
  union(colnames(Ach_mat_55k1)) %>% 
  union(colnames(Ach_mat_55k2))

#create merged dataframe (average LFC for same cell lines hairpins where repeating)
Ach_merged <- matrix(0, nrow = length(all_hps), ncol = length(all_CLs)) %>% 
  set_rownames(all_hps) %>% 
  set_colnames(all_CLs)
Ach_counts <- matrix(0, nrow = length(all_hps), ncol = length(all_CLs)) %>% 
  set_rownames(all_hps) %>% 
  set_colnames(all_CLs)
uu <- !is.na(Ach_mat_98k)
Ach_merged[rownames(Ach_mat_98k), colnames(Ach_mat_98k)][uu] %<>% 
  magrittr::add(Ach_mat_98k[uu])
Ach_counts[rownames(Ach_mat_98k), colnames(Ach_mat_98k)][uu] %<>% 
  magrittr::add(1)

uu <- !is.na(Ach_mat_55k1)
Ach_merged[rownames(Ach_mat_55k1), colnames(Ach_mat_55k1)][uu] %<>% 
  magrittr::add(Ach_mat_55k1[uu])
Ach_counts[rownames(Ach_mat_55k1), colnames(Ach_mat_55k1)][uu] %<>% 
  magrittr::add(1)

uu <- !is.na(Ach_mat_55k2)
Ach_merged[rownames(Ach_mat_55k2), colnames(Ach_mat_55k2)][uu] %<>% 
  magrittr::add(Ach_mat_55k2[uu])
Ach_counts[rownames(Ach_mat_55k2), colnames(Ach_mat_55k2)][uu] %<>% 
  magrittr::add(1)

Ach_merged %<>% magrittr::divide_by(Ach_counts)
Ach_merged[Ach_counts == 0] <- NA
Ach_merged %<>% t()

sh_targets %<>% 
  filter(`Barcode Sequence` %in% all_hps, 
         !grepl('NO_CURRENT', `Gene ID`), 
         `Gene ID` %in% Ach_used_genes,
         grepl('[0-9]', `Gene ID`))

all_genes <- sh_targets %>% 
  dplyr::select(`Gene Symbol`, `Gene ID`) %>% 
  dplyr::distinct(`Gene ID`, .keep_all=T)

res <- ddply(all_genes, 1, function(grow) {
  target_gene <- grow[['Gene ID']]
  cur_seqs <- sh_targets %>% 
    filter(`Gene ID` == target_gene) %>% 
    .[['Barcode Sequence']]
  
  temp <- cor(Ach_merged[,cur_seqs], use = 'pairwise.complete.obs')
  diag(temp) <- NA
  cor_stats <- adply(temp, 1, function(row) {
    data.frame(
      best_cor = max(row, na.rm=T),
      up_quant_cor = quantile(row, 0.75, na.rm=T)
    )
  }, .id = 'hp')
  cor_stats %<>% mutate(`Gene ID` = grow[['Gene ID']],
                        `Gene Symbol` = grow[['Gene Symbol']])
  return(cor_stats)
}, .parallel = TRUE)

other_data <- read.csv('~/CPDS/demeter2/kube_results/Ach_final/1/hp_R2_post.csv', stringsAsFactors = F, check.names = F)
comb <- full_join(res, other_data, by = 'hp')


