library(plyr); 
library(magrittr);  
library(reshape2)
library(tibble)
library(dplyr)
library(tidyr)
library(taigr);
library(ggplot2)
library(ggrepel)
library(dependr)
library(readr)
library(stringr)
library(cowplot)
library(cdsr)
source('~/CPDS/packages/demeter2_pub/scripts/benchmark_helpers.R')

min_hps_per_gene <- 4

true_gene_deps <- load.from.taiga(data.name='demeter2-achilles-5386', data.version=12, data.file='gene_means_proc')
rownames(true_gene_deps) <- str_match(rownames(true_gene_deps), '\\((.+)\\)')[,2]

sh_targets <- load.from.taiga(data.name='gpp-shrna-mapping-8759', data.version=6, data.file = 'shmap_19mer_noXLOC_Entrezonly')
hart_ess <- load.from.taiga(data.name='demeter2-pos-neg-controls-a5c6', data.version=1, data.file='hart_pos_controls')$Gene_ID
hart_non_ess <- load.from.taiga(data.name='demeter2-pos-neg-controls-a5c6', data.version=1, data.file='hart_neg_controls')$Gene_ID

ground_truth_dir <- '~/CPDS/data/d2_testing/simulation/'
ground_truth_file <- 'ground_truth_LFCs_0.csv'
noise_SD_file <- 'noise_SDs_0.csv'

ground_truth_LFCs <- jmmBase::load_matrix_fast(file.path(ground_truth_dir, ground_truth_file))
noise_SDs <- read.csv(file.path(ground_truth_dir, noise_SD_file), check.names = F, stringsAsFactors = F, row.names = 1) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = 'CCLE_ID')

sh_targets %<>% filter(`Barcode Sequence` %in% rownames(ground_truth_LFCs), !grepl('NO_CURRENT', `Gene ID`))

hps_per_gene <- sh_targets %>% 
  group_by(`Gene ID`) %>% 
  summarise(n = n())

poss_genes <- hps_per_gene %>% 
  filter(n >= min_hps_per_gene) %>% 
  .[['Gene ID']]
poss_CLs <- colnames(ground_truth_LFCs)

# n_CL_set <- c(2, 5, 10, 25, 50, 100, 200)
# n_genes_set <- c(Inf, Inf, Inf, Inf, Inf, Inf, Inf)
# n_iters_set <- c(50, 25, 20, 10, 10, 5, 5)
# out_file_path <- '~/CPDS/demeter2/results/GA_sim_results_n_CLs.csv'

# n_CL_set <- c(200, 200, 200, 200, 200)
# n_genes_set <- c(500, 1000, 3500, 7500, 16000)
# n_iters_set <- c(50, 25, 10, 5, 2)
out_file_path <- '~/CPDS/demeter2/results/GA_sim_results_n_genes.csv'

all_res <- ldply(seq(length(n_CL_set)), function(jj) {
  n_iters <- n_iters_set[jj]
  n_CLs <- n_CL_set[jj]
  n_genes <- n_genes_set[jj]
  print(jj)
  ldply(seq(n_iters), function(ii) {
    CL_samples <- sample(poss_CLs, size = n_CLs, replace = FALSE)
    CL_noise_samples <- noise_SDs[CL_samples, 'noise_SD'] %>% set_names(CL_samples)
    gene_samples <- sample(poss_genes, size = pmin(n_genes, length(poss_genes)), replace = FALSE)
    used_hps <- sh_targets %>% 
      filter(`Gene ID` %in% gene_samples) %>% 
      .[['Barcode Sequence']]
    sprintf('%d hps, %d genes, %d poss cons, %d neg cons',
            length(used_hps),
            pmin(n_genes, length(poss_genes)),
            length(intersect(hart_ess, gene_samples)), 
            length(intersect(hart_non_ess, gene_samples)))
    
    sampled_gt_LFCs <- ground_truth_LFCs[used_hps, CL_samples]
    noise <- laply(CL_samples, function(CL) {
      rnorm(nrow(sampled_gt_LFCs), mean = 0, sd = CL_noise_samples[CL])
    }) %>% t()
    sampled_LFCs <- sampled_gt_LFCs + noise
    
    # write.csv(sampled_LFCs, '~/CPDS/data/d2_testing/cur_sim/sim_LFCs.csv')
    
    LFCs_norm <- sampled_LFCs %>% 
      scale(center = T, scale = T)
    
    GS <- LFCs_norm %>% 
      melt() %>% 
      set_colnames(c('SEQ', 'CCLE_ID', 'LFC')) %>% 
      left_join(sh_targets %>% dplyr::select(SEQ = `Barcode Sequence`, `Gene ID`), by = 'SEQ') %>% 
      group_by(`Gene ID`, CCLE_ID) %>% 
      summarise(avg_LFC = mean(LFC, na.rm=T)) %>% 
      ungroup() %>% 
      acast(`Gene ID` ~ CCLE_ID, value.var = 'avg_LFC')
    
    #calc SSMD
    avg_GS <- rowMeans(GS, na.rm=T)
    gene_types <- get_gene_type(names(avg_GS), hart_ess, hart_non_ess)
    cur_NE_calls <- get_gene_calls_NE(gene_types)
    posneg_SSMD <- get_SSMD(avg_GS, cur_NE_calls)
    
    intersect_genes <- intersect(rownames(GS), rownames(true_gene_deps))
    intersect_CLs <- intersect(colnames(GS), colnames(true_gene_deps))
    
    ov_cor <- cor(as.vector(GS[intersect_genes, intersect_CLs]),
                  as.vector(true_gene_deps[intersect_genes, intersect_CLs]),
                  use = 'pairwise.complete.obs')
    
    true_ms <- scale(true_gene_deps[intersect_genes, intersect_CLs] %>% t(), center = T, scale = F) %>% t()
    fit_ms <- scale(GS[intersect_genes, intersect_CLs] %>% t(), center = T, scale = F) %>% t()
    ov_cor_ms <- cor(as.vector(true_ms),
                     as.vector(fit_ms),
                     use = 'pairwise.complete.obs')
    
    results <- data.frame(n_genes = n_genes, n_CLs = n_CLs,
                          ov_cor = ov_cor, ov_cor_ms = ov_cor_ms,
                          posneg_SSMD = posneg_SSMD, iter = ii)
    return(results)
  })
})

write.csv(all_res, out_file_path, row.names = F)


