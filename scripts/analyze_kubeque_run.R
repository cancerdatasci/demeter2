library(plyr); library(magrittr);  
library(reshape2)
library(tibble)
library(dplyr)
library(tidyr)
library(taigr);
library(jsonlite)
library(ggplot2)
library(ggrepel)
library(dependr)
library(readr)
library(stringr)
library(cowplot)
library(weights)
library(psych)
library(knitr)
library(cdsr)
library(PRROC)
library(GGally)
library(limma)
library(RColorBrewer)
library(latex2exp)

source('~/CPDS/demeter2/scripts/benchmark_helpers.R')

results_dir <- '~/CPDS/demeter2/kube_results/Ach_rel_scan'
fig_dir <- '~/CPDS/demeter2/results/rev_figures'
# results_dir <- '/Volumes/xchip_cga_home/jmmcfarl/demeter2_models/Ach_rel_scan/'

CERES = load.from.taiga(
  data.name = 'avana-public-tentative-18q1-92d9',
  data.version = 5,
  data.file = 'gene_effect')
colnames(CERES) <- str_match(colnames(CERES), '\\((.+)\\)')[,2]

feature_datasets <- list(
  GE = list(taiga_name = 'ccle-rnaseq-logrpkm-protein-coding-2700',
            taiga_version = 1,
            transpose = T),
  CN = list(taiga_name = 'depmap-wes-cn-data-97cc', 
            taiga_version = 5, 
            data.file = 'public_18q1_gene_cn')
  # MUT_HOT = list(taiga_name='hotspot-924f', taiga_version=1)
)
feature_dnames <- names(feature_datasets) %>% set_names(names(feature_datasets))
feature_data <- load_all_feature_data(feature_dnames, feature_datasets)

hart_ess <- read_tsv('~/CPDS/demeter2/data/pos_con_IDs.tsv')$pos_con_IDs
hart_non_ess <- read_tsv('~/CPDS/demeter2/data/neg_con_IDs.tsv')$neg_con_IDs

CRISPR_gene_avgs <- data.frame(Gene = colnames(CERES), 
                               CRISPR_avgs = colMeans(CERES, na.rm=T))

sample_info <- load.from.taiga(data.name='cell-lines-found-by-dataset-b427', data.version=7)
colnames(sample_info)[1] <- 'CCLE_ID'
sample_info %<>% dplyr::select(-lineage) %>% 
  dplyr::rename(lineage = `Lineage (for CDS)`) %>% 
  mutate(lineage = ifelse(lineage == '', NA, lineage))

all_res_df <- data.frame()
run_dirs <- list.dirs(results_dir, recursive = FALSE)

for (cur_dir in run_dirs) {
  print(cur_dir)
  # cur_dir <- run_dirs[1]
  cur_d2_mod <- load_demeter2_model(cur_dir, use_bayes = TRUE)
  
  cur_pos_cons <- intersect(hart_ess, colnames(cur_d2_mod$gene_scores))
  cur_neg_cons <- intersect(hart_non_ess, colnames(cur_d2_mod$gene_scores))
  pos_median <- median(cur_d2_mod$gene_scores[, cur_pos_cons], na.rm=T)
  neg_median <- median(cur_d2_mod$gene_scores[, cur_neg_cons], na.rm=T)
  if (!is.null(cur_d2_mod$gene_SDs)) {
    cur_d2_mod$gene_SDs <- cur_d2_mod$gene_SDs / (neg_median - pos_median)
  }
  cur_d2_mod$gene_scores <- (cur_d2_mod$gene_scores - neg_median) / (neg_median - pos_median)
  
  cur_gene_avgs <- sapply(seq(ncol(cur_d2_mod$gene_scores)), function(colind) {
    wtd.mean(cur_d2_mod$gene_scores[, colind],
             1/cur_d2_mod$gene_SDs[, colind]^2)
  })
  gene_vars <- colMeans(cur_d2_mod$gene_SDs^2, na.rm=T)
  
  gene_avg_df <- data.frame(Gene = colnames(cur_d2_mod$gene_scores),
                            avg_score = cur_gene_avgs,
                            avg_var = gene_vars) %>% 
    left_join(CRISPR_gene_avgs, by = 'Gene') %>% 
    mutate(gene_type = get_gene_type(Gene, hart_ess, hart_non_ess))
  
  
  #calculate pos-neg sep for gene averages
  # cur_SSMD_rob <- with(gene_avg_df %>% filter(gene_type %in% c('essential', 'non_essential')), 
  # get_SSMD(avg_score, gene_type == 'essential', robust = TRUE))
  cur_SSMD <- with(gene_avg_df %>% filter(gene_type %in% c('essential', 'non_essential')), 
                   get_SSMD(avg_score, gene_type == 'essential', weights = 1/avg_var, robust = FALSE))
  
  #CRISPR agreement for gene averages
  CRISPR_cor <- with(gene_avg_df, wtd.cors(CRISPR_avgs, avg_score, weight = 1/avg_var))[1,1]
  
  #calculate pos-neg sep for each individual cl
  
  #fraction of variance accounted for by within- vs between-gene variance
  within_var <- mean(apply(cur_d2_mod$gene_scores, 2, var, na.rm=T), na.rm=T)
  between_var <- var(colMeans(cur_d2_mod$gene_scores, na.rm=T), na.rm=T)
  within_gene_vf <- within_var/(within_var+between_var)
  
  #AGO2 correlation for pos-cons
  vector <- feature_data$GE[, '27161'] %>% set_names(rownames(feature_data$GE))
  common_CLs <- intersect(rownames(cur_d2_mod$gene_scores), names(vector))
  cur_pos_cons <- intersect(hart_ess, colnames(cur_d2_mod$gene_scores))
  CL_avg_vars <- rowMeans(cur_d2_mod$gene_SDs^2, na.rm=T)
  AGO2_cors <- wtd.cors(cur_d2_mod$gene_scores[common_CLs, cur_pos_cons], 
                        feature_data$GE[common_CLs, '27161', drop = FALSE], 
                        weight = 1/CL_avg_vars[common_CLs])
  avg_AGO2_cor <- mean(AGO2_cors, na.rm=T)
  
  #CRISPR agreement for POSCONS
  common_CLs <- intersect(rownames(cur_d2_mod$gene_scores), rownames(CERES))
  common_pos_cons <- intersect(cur_pos_cons, colnames(CERES))
  CRISPR_poscon_cors <- sapply(common_pos_cons, function(cur_gene) {
    wtd.cors(cur_d2_mod$gene_scores[common_CLs, cur_gene], 
             CERES[common_CLs, cur_gene], 
             weight = 1/CL_avg_vars[common_CLs])
  })
  avg_CRISPR_cor <- mean(CRISPR_poscon_cors, na.rm=T)
  
  #self-CN cor for POSCONS
  common_CLs <- intersect(rownames(cur_d2_mod$gene_scores), rownames(feature_data$CN))
  common_pos_cons <- intersect(cur_pos_cons, colnames(feature_data$CN))
  CN_poscon_cors <- sapply(common_pos_cons, function(cur_gene) {
    wtd.cors(cur_d2_mod$gene_scores[common_CLs, cur_gene], 
             feature_data$CN[common_CLs, cur_gene], 
             weight = 1/CL_avg_vars[common_CLs])
    
  })
  avg_CN_cor <- mean(CN_poscon_cors, na.rm=T)
  
  # #curated set separation
  # ##KRAS
  # df <- data.frame(CCLE_ID = rownames(cur_d2_mod$gene_scores),
  #                  response = cur_d2_mod$gene_scores[, '3845']) %>% 
  #     inner_join(
  #         data.frame(CCLE_ID = rownames(feature_data$MUT_HOT),
  #                   variable = feature_data$MUT_HOT[, '3845'] > 0), 
  #         by = 'CCLE_ID')
  # 
  # m <- with(df, t.test(response[variable], response[!variable]))
  # KRAS_MUT_mdiff <- diff(m$estimate)
  # 
  # ##BRAF
  # df <- data.frame(CCLE_ID = rownames(cur_d2_mod$gene_scores),
  #                  response = cur_d2_mod$gene_scores[, '673']) %>% 
  #     inner_join(
  #         data.frame(CCLE_ID = rownames(feature_data$MUT_HOT),
  #                    variable = feature_data$MUT_HOT[, '673'] > 0), 
  #         by = 'CCLE_ID')
  # 
  # m <- with(df, t.test(response[variable], response[!variable]))
  # BRAF_MUT_mdiff <- diff(m$estimate)
  # 
  # ##SOX10
  # df <- data.frame(CCLE_ID = rownames(cur_d2_mod$gene_scores),
  #                  response = cur_d2_mod$gene_scores[, '6663']) %>% 
  #     inner_join(sample_info, 
  #         by = 'CCLE_ID') %>% 
  #     mutate(variable = lineage == 'skin')
  # m <- with(df, t.test(response[variable], response[!variable]))
  # SOX10_skin_mdiff <- diff(m$estimate)
  
  cur_df <- data.frame(mod_path = cur_dir,
                       SSMD = cur_SSMD,
                       # rob_SSMD = cur_SSMD_rob,
                       CRISPR_cor = CRISPR_cor,
                       within_gene_vf = within_gene_vf,
                       avg_AGO2_cor = avg_AGO2_cor,
                       avg_CRISPR_cor = avg_CRISPR_cor,
                       avg_CN_cor = avg_CN_cor
                       # KRAS_MUT_mdiff = KRAS_MUT_mdiff,
                       # BRAF_MUT_mdiff = BRAF_MUT_mdiff,
                       # SOX10_skin_mdiff = SOX10_skin_mdiff
  ) %>% 
    cbind(cur_d2_mod$other_dat$reg_params %>% as.data.frame())
  
  #R2?
  if (!is.null(cur_d2_mod$other_dat$R2_vals)) {
    R2_test <- cur_d2_mod$other_dat$R2_vals$test %>% unlist()
    R2_test <- R2_test[length(R2_test)]
    
    cur_df$R2_test <- R2_test
    
    R2_train <- cur_d2_mod$other_dat$R2_vals$train %>% unlist()
    R2_train <- R2_train[length(R2_train)]
    cur_df$R2_train <- R2_train
  }
  
  all_res_df %<>% rbind(cur_df)
}
write_rds(all_res_df, '~/CPDS/demeter2/results/new_kube_res.rds')


# ggplot(all_res_df, aes(log10(rel_seed_l2_lambda), log10(rel_gene_l2_lambda), fill = SSMD)) + 
#   geom_raster()
# ggsave(file.path(fig_dir, 'rel_lambdas_vs_SSMD.png'), width = 6, height = 6)

# ggplot(all_res_df, aes(log10(rel_seed_l2_lambda), log10(rel_gene_l2_lambda), fill = rob_SSMD)) + geom_raster()
# ggplot(all_res_df, aes(log10(rel_seed_l2_lambda), log10(rel_gene_l2_lambda), fill = CRISPR_cor)) + 
#   geom_raster()

ggplot(all_res_df, aes(log10(rel_seed_l2_lambda), log10(rel_gene_l2_lambda), fill = within_gene_vf)) + 
  geom_raster() +
  xlab(TeX('$log10(\\lambda_b)$')) +
  ylab(TeX('$log10(\\lambda_g)$')) +
  guides(fill = guide_colorbar(title = 'Within-gene\nvariance fraction')) +
  theme_Publication() +
  theme(legend.key.width = unit(0.25, 'inches'), legend.text = element_text(size = 10))
ggsave(file.path(fig_dir, 'rel_lambdas_vs_var_fraction.png'), width = 3.5, height = 3.5)

# ggplot(all_res_df, aes(log10(rel_seed_l2_lambda), log10(rel_gene_l2_lambda), fill = KRAS_MUT_mdiff)) + geom_raster()
# ggplot(all_res_df, aes(log10(rel_seed_l2_lambda), log10(rel_gene_l2_lambda), fill = BRAF_MUT_mdiff)) + geom_raster()
# ggplot(all_res_df, aes(log10(rel_seed_l2_lambda), log10(rel_gene_l2_lambda), fill = SOX10_skin_mdiff)) + geom_raster()

# ggplot(all_res_df, aes(log10(rel_seed_l2_lambda), log10(rel_gene_l2_lambda), fill = avg_CN_cor)) + 
#   geom_raster()

ggplot(all_res_df, aes(log10(rel_seed_l2_lambda), log10(rel_gene_l2_lambda), fill = avg_CRISPR_cor)) +
  geom_raster() +
  xlab(TeX('$log10(\\lambda_b)$')) +
  ylab(TeX('$log10(\\lambda_g)$')) +
  scale_fill_gradient(limits = c(0.215, 0.24), breaks = c(0.22, 0.24)) +
  guides(fill = guide_colorbar(title = 'Corr. with CRISPR')) +
  theme_Publication() +
  theme(legend.key.width = unit(0.25, 'inches'), legend.text = element_text(size = 10))
ggsave(file.path(fig_dir, 'rel_lambdas_vs_CRISPR_cor.png'), width = 3.5, height = 3.5)

# ggplot(all_res_df, aes(log10(rel_seed_l2_lambda), log10(rel_gene_l2_lambda), fill = avg_AGO2_cor)) + 
#   geom_raster()

ggplot(all_res_df, aes(log10(rel_seed_l2_lambda), log10(rel_gene_l2_lambda), fill = R2_test)) + 
  geom_raster() +
  xlab(TeX('$log10(\\lambda_b)$')) +
  ylab(TeX('$log10(\\lambda_g)$')) +
  guides(fill = guide_colorbar(title = TeX('Test $R^2$'))) +
  scale_fill_gradient(breaks = c(0.74, 0.78)) +
  theme_Publication() +
  theme(legend.key.width = unit(0.25, 'inches'), legend.text = element_text(size = 10))
ggsave(file.path(fig_dir, 'rel_lambdas_vs_R2test.png'), width = 3.5, height = 3.5)

ggplot(all_res_df, aes(log10(rel_seed_l2_lambda), log10(rel_gene_l2_lambda), fill = R2_train)) + 
  geom_raster() +
  xlab(TeX('$log10(\\lambda_b)$')) +
  ylab(TeX('$log10(\\lambda_g)$')) +
  guides(fill = guide_colorbar(title = TeX('Train $R^2$'))) +
  theme_Publication() +
  theme(legend.key.width = unit(0.25, 'inches'), legend.text = element_text(size = 10))
ggsave(file.path(fig_dir, 'rel_lambdas_vs_R2train.png'), width = 3.5, height = 3.5)
