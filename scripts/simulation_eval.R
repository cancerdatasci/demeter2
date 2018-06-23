library(plyr); 
library(magrittr);  
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
library(cdsr)

source('~/CPDS/packages/demeter2_pub/scripts/benchmark_helpers.R')
normalization <- 'pos_neg' #type of normalization [pos_neg, global_z] (absolute GS only)
include_gene_families <- FALSE #whether to include gene-families in analysis
use_bayes <- TRUE

min_avg_Geff <- 0.2 #minimum average Geff across hairpins to include a gene
min_sum_Geff <- 1.5 #minimum sum Geff across hairpins to include a gene
min_hps_per_gene <- 3 #minimum hairpins per gene for inclusion

fig_dir <- '~/CPDS/demeter2/results/rev_figures'
fig_width <- 3.5
fig_height <- 3

base_data_dir <- '~/CPDS/demeter2/kube_results/'
dep_datasets <- list(
  fit = list(mod_type = 'D2', 
             # path = file.path(base_data_dir, 'sim_fit_noeff_sp/1/'),
             path = file.path(base_data_dir, 'sim_fit/2/'),
             dataset = 'Ach', 
             type = 'abs'),
  true = list(mod_type = 'D2', 
              path = file.path(base_data_dir, 'Ach_final/1/'), 
              dataset = 'Ach', type = 'abs')
)
hart_ess <- load.from.taiga(data.name='demeter2-pos-neg-controls-a5c6', data.version=1, data.file='hart_pos_controls')$Gene_ID
hart_non_ess <- load.from.taiga(data.name='demeter2-pos-neg-controls-a5c6', data.version=1, data.file='hart_neg_controls')$Gene_ID
dep_dnames <- names(dep_datasets) %>% set_names(names(dep_datasets))
dep_data <- load_all_dep_data(dep_dnames, dep_datasets, include_gene_families = include_gene_families, use_bayes = use_bayes, name_map = NULL)

sh_targets <- load.from.taiga(data.name='gpp-shrna-mapping-8759', data.version=6, data.file = 'shmap_19mer_noXLOC_Entrezonly')

#remove genes with poor quality reagenets
avg_Geffs <- data.frame()
for (cur_dset in names(dep_datasets)) {
  cur_avg_Geff <- read.csv(file.path(dep_data[[cur_dset]]$path, 'hp_data.csv'), stringsAsFactors = F, check.names=F) %>% 
    dplyr::select(Geff, Seff, hp) %>% 
    left_join(sh_targets %>% dplyr::select(hp = `Barcode Sequence`, Gene_ID = `Gene ID`), by = 'hp') %>% 
    group_by(Gene_ID) %>% 
    summarise(mean_Geff = mean(Geff), sum_Geff = sum(Geff), n_hps = n())
  
  #eliminate any genes with either low mean Geff, low sum Geff, OR too few hps
  bad_genes <- cur_avg_Geff %>% 
    filter(Gene_ID %in% colnames(dep_data[[cur_dset]]$gene_scores),
           mean_Geff < min_avg_Geff | sum_Geff < min_sum_Geff | n_hps < min_hps_per_gene) %>% 
    .[['Gene_ID']]
  print(sprintf('%s: removing %d/%d bad genes', cur_dset, length(bad_genes), ncol(dep_data[[cur_dset]]$gene_scores)))
  dep_data[[cur_dset]]$gene_scores[, bad_genes] <- NA
  dep_data[[cur_dset]]$gene_score_SDs[, bad_genes] <- NA
  
  avg_Geffs %<>% rbind(cur_avg_Geff %>% mutate(dset = cur_dset))
}
avg_Geffs %<>% filter(!is.na(Gene_ID))

#normalize each dataset
cur_gene_avgs <- get_gene_avgs(dep_data, names(dep_data), hart_ess, hart_non_ess, use_bayes = use_bayes)
for (cur_dset in names(dep_datasets)) {
  cur_GA <- cur_gene_avgs %>% filter(dset == cur_dset)
  pos_median <- cur_GA %>% filter(gene_type == 'essential') %>% .[['avg_score']] %>% median(na.rm=T)
  neg_median <- cur_GA %>% filter(gene_type == 'non_essential') %>% .[['avg_score']] %>% median(na.rm=T)
  if (!is.null(dep_data[[cur_dset]]$gene_score_SDs)) {
    dep_data[[cur_dset]]$gene_score_SDs <- dep_data[[cur_dset]]$gene_score_SDs / (neg_median - pos_median)
  }
  dep_data[[cur_dset]]$gene_scores <- (dep_data[[cur_dset]]$gene_scores - neg_median) / (neg_median - pos_median)
}


#compare overall gene averages
gene_avgs <- get_gene_avgs(dep_data, names(dep_datasets), hart_ess, hart_non_ess, use_bayes = use_bayes) %>% 
  mutate(Gene = as.character(Gene))
gene_avgs_wide <- spread(gene_avgs, dset, avg_score)
plot_colorByDensity(gene_avgs_wide, 'true', 'fit') +
  geom_abline() +
  xlab('True avg dependency score') +
  ylab('Fit avg dependency score') +
  coord_equal() + 
  theme_Publication()
ggsave(file.path(fig_dir, 'gene_avg_sim_compare.png'), width = fig_width, height = fig_height, dpi = 350)

cur_SSMD <- with(gene_avgs, ssmdMM(avg_score[gene_type == 'non_essential'],
                                   avg_score[gene_type == 'essential']))

#CL scaling parameters
true <- read.csv(file.path(dep_datasets$true$path, 'CL_batch_data.csv'), stringsAsFactors = F, check.names = F) 
fit <- read.csv(file.path(dep_datasets$fit$path, 'CL_batch_data.csv'), stringsAsFactors = F, check.names = F) 
CL_data <- full_join(true, fit, by = "CCLE_ID", suffix = c('_true', '_fit'))
true <- read.csv(file.path(dep_datasets$true$path, 'CL_data.csv'), stringsAsFactors = F, check.names = F) 
fit <- read.csv(file.path(dep_datasets$fit$path, 'CL_data.csv'), stringsAsFactors = F, check.names = F) 
CL_data %<>% full_join(full_join(true, fit, by = "CCLE_ID", suffix = c('_true', '_fit')), by = "CCLE_ID")

gene_slope_cor <- with(CL_data, cor(gene_slope_true, gene_slope_fit))


plot_colorByDensity(CL_data, 'gene_slope_true', 'gene_slope_fit', size = 1.5) +
  geom_abline() +
  xlab('True screen signal') + 
  ylab('Fit screen signal') +
  coord_equal() + 
  theme_Publication()
ggsave(file.path(fig_dir, 'gene_slope_sim_compare.png'), width = fig_width, height = fig_height, dpi = 350)

plot_colorByDensity(CL_data, 'CL_slope_true', 'CL_slope_fit', size = 1.5) +
  geom_abline() +
  xlab('True overall CL slope') + 
  ylab('Fit overall CL slope') +
  coord_equal() + 
  theme_Publication()
ggsave(file.path(fig_dir, 'CL_slope_sim_compare.png'), width = fig_width, height = fig_height, dpi = 350)


df <- data.frame(true = as.vector(dep_data$true$gene_scores), fit = as.vector(dep_data$fit$gene_scores))
plot_colorByDensity(df %>% sample_n(100000), 'true', 'fit') +
  geom_abline() +
  xlab('True gene dependency') + 
  ylab('Fit gene dependency') +
  coord_equal() + 
  theme_Publication()
ggsave(file.path(fig_dir, 'gene_dep_sim_compare.png'), width = fig_width, height = fig_height, dpi = 350)


df <- data.frame(true = as.vector(dep_data$true$gene_scores %>% scale(center = T, scale = F)), 
                 fit = as.vector(dep_data$fit$gene_scores %>% scale(center = T, scale = F)))
plot_colorByDensity(df %>% sample_n(100000), 'true', 'fit') +
  geom_abline() +
  coord_equal() + 
  xlab('True relative gene dependency') + 
  ylab('Fit relative gene dependency') +
  theme_Publication()
ggsave(file.path(fig_dir, 'gene_dep_ms_sim_compare.png'), width = fig_width, height = fig_height, dpi = 350)


# ggplot(CL_data, aes(CL_slope_true, CL_slope_fit)) + 
#   geom_point(alpha = 0.5) + 
#   geom_abline() +
#   scale_x_log10() +
#   scale_y_log10()

# ggplot(CL_data, aes(noise_vars_true, noise_vars_fit)) + 
#   geom_point(alpha = 0.5) + 
#   geom_abline() 

# get_R2 <- function(pred,true) {
#   SSE = sum((pred-true)^2, na.rm=T)
#   SST = sum((true - mean(true, na.rm=T))^2, na.rm=T)
#   return(1 - SSE/SST)
# }
# CL_stats <- ldply(rownames(dep_data$fit$gene_scores), function(CL) {
#   cc_p <- cor(dep_data$fit$gene_scores[CL,],
#                  dep_data$true$gene_scores[CL,],
#                  use = 'pairwise.complete.obs')
#   data.frame(cc_p = cc_p,
#              CCLE_ID = CL)
# })
# CL_data %<>% full_join(CL_stats, by = "CCLE_ID")
# ggplot(CL_data, aes(gene_slope_true, cc_p^2)) + 
#   geom_point() +
#   scale_x_log10()


gene_stats <- ldply(colnames(dep_data$fit$gene_scores), function(gene) {
  cc_p <- cor(dep_data$fit$gene_scores[,gene],
              dep_data$true$gene_scores[,gene],
              use = 'pairwise.complete.obs')
  gene_sd <- sd(dep_data$true$gene_scores[,gene], na.rm=T)
  data.frame(
    Gene = gene,
    cc_p = cc_p,
    gene_sd = gene_sd
  )
})
ggplot(gene_stats, aes(gene_sd, cc_p^2)) + 
  geom_point(alpha = 0.5) +
  geom_density_2d()

ov_R2 <- cor(as.vector(dep_data$fit$gene_scores), 
             as.vector(dep_data$true$gene_scores),
             use = 'pairwise.complete.obs')^2
# ov_wtd_R2 <- wtd.cor(as.vector(dep_data$fit$gene_scores),
#              as.vector(dep_data$true$gene_scores),
#              weight = 1/as.vector(dep_data$fit$gene_score_SDs)^2)[1, 'correlation']^2

ov_R2_ms <- cor(as.vector(dep_data$fit$gene_scores %>% scale(center = T, scale = F)), 
                as.vector(dep_data$true$gene_scores %>% scale(center = T, scale = F)),
                use = 'pairwise.complete.obs')^2
# ov_wtd_R2_ms <- wtd.cor(as.vector(dep_data$fit$gene_scores %>% scale(center = T, scale = F)), 
#                      as.vector(dep_data$true$gene_scores %>% scale(center = T, scale = F)),
#                      weight = 1/as.vector(dep_data$fit$gene_score_SDs)^2)[1, 'correlation']^2

all_results <- data.frame(ov_cor = sqrt(ov_R2), ov_cor_ms = sqrt(ov_R2_ms), posneg_SSMD = cur_SSMD, gene_slope_cor = gene_slope_cor,
                          n_CLs = nrow(dep_data$fit$gene_scores), n_genes = ncol(dep_data$fit$gene_scores))
write.csv(all_results, '~/CPDS/demeter2/kube_results/sim_run_n_CLs/full_sim/sim_results.csv', row.names = T)

# data.frame(true = orig_mod$CL_data$gene_slope, 
#            fit = new_mod$CL_data$gene_slope) %>% 
#   ggplot(aes(true, fit)) + geom_point() + geom_abline() + scale_x_log10() + 
#   scale_y_log10()
# 
# orig_GS_ms <- orig_mod$gene_scores %>% 
#   scale(center = T, scale = F)
# new_GS_ms <- new_mod$gene_scores %>% 
#   scale(center = T, scale = F)
# 
# # data.frame(true = orig_GS_ms['EFO21_OVARY',],
# #            fit = new_GS_ms['EFO21_OVARY',]) %>% 
# #   ggplot(aes(true, fit)) + geom_point(alpha = 0.5) + geom_abline()
# data.frame(true = orig_GS_ms['EFO21_OVARY',],
#            fit = new_GS_ms['EFO21_OVARY',]) %>% 
#   ggplot(aes(true, fit)) + geom_point(alpha = 0.5) + geom_abline()


#uncertainty estimates
fit_err <- dep_data$fit$gene_scores - dep_data$true$gene_scores
fit_t <- fit_err / dep_data$fit$gene_score_SDs
# ggplot() + 
#   geom_density(data = data.frame(fit_t = as.vector(fit_t)), aes(fit_t, fill = 'fit', color = 'fit'), lwd = 1, alpha = 0.1) +
#   geom_density(data = data.frame(norm = rnorm(10000)), aes(norm, fill = 'expected', color = 'expected'), lwd = 1, alpha = 0.1) +
#   xlim(-4,4) +
#   xlab('Gene dependency\nestimate error (t-stat)') +
#   guides(fill = FALSE, color = guide_legend(title = element_blank())) +
#   theme_Publication() +
#   scale_colour_Publication() + 
#   scale_fill_Publication()
ggplot() +
  geom_qq(data = data.frame(fit_t = as.vector(fit_t)) %>% 
            sample_n(500000), aes(sample = fit_t)) +
  geom_abline() +
  xlab('Standard normal quantiles') +
  ylab('Model standardized\nerror quantiles') +
  theme_Publication()
ggsave(file.path(fig_dir, 'gene_dep_t_stat_dist.png'), width = fig_width, height = fig_height, dpi = 350)




# df <- data.frame(fit_err = abs(as.vector(fit_err)), fit_unc = as.vector(dep_data$fit$gene_score_SDs))
# plot_colorByDensity(df %>% sample_n(50000), 'fit_err', 'fit_unc') +
# geom_smooth()
# ggplot() + 
#   geom_density(data = data.frame(fit_t = fit_t[50,]), aes(fit_t, fill = 'fit'), alpha = 0.25) + 
#   geom_density(data = data.frame(norm = rnorm(10000)), aes(norm, fill = 'norm'), alpha = 0.25) +
#   geom_vline(xintercept = 0)

# unc_biases <- adply(fit_t, 1, mad, na.rm=T)


