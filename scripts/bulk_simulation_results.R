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

fig_dir <- '~/CPDS/demeter2/results/figures'
fig_width <- 3.5
fig_height <- 3.0

sim_folder <- '~/CPDS/demeter2/kube_results/sim_run_n_CLs'
res_dirs <- list.files(sim_folder, include.dirs = TRUE)
all_res <- ldply(res_dirs, function(cur_res_dir) {
  print(cur_res_dir)
  read.csv(file.path(sim_folder, cur_res_dir, 'sim_results.csv'), stringsAsFactors = F, check.names = F, row.names = 1)
})

sim_folder <- '~/CPDS/demeter2/kube_results/sim_run_n_CLs2'
res_dirs <- list.files(sim_folder, include.dirs = TRUE)
all_res %<>% rbind(
  ldply(res_dirs, function(cur_res_dir) {
  print(cur_res_dir)
  read.csv(file.path(sim_folder, cur_res_dir, 'sim_results.csv'), stringsAsFactors = F, check.names = F, row.names = 1)
})
)

GA_res <- read.csv('~/CPDS/demeter2/results/GA_sim_results_n_CLs.csv', stringsAsFactors = F, check.names = F)
all_res %<>%
  mutate(model = 'D2') %>% 
  rbind.fill(GA_res %>% mutate(model = 'GA'))


all_res_g <- all_res %>% 
  mutate(gene_slope_R2 = gene_slope_cor^2,
         ov_R2 = ov_cor^2,
         ov_R2_ms = ov_cor_ms^2) %>% 
  group_by(n_CLs, n_genes, model) %>% 
  summarise_all(funs(mean, sd))

# all_res_g %>% 
#   ggplot(aes(n_CLs, CL_slope_cor_mean)) +
#   geom_point() + 
#   geom_line() +
#   geom_errorbar(aes(ymax = CL_slope_cor_mean + CL_slope_cor_sd, ymin = CL_slope_cor_mean - CL_slope_cor_sd))

## N_CLs vs screen signal R2
all_res_g %>% 
  filter(model == 'D2') %>% 
  filter(n_CLs > 2) %>% 
  ggplot(aes(n_CLs, gene_slope_R2_mean)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(aes(ymax = gene_slope_R2_mean + gene_slope_R2_sd, ymin = gene_slope_R2_mean - gene_slope_R2_sd)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_x_log10() + 
  theme_Publication() +
  xlab('N Cell Lines') + 
  ylab('Screen signal fit R2')
ggsave(file.path(fig_dir, 'gene_slope_R2_vs_n_CLs.png'), width = fig_width, height = fig_height, dpi = 350)


## N_CLs vs posneg SSMD
all_res_g %>% 
  ggplot(aes(n_CLs, posneg_SSMD_mean, color = model, group = model)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(aes(ymax = posneg_SSMD_mean + posneg_SSMD_sd, ymin = posneg_SSMD_mean - posneg_SSMD_sd)) +
  scale_x_log10() + 
  theme_Publication() +
  xlab('N Cell Lines') + 
  ylab('Positive-negative control\nseparation (SSMD)') +
  scale_color_manual(values = mod_type_pal)
ggsave(file.path(fig_dir, 'SSMD_vs_n_CLs.png'), width = fig_width, height = fig_height, dpi = 350)


## N_Cls vs overall gene dep R2
all_res_g %>% 
  ggplot(aes(n_CLs, ov_R2_mean, color = model, group = model)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(aes(ymax = ov_R2_mean + ov_R2_sd, ymin = ov_R2_mean - ov_R2_sd)) +
  scale_x_log10() + 
  theme_Publication() +
  xlab('N Cell Lines') + 
  ylab('Gene dependency fit (R2)') +
  scale_color_manual(values = mod_type_pal)
ggsave(file.path(fig_dir, 'gene_dep_ov_R2_vs_n_CLs.png'), width = fig_width, height = fig_height, dpi = 350)


all_res_g %>% 
  ggplot(aes(n_CLs, ov_R2_ms_mean, color = model, group = model)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(aes(ymax = ov_R2_ms_mean + ov_R2_ms_sd, ymin = ov_R2_ms_mean - ov_R2_ms_sd)) +
  scale_x_log10() + 
  theme_Publication() +
  xlab('N Cell Lines') + 
  ylab('Relative gene\ndependency fit (R2)') +
  scale_color_manual(values = mod_type_pal)
ggsave(file.path(fig_dir, 'gene_dep_rel_R2_vs_n_CLs.png'), width = fig_width, height = fig_height, dpi = 350)



####### FOR N_GENES
sim_folder <- '~/CPDS/demeter2/kube_results/sim_run_n_genes'
res_dirs <- list.files(sim_folder, include.dirs = TRUE)

all_res <- ldply(res_dirs, function(cur_res_dir) {
  print(cur_res_dir)
  read.csv(file.path(sim_folder, cur_res_dir, 'sim_results.csv'), stringsAsFactors = F, check.names = F, row.names = 1)
})

GA_res <- read.csv('~/CPDS/demeter2/results/GA_sim_results_n_genes.csv', stringsAsFactors = F, check.names = F)
all_res %<>%
  mutate(model = 'D2') %>% 
  rbind.fill(GA_res %>% mutate(model = 'GA'))

all_res %<>% mutate(n_genes = pmin(n_genes, 16000))
all_res_g <- all_res %>% 
  mutate(gene_slope_R2 = gene_slope_cor^2,
         ov_R2 = ov_cor^2,
         ov_R2_ms = ov_cor_ms^2) %>% 
  group_by(n_CLs, n_genes, model) %>% 
  summarise_all(funs(mean, sd))

# all_res_g %>% 
#   ggplot(aes(n_genes, CL_slope_cor_mean)) +
#   geom_point() + 
#   geom_line() +
#   geom_errorbar(aes(ymax = CL_slope_cor_mean + CL_slope_cor_sd, ymin = CL_slope_cor_mean - CL_slope_cor_sd))


all_res_g %>% 
  filter(model == 'D2') %>% 
  ggplot(aes(n_genes, gene_slope_R2_mean)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(aes(ymax = gene_slope_R2_mean + gene_slope_R2_sd, ymin = gene_slope_R2_mean - gene_slope_R2_sd)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_x_log10() +
  theme_Publication() +
  xlab('N genes') + 
  ylab('Screen signal fit R2')
ggsave(file.path(fig_dir, 'gene_slope_R2_vs_n_genes.png'), width = fig_width, height = fig_height, dpi = 350)


all_res_g %>% 
  ggplot(aes(n_genes, posneg_SSMD_mean, color = model, group = model)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(aes(ymax = posneg_SSMD_mean + posneg_SSMD_sd, ymin = posneg_SSMD_mean - posneg_SSMD_sd)) +
  scale_x_log10() +
  theme_Publication() +
  xlab('N genes') + 
  ylab('Positive-negative control\nseparation (SSMD)') +
  scale_color_manual(values = mod_type_pal)
ggsave(file.path(fig_dir, 'SSMD_vs_n_genes.png'), width = fig_width, height = fig_height, dpi = 350)


all_res_g %>% 
  ggplot(aes(n_genes, ov_R2_mean, color = model, group = model)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(aes(ymax = ov_R2_mean + ov_R2_sd, ymin = ov_R2_mean - ov_R2_sd)) +
  scale_x_log10() +
  theme_Publication() +
  xlab('N genes') + 
  ylab('Gene dependency fit R2') +
  scale_color_manual(values = mod_type_pal)
ggsave(file.path(fig_dir, 'gene_dep_ov_R2_vs_n_genes.png'), width = fig_width, height = fig_height, dpi = 350)

all_res_g %>% 
  ggplot(aes(n_genes, ov_R2_ms_mean, color = model, group = model)) +
  geom_point() + 
  geom_line() +
  geom_errorbar(aes(ymax = ov_R2_ms_mean + ov_R2_ms_sd, ymin = ov_R2_ms_mean - ov_R2_ms_sd)) +
  scale_x_log10() +
  theme_Publication() +
  xlab('N genes') + 
  ylab('Relative gene\ndependency fit (R2)') +
  scale_color_manual(values = mod_type_pal)
ggsave(file.path(fig_dir, 'gene_dep_rel_R2_vs_n_genes.png'), width = fig_width, height = fig_height, dpi = 350)

