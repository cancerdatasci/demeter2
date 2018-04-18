library(plyr); 
library(magrittr);  
library(reshape2)
library(stringr)
library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)
library(taigr)

dat <- load.all.from.taiga(list(
  Ach_GS = list( 
    data.name='demeter2-achilles-5386',
    data.version=10,
    data.file = 'gene_means_proc',
    transpose = T),
  Ach_CL_data = list( 
      data.name='demeter2-achilles-5386',
      data.version=10,
      data.file = 'CL_data_comb',
      transpose = F),
  Ach_SD = list( 
    data.name='demeter2-achilles-5386',
    data.version=10,
    data.file = 'gene_SDs_proc',
    transpose = T),
  Ach_D1_old  = list(data.name='demeter-2-20-2-z-scored-gene-solutions-expanded-gene-families-', 
                     data.version=6, 
                     data.file='data',
                     transpose = T)
))


#mean subtract per gene
dat$Ach_GS %<>% scale(center = T, scale = F)
dat$Ach_D1_old_r %<>% scale(center = T, scale = F)

#global z-score normalization
dat$Ach_SD  <- dat$Ach_SD / sd(dat$Ach_GS, na.rm=T)
dat$Ach_GS  <- dat$Ach_GS / sd(dat$Ach_GS, na.rm=T)
dat$Ach_D1_old  <- dat$Ach_D1_old / sd(dat$Ach_D1_old, na.rm=T)

gene_stats <- load.from.taiga(data.name='gene-distribution-stats-137d', data.version=1) %>% 
  mutate(Gene = str_match(Gene, '(.+) \\(')[,2])

dmap_genes <- read.csv('~/CPDS/data/Achilles/shRNA/RNAi_Cell_gene_dep_table.csv', stringsAsFactors = F, check.names = F)

dmap_six_sig <- dmap_genes %>% 
  left_join(gene_stats %>% filter(dset == 'RNAi_Ach') %>% 
              dplyr::select(`Gene dependency` = Gene, 
                            RNAi_Ach_normLRT = normLRT, 
                            RNAi_Ach_skew = skewness, 
                            RNAi_Ach_mean = mean)) %>% 
    left_join(data.frame(`Gene dependency` = str_match(colnames(dat$Ach_GS), '(.+) \\(')[,2],
                         D2_Ach_min_score = apply(dat$Ach_GS, 2, min, na.rm=T), 
                         D2_Ach_n_six_sig = colSums(dat$Ach_GS < -6, na.rm=T),
                         D2_Ach_n_five_sig = colSums(dat$Ach_GS < -5, na.rm=T),
                         check.names = FALSE)) %>% 
    left_join(data.frame(`Gene dependency` = str_match(colnames(dat$Ach_D1_old), '(.+) \\(')[,2],
                       D1o_Ach_min_score = apply(dat$Ach_D1_old, 2, min, na.rm=T), 
                       D1o_Ach_n_six_sig = colSums(dat$Ach_D1_old < -6, na.rm=T),
                       check.names = FALSE)) 
   
    

frac_D1o_six_sig <- with(dmap_six_sig, sum(D1o_Ach_min_score < -6, na.rm=T) / sum(!is.na(D1o_Ach_min_score)))
equiv_D2_thresh <- with(dmap_six_sig %>% filter(!is.na(D2_Ach_min_score)), quantile(D2_Ach_min_score, frac_D1o_six_sig))

with(dmap_six_sig %>% filter(is.six.sigma), mean(D2_Ach_min_score < equiv_D2_thresh, na.rm=T))
with(dmap_six_sig %>% filter(is.six.sigma, D1o_Ach_n_six_sig >= 2), mean(D2_Ach_min_score < equiv_D2_thresh, na.rm=T))
with(dmap_six_sig %>% filter(is.six.sigma, D1o_Ach_n_six_sig >= 3), mean(D2_Ach_min_score < equiv_D2_thresh, na.rm=T))

ggplot(dmap_six_sig, aes(D1o_Ach_min_score, D2_Ach_min_score, color = is.six.sigma)) +
    geom_point(alpha = 0.5) + 
    # geom_abline() + 
    geom_hline(yintercept = equiv_D2_thresh, linetype = 'dashed') + 
    geom_vline(xintercept = -6, linetype = 'dashed') +
    xlab('D1 (original) strongest dependency (z)') +
    ylab('D2 strongest dependency (z)')
ggsave('~/CPDS/demeter2/results/figures/strongest_dep_scatter.png', width = 6, height = 4)


CL_qual <- read.csv('~/CPDS/demeter2/kube_results/Ach_final/1/CL_data_comb.csv', stringsAsFactors = F, check.names = F)


CL_w_min <- rownames(dat$Ach_D1_old)[apply(dat$Ach_D1_old, 2, function(col) {
    which.min(col)
})] 
D1_min <- data.frame(Gene = colnames(dat$Ach_D1_old),
                 CCLE_ID = CL_w_min) %>% 
    mutate(Gene = str_match(Gene, '(.+) \\(')[,2]) %>% 
    left_join(CL_qual) 

CL_w_min <- rownames(dat$Ach_GS)[apply(dat$Ach_GS, 2, function(col) {
    which.min(col)
})] 
D2_min <- data.frame(Gene = colnames(dat$Ach_GS),
                 CCLE_ID = CL_w_min) %>% 
    mutate(Gene = str_match(Gene, '(.+) \\(')[,2]) %>% 
    left_join(CL_qual) 

dmap_six_sig %<>% left_join(D1_min %>% 
                                dplyr::select(gene_slope_D1 = gene_slope, CCLE_ID_D1 = CCLE_ID, Gene), 
                            by = c('Gene dependency' = 'Gene')) %>% 
    left_join(D2_min %>% dplyr::select(gene_slope_D2 = gene_slope, CCLE_ID_D2 = CCLE_ID, Gene), by = c('Gene dependency' = 'Gene'))

ggplot(dmap_six_sig, aes(D1o_Ach_min_score, D2_Ach_min_score, color = gene_slope_D1)) +
    geom_point(alpha = 0.75) + 
    # geom_abline() + 
    geom_hline(yintercept = equiv_D2_thresh, linetype = 'dashed') + 
    geom_vline(xintercept = -6, linetype = 'dashed') +
    scale_color_gradient2(midpoint = 1, limits = c(0.5, 2), oob = scales::squish)

ggplot(dmap_six_sig, aes(D1o_Ach_min_score, D2_Ach_min_score, color = gene_slope_D2)) +
    geom_point(alpha = 0.75) + 
    # geom_abline() + 
    geom_hline(yintercept = equiv_D2_thresh, linetype = 'dashed') + 
    geom_vline(xintercept = -6, linetype = 'dashed') +
    scale_color_gradient2(midpoint = 1, limits = c(0.5, 2), oob = scales::squish)

with(dmap_six_sig %>% filter(is.six.sigma), median(gene_slope_D1, na.rm=T))
with(dmap_six_sig %>% filter(is.six.sigma), median(gene_slope_D2, na.rm=T))

dmap_six_sig %>% 
    # filter(is.six.sigma) %>% 
    ggplot() + 
    geom_density(data = filter(dmap_six_sig, is.six.sigma), aes(gene_slope_D1, color = 'D1')) +
    geom_density(data = filter(dmap_six_sig, is.six.sigma), aes(gene_slope_D2, color = 'D2')) +
    geom_density(data = CL_qual, aes(gene_slope, color = 'overall')) +
    xlab('Most dependent CL\nscreen signal') +
    scale_x_log10(breaks = c(0.5, 1, 2)) +
    guides(color = guide_legend(title = element_blank())) +
    geom_vline(xintercept = 1, linetype = 'dashed')
ggsave('~/CPDS/demeter2/results/figures/outlier_CL_screen_signal.png', width = 6, height = 4)


dmap_six_sig %>% filter(is.six.sigma) %>% 
    group_by(CCLE_ID_D1) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n)) %>% 
    head(10) %>% 
    .[['n']] %>% 
    sum()

dmap_six_sig %>% filter(is.six.sigma) %>% 
    group_by(CCLE_ID_D2) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n)) %>% 
    head(10) %>% 
    .[['n']] %>% 
    sum()


with(dmap_six_sig %>% filter(is.six.sigma), wilcox.test(gene_slope_D1-1, na.rm=T))
with(dmap_six_sig %>% filter(is.six.sigma), wilcox.test(gene_slope_D2-1, na.rm=T))
with(dmap_six_sig %>% filter(!is.six.sigma, D2_Ach_min_score < equiv_D2_thresh), median(gene_slope_D1, na.rm=T))
with(dmap_six_sig %>% filter(is.six.sigma, D2_Ach_min_score > equiv_D2_thresh), median(gene_slope_D1, na.rm=T))

with(dmap_six_sig %>% filter(is.six.sigma), wilcox.test(gene_slope_D1[D2_Ach_min_score > equiv_D2_thresh], gene_slope_D1[D2_Ach_min_score < equiv_D2_thresh], na.rm=T))


with(dmap_six_sig %>% filter(is.six.sigma), median(gene_slope_D2, na.rm=T))
with(dmap_six_sig %>% filter(!is.six.sigma, D2_Ach_min_score < equiv_D2_thresh), median(gene_slope_D2, na.rm=T))
with(dmap_six_sig %>% filter(is.six.sigma, D2_Ach_min_score > equiv_D2_thresh), median(gene_slope_D2, na.rm=T))
