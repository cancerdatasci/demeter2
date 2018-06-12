library(plyr); 
library(magrittr); 
library(reshape2)
library(tibble)
library(dplyr)
library(taigr);
library(jsonlite)
library(ggplot2)
library(cowplot)
library(dependr)
library(readr)
library(stringr)

fig_dir <- '~/CPDS/demeter2/results/rev_figures'


# target_dir <- '~/CPDS/demeter2/kube_results/20170703-152054-5115'
# target_dir <- '~/CPDS/demeter2/kube_results/20170714-115027-86ad'
target_dir <- '~/CPDS/demeter2/kube_results/Ach_abs_scan'
model_dirs <- list.dirs(path = target_dir, full.names = TRUE, recursive = FALSE)

R2_res <- ldply(model_dirs, function(mod_dir) {
    print(mod_dir)
    mod_res <- jsonlite::read_json(file.path(mod_dir, 'other_info.json'))
    df <- data.frame(
        # test_R2_ms = mod_res$R2_vals$test_ms %>% unlist %>% tail(1),
        # test_R2 = mod_res$R2_vals$test %>% unlist %>% tail(1),
        # train_R2_ms = mod_res$R2_vals$train_ms %>% unlist %>% tail(1),
        train_R2 = mod_res$R2_vals$train %>% unlist %>% tail(1),
        train_SSMD = mod_res$SSMD$train %>% unlist %>% tail(1)
    )
    df %<>% cbind(data.frame(mod_res$reg_params))
    if (!is.null(mod_res$SSMD$test)) {
        df$test_SSMD <-  mod_res$SSMD$test %>% unlist %>% tail(1)
    }
    return(df)
})


cat('Optimal SSMD:\n')
best_params <- R2_res %>% arrange(desc(train_SSMD)) %>% head(1) %>% 
    dplyr::select(train_SSMD, gene_l2_lambda, seed_l2_lambda, hairpin_l2_lambda, hp_unpred_l2_lambda)
print(best_params)

# R2_res %>% group_by(gene_l2_lambda) %>% 
#   summarise(avg_SSMD = mean(train_SSMD)) %>% 
#   full_join(R2_res, by = 'gene_l2_lambda') %>% 
#   ggplot(aes(gene_l2_lambda, train_SSMD)) + 
#   geom_boxplot(aes(group = gene_l2_lambda)) + 
#   geom_point(aes(y = avg_SSMD), size = 3) + 
#   geom_line(aes(y = avg_SSMD)) + 
#   scale_x_log10()

g1 <- R2_res %>% 
  ggplot(aes(gene_l2_lambda, train_SSMD, group = gene_l2_lambda)) + 
  geom_boxplot() +
  scale_x_log10() +
  xlab('Avg gene-effect lambda') +
  ylab('Pos-neg control\nseparation (SSMD)')
g2 <- R2_res %>% 
  ggplot(aes(seed_l2_lambda, train_SSMD, group = seed_l2_lambda)) + 
  geom_boxplot() +
  scale_x_log10() +
  xlab('Avg seed-effect lambda') +
  ylab('Pos-neg control\nseparation (SSMD)')
g3 <- R2_res %>% 
  ggplot(aes(hairpin_l2_lambda, train_SSMD, group = hairpin_l2_lambda)) + 
  geom_boxplot() +
  scale_x_log10() +
  xlab('shRNA pDNA-offset lambda') +
  ylab('Pos-neg control\nseparation (SSMD)')
g4 <- R2_res %>% 
  ggplot(aes(hp_unpred_l2_lambda, train_SSMD, group = hp_unpred_l2_lambda)) + 
  geom_boxplot() +
  scale_x_log10() +
  xlab('shRNA overall offset lambda') +
  ylab('Pos-neg control\nseparation (SSMD)')

# g1 <- R2_res %>% group_by(gene_l2_lambda) %>% 
#     summarise(avg_SSMD = mean(train_SSMD)) %>% 
#     ggplot(aes(gene_l2_lambda, avg_SSMD)) + 
#     geom_line() + geom_point() +
#     scale_x_log10()
# g2 <- R2_res %>% group_by(seed_l2_lambda) %>% 
#     summarise(avg_SSMD = mean(train_SSMD)) %>% 
#     ggplot(aes(seed_l2_lambda, avg_SSMD)) + 
#     geom_line() + geom_point() +
#     scale_x_log10()
# g3 <- R2_res %>% group_by(hairpin_l2_lambda) %>% 
#     summarise(avg_SSMD = mean(train_SSMD)) %>% 
#     ggplot(aes(hairpin_l2_lambda, avg_SSMD)) + 
#     geom_line() + geom_point() +
#     scale_x_log10()
# g4 <- R2_res %>% group_by(hp_unpred_l2_lambda) %>% 
#     summarise(avg_SSMD = mean(train_SSMD)) %>% 
#     ggplot(aes(hp_unpred_l2_lambda, avg_SSMD)) + 
#     geom_line() + geom_point() +
#     scale_x_log10()
plot_grid(g1, g2, g3, g4, ncol =2)

ggsave(file.path(fig_dir, 'abs_hyperparam_boxplot.png'), width = 7, height = 7, dpi = 350)


R2_res[R2_res$gene_l2_lambda == best_params$gene_l2_lambda & R2_res$seed_l2_lambda == best_params$seed_l2_lambda, ] %>% 
    ggplot(aes(log10(hairpin_l2_lambda), log10(hp_unpred_l2_lambda), fill = train_SSMD)) +
    geom_raster()


R2_res[R2_res$hairpin_l2_lambda == best_params$hairpin_l2_lambda & R2_res$hp_unpred_l2_lambda == best_params$hp_unpred_l2_lambda, ] %>% 
    ggplot(aes(log10(gene_l2_lambda), log10(seed_l2_lambda), fill = train_SSMD)) +
    geom_raster()

# #compare results of R2 (on training sample)
# yr <- c(min(R2_res$train_R2), max(R2_res$train_R2))
# g1 <- R2_res %>% 
#   ggplot(aes(gene_l2_lambda, train_R2, group = gene_l2_lambda)) + 
#   geom_boxplot() +
#   scale_x_log10() + 
#   ylim(yr)
# g2 <- R2_res %>% 
#   ggplot(aes(seed_l2_lambda, train_R2, group = seed_l2_lambda)) + 
#   geom_boxplot() +
#   scale_x_log10()+ 
#   ylim(yr)
# g3 <- R2_res %>% 
#   ggplot(aes(hairpin_l2_lambda, train_R2, group = hairpin_l2_lambda)) + 
#   geom_boxplot() +
#   scale_x_log10()+ 
#   ylim(yr)
# g4 <- R2_res %>% 
#   ggplot(aes(hp_unpred_l2_lambda, train_R2, group = hp_unpred_l2_lambda)) + 
#   geom_boxplot() +
#   scale_x_log10()+ 
#   ylim(yr)
# g1 <- R2_res %>% group_by(gene_l2_lambda) %>% 
#   summarise(avg_R2 = mean(train_R2)) %>% 
#   ggplot(aes(gene_l2_lambda, avg_R2)) + 
#   geom_line() + geom_point() +
#   ylim(yr) +
#   scale_x_log10()
# g2 <- R2_res %>% group_by(seed_l2_lambda) %>% 
#   summarise(avg_R2 = mean(train_R2)) %>% 
#   ggplot(aes(seed_l2_lambda, avg_R2)) + 
#   geom_line() + geom_point() +
#   ylim(yr) +
#   scale_x_log10()
# g3 <- R2_res %>% group_by(hairpin_l2_lambda) %>% 
#   summarise(avg_R2 = mean(train_R2)) %>% 
#   ggplot(aes(hairpin_l2_lambda, avg_R2)) + 
#   geom_line() + geom_point() +
#   ylim(yr) +
#   scale_x_log10()
# g4 <- R2_res %>% group_by(hp_unpred_l2_lambda) %>% 
#   summarise(avg_R2 = mean(train_R2)) %>% 
#   ggplot(aes(hp_unpred_l2_lambda, avg_R2)) + 
#   geom_line() + geom_point() +
#   ylim(yr) +
#   scale_x_log10()

# plot_grid(g1, g2, g3, g4, ncol =2)

