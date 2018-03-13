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
        # train_R2 = mod_res$R2_vals$train %>% unlist %>% tail(1),
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

g1 <- R2_res %>% group_by(gene_l2_lambda) %>% 
    summarise(avg_SSMD = mean(train_SSMD)) %>% 
    ggplot(aes(gene_l2_lambda, avg_SSMD)) + 
    geom_line() + geom_point() +
    scale_x_log10()
g2 <- R2_res %>% group_by(seed_l2_lambda) %>% 
    summarise(avg_SSMD = mean(train_SSMD)) %>% 
    ggplot(aes(seed_l2_lambda, avg_SSMD)) + 
    geom_line() + geom_point() +
    scale_x_log10()
g3 <- R2_res %>% group_by(hairpin_l2_lambda) %>% 
    summarise(avg_SSMD = mean(train_SSMD)) %>% 
    ggplot(aes(hairpin_l2_lambda, avg_SSMD)) + 
    geom_line() + geom_point() +
    scale_x_log10()
g4 <- R2_res %>% group_by(hp_unpred_l2_lambda) %>% 
    summarise(avg_SSMD = mean(train_SSMD)) %>% 
    ggplot(aes(hp_unpred_l2_lambda, avg_SSMD)) + 
    geom_line() + geom_point() +
    scale_x_log10()

plot_grid(g1, g2, g3, g4, ncol =2)


R2_res[R2_res$gene_l2_lambda == best_params$gene_l2_lambda & R2_res$seed_l2_lambda == best_params$seed_l2_lambda, ] %>% 
    ggplot(aes(log10(hairpin_l2_lambda), log10(hp_unpred_l2_lambda), fill = train_SSMD)) +
    geom_raster()


R2_res[R2_res$hairpin_l2_lambda == best_params$hairpin_l2_lambda & R2_res$hp_unpred_l2_lambda == best_params$hp_unpred_l2_lambda, ] %>% 
    ggplot(aes(log10(gene_l2_lambda), log10(seed_l2_lambda), fill = train_SSMD)) +
    geom_raster()

