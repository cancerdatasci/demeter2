library(plyr); 
library(magrittr); 
library(reshape2)
library(tibble)
library(dplyr)
library(taigr);
library(jsonlite)
library(ggplot2)
library(dependr)
library(readr)
library(stringr)

target_dir <- '~/CPDS/demeter2/kube_results/20170715-171303-259b'
model_dirs <- list.dirs(path = target_dir, full.names = TRUE, recursive = FALSE)

R2_res <- ldply(model_dirs, function(mod_dir) {
    mod_res <- jsonlite::read_json(file.path(mod_dir, 'other_info.json'))
    df <- data.frame(
        test_R2_ms = mod_res$R2_vals$test_ms %>% unlist %>% tail(1),
        test_R2 = mod_res$R2_vals$test %>% unlist %>% tail(1),
        train_R2_ms = mod_res$R2_vals$train_ms %>% unlist %>% tail(1),
        train_R2 = mod_res$R2_vals$train %>% unlist %>% tail(1),
        train_SSMD = mod_res$SSMD$train %>% unlist %>% tail(1)
    )
    df %<>% cbind(data.frame(mod_res$reg_params))
    if (!is.null(mod_res$SSMD$test)) {
       df$test_SSMD <-  mod_res$SSMD$test %>% unlist %>% tail(1)
    }
    return(df)
})

cat('Optimal test R2:\n')
R2_res %>% arrange(desc(test_R2)) %>% head(1) %>% 
    dplyr::select(test_R2, rel_gene_l2_lambda, rel_seed_l2_lambda, gene_l2_lambda, rel_gene_l2_lambda)
ggplot(R2_res, aes(log10(rel_gene_l2_lambda), log10(rel_seed_l2_lambda), fill = test_R2)) + 
    geom_raster()

cat('Optimal test R2 mean-sub:\n')
R2_res %>% arrange(desc(test_R2_ms)) %>% head(1) %>% 
    dplyr::select(test_R2_ms, rel_gene_l2_lambda, rel_seed_l2_lambda, gene_l2_lambda, rel_gene_l2_lambda)
ggplot(R2_res, aes(log10(rel_gene_l2_lambda), log10(rel_seed_l2_lambda), fill = test_R2_ms)) + 
    geom_raster()


cat('Optimal SSMD:\n')
R2_res %>% arrange(desc(train_SSMD)) %>% head(1) %>% 
    dplyr::select(train_SSMD, rel_gene_l2_lambda, rel_seed_l2_lambda, gene_l2_lambda, rel_gene_l2_lambda)
ggplot(R2_res, aes(log10(rel_gene_l2_lambda), log10(rel_seed_l2_lambda), fill = train_SSMD)) + 
    geom_raster()

