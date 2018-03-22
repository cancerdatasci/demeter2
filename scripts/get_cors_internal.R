library(readr)
library(weights)
library(magrittr)
input_cache_file <- '~/CPDS/demeter2/results/temp_input_cache.rds'
output_res_file <- '~/CPDS/demeter2/results/temp_output_cache.rds'

inputs <- read_rds(input_cache_file)
res <- ldply(inputs$genes, function(gene) {
    if (inputs$cur_use_bayes) {
        cur_cor <- wtd.cors(inputs$ref,
                            inputs$dep[, gene],
                            1/inputs$SD[, gene]^2)
    } else {
        cur_cor <- cor(inputs$ref,
                       inputs$dep[, gene],
                       use = 'pairwise.complete.obs')
    }
    cur_cor <- as.vector(cur_cor)
    cur_cor[is.infinite(cur_cor)] <- NA
    best_ind <- which.max(abs(cur_cor))
    cor_df <- data.frame(cor = cur_cor[best_ind],
                         feat_gene = colnames(inputs$ref)[best_ind],
                         targ_gene = gene)
    return(cor_df)
}, .parallel=TRUE) %>% mutate(dset = inputs$cur_dset)

write_rds(res, output_res_file)
