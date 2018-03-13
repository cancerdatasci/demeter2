library(readr); library(stringr)
library(plyr); library(magrittr); 
library(reshape2)
library(taigr)
library(tibble)
library(dplyr)

rel_gene_L2 <- 10 ^ seq(from = log10(0.1), to = log10(10), length.out = 7)
rel_seed_L2 <- 10 ^ seq(from = log10(0.1), to = log10(10), length.out = 7)
# rel_gene_L2 <- c(0.01, 1)
# rel_seed_L2 <- c(0.01, 1)
gene_L2 <- 1
seed_L2 <- 1
hp_L2 <- 10
hp_unpred_L2 <- 10

df <- expand.grid(gene_rel_L2 = rel_gene_L2, 
                  seed_rel_L2 = rel_seed_L2,
                  gene_L2 = gene_L2,
                  seed_L2 = seed_L2,
                  hp_offset_L2 = hp_L2,
                  hp_unpred_L2 = hp_unpred_L2) %>% 
    rownames_to_column(var = 'Model_num') %>% 
    mutate(Model_num = paste0('Mod_', Model_num))

write.csv(df, '~/CPDS/demeter2/params.csv', row.names=F)
write.csv(df %>% head(1), '~/CPDS/demeter2/params_test.csv', row.names=F)


#scan of other params
rel_gene_L2 <- c(1.0)
rel_seed_L2 <- c(2.0)
abs_L2 <- 10 ^ seq(from = 0, to = 2, length.out = 3)
hp_L2 <- 10 ^ seq(from = 0, to = 2, length.out = 3)

df <- expand.grid(gene_rel_L2 = rel_gene_L2, 
                  seed_rel_L2 = rel_seed_L2,
                  gene_L2 = abs_L2,
                  seed_L2 = abs_L2,
                  hp_offset_L2 = hp_L2,
                  hp_unpred_L2 = hp_L2) %>% 
    rownames_to_column(var = 'Model_num') %>% 
    mutate(Model_num = paste0('Mod_', Model_num))

write.csv(df, '~/CPDS/demeter2/abs_params.csv', row.names=F)
write.csv(df %>% head(1), '~/CPDS/demeter2/abs_params_test.csv', row.names=F)

