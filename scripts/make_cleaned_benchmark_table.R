library(biomaRt)
library(plyr); library(magrittr);  
library(reshape2)
library(tibble)
library(dplyr)
library(taigr);
library(ggplot2)
library(ggrepel)
library(dependr)
library(readr)
library(stringr)

benchmark_rel <- read.csv('~/CPDS/data/general/benchmark_list_cleaned.csv', 
                          stringsAsFactors = F, check.names = F) %>% 
    filter(Feature_matrix %in% c('CN', 'EXP', 'MUT')) %>% 
    dplyr::select(Feature, Dependency, Feature_type = Feature_matrix) %>% 
    distinct(Feature, Dependency, Feature_type, .keep_all=T)

my_gene_names <- union(unique(benchmark_rel$Feature), benchmark_rel$Dependency)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
mapTab <- getBM(attributes = c("hgnc_symbol", "entrezgene"), 
                filters = "hgnc_symbol", 
                values = my_gene_names, 
                mart = ensembl, 
                uniqueRows=FALSE) %>% 
    distinct(hgnc_symbol, entrezgene)
stopifnot(sum(duplicated(mapTab$hgnc_symbol)) == 0)
symbol_to_ID <- as.character(mapTab$entrezgene) %>% set_names(mapTab$hgnc_symbol)

benchmark_rel %<>% mutate(Feature_symbol = Feature,
                          Dependency_symbol = Dependency,
                          Feature = symbol_to_ID[Feature],
                          Dependency = symbol_to_ID[Dependency])

write.csv(benchmark_rel, '~/CPDS/data/general/benchmark_list_cleaned_ID.csv', row.names = F) 
                          