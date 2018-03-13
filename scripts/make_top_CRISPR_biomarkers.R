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
library(cdsr)
library(doMC)

CERES_CRISPR = load.from.taiga(
    data.name = 'avana-public-tentative-18q1-92d9',
    data.version = 5,
    data.file = 'gene_effect',
    transpose = F)

feature_datasets <- list(
    GE = list(data.name = 'ccle-rnaseq-logrpkm-protein-coding-2700',
              data.version = 2, transpose = T),
    CN = list(data.name = 'depmap-wes-cn-data-97cc', 
              data.version = 5, 
              data.file = 'public_18q1_gene_cn'),
    MUT = list(data.name = 'ccle-mutation-data-full-', data.version = 2)
)

min_mutants <- 5
min_var <- 0.01
top_N <- 1000

features <- load.all.from.taiga(feature_datasets)

#parse hot-spot misense and damaging mutation datasets
features$MUT %<>% mutate(Gene = paste0(Hugo_Symbol, ' (', Entrez_Gene_Id, ')'),
                         cell_line = Tumor_Sample_Barcode)

features$MUT_HOT <- features$MUT %>% 
    filter(
        Variant_Classification == 'Missense_Mutation',
        (isTCGAhotspot) | (isCOSMIChotspot)) %>% 
    dplyr::select(cell_line, Gene) %>% 
    mutate(mutation = 1) %>% 
    unique() %>% 
    reshape2::acast(cell_line ~ Gene, value.var = 'mutation', fill = 0)
features$MUT_DAM <- features$MUT %>% 
    filter(isDeleterious) %>% 
    dplyr::select(cell_line, Gene) %>% 
    mutate(mutation = 1) %>% 
    unique() %>% 
    reshape2::acast(cell_line ~ Gene, value.var = 'mutation', fill = 0)

# rownames(features$CN) <- str_match(rownames(features$CN), '[:alnum:]+_(.+)')[,2]

features <- features[setdiff(names(features), 'MUT')]

#keep only features that have at least minimum variance over common cell lines
features <- lapply(features, function(feature) {
    new_feature <- feature[intersect(rownames(feature), rownames(CERES_CRISPR)), ]
    new_feature <- new_feature[, apply(new_feature, 2, var, na.rm=T) > min_var]
    return(new_feature)
})

#for mutation features use only mutations with min number mutant lines
features$MUT_DAM <- features$MUT_DAM[, colSums(features$MUT_DAM) >= min_mutants]
features$MUT_HOT <- features$MUT_HOT[, colSums(features$MUT_HOT) >= min_mutants]

#get top top_N feature-dep correlations, one per dep
all_top_cors <- ldply(names(features), function(cur_feat_type) {
    cur_features <- features[[cur_feat_type]]
    cur_CLs <- intersect(rownames(cur_features), rownames(CERES_CRISPR))
    cors <- cor(CERES_CRISPR[cur_CLs,], cur_features[cur_CLs,])
    melt(cors) %>% 
        set_colnames(c('Dep_gene', 'Feature_gene', 'cor')) %>% 
        group_by(Dep_gene) %>% 
        top_n(n = 1, wt = abs(cor)) %>% 
        mutate(feat_type = cur_feat_type) %>% 
        arrange(desc(abs(cor))) %>% 
        head(top_N)
})

write.csv(all_top_cors, '~/CPDS/demeter2/results/top_CRISPR_biomarkers.csv', row.names = F)

