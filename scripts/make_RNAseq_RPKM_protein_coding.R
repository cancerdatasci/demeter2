library(plyr);
library(magrittr);
library(reshape2)
library(dplyr)
library(taigr);
library(dependr)
library(readr)
library(stringr)
library(cdsr)
library(biomaRt)

GE <- read.gct('~/CPDS/data/CCLE/CCLE_DepMap_18Q1_RNAseq_RPKM_20180214.gct')
ensembl_ids <- rownames(GE)

ensembl_ids_nv <- str_match(ensembl_ids, '(.+)\\.')[,2]

att <- c('entrezgene', 'ensembl_peptide_id', 'ensembl_gene_id', 'hgnc_symbol')
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_features <- getBM(att, filters = 'ensembl_gene_id', values = ensembl_ids_nv, mart = ensembl)
used_gene_features <- gene_features %>%
  filter(ensembl_peptide_id != '') %>%
  distinct(ensembl_gene_id, .keep_all=T)
used_gene_features %<>% left_join(data.frame(ensembl_gene_id = ensembl_ids_nv, full_ensembl_id = ensembl_ids, stringsAsFactors = F), by = 'ensembl_gene_id')
used_gene_features %<>% mutate(gene_name = paste0(hgnc_symbol, ' (', entrezgene, ')'))
used_gene_features %<>% filter(full_ensembl_id %in% rownames(GE))
GE_used <- GE[used_gene_features$full_ensembl_id,] %>%
  set_rownames(used_gene_features$gene_name)

GE_used <- log2(GE_used + 0.001)
entrez_ID <- str_match(rownames(GE_used), ' \\((.+)\\)')[, 2]
GE_used <- GE_used[!duplicated(entrez_ID),]
write.csv(GE_used, '~/CPDS/data/CCLE/CCLE_DepMap_18Q1_RNAseq_log2RPKM_protein_coding.csv', row.names=TRUE)
