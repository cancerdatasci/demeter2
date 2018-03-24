library(plyr)
library(dplyr)
library(magrittr)
library(readr)
library(taigr)
library(readr)
library(cdsr)

inWeb <- read_tsv('~/CPDS/data/InBio_Map_core_2016_09_12/core.psimitab', 
                  col_names = c('Interactor_A', 'Interactor_B', 
                                'Alt_name_A', 'Alt_name_B',
                                'Alias_A', 'Alias_B',
                                'Detection_methods',
                                'First_Author', 
                                'publication',
                                'NCBI_Taxon_A', 'NCBI_Taxon_B',
                                'Interaction_types', 
                                'Source_database', 
                                'Interaction_ID',
                                'Confidence_score',
                                'Complex_expansion')) %>% 
  mutate(conf = as.numeric(str_match(Confidence_score, '(.+)\\|')[,2]),
         lit_conf = as.numeric(str_match(Confidence_score, '\\|(.+)')[,2]))

inWeb %<>% mutate(gene_A = str_match(Alt_name_A, 'ensembl:(ENSG[0-9]+)\\|')[,2],
                  gene_B = str_match(Alt_name_B, 'ensembl:(ENSG[0-9]+)\\|')[,2])
  
all_ensembl_genes <- union(inWeb$gene_A, inWeb$gene_B)
ensembl_to_entrez <- cdsr::map_genes(all_ensembl_genes, 'ensembl', 'entrez')
ensembl_to_entrez <- ensembl_to_entrez$mapping$entrez %>% set_names(ensembl_to_entrez$mapping$ensembl)

inWeb %<>% mutate(Entrez_A = as.character(ensembl_to_entrez[gene_A]),
                  Entrez_B = as.character(ensembl_to_entrez[gene_B])) %>% 
  filter(!is.na(Entrez_A), !is.na(Entrez_B))
inWeb %>% 
  dplyr::select(Entrez_A, Entrez_B, conf) %>% 
  write.csv('~/CPDS/data/InBio_Map_core_2016_09_12/Inweb_pairs_Entrez_mapped.csv', row.names = F)

