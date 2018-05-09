library(plyr)
library(dplyr)
library(magrittr)
library(readr)
library(taigr)

CCLE_name_correction = read.csv('~/CPDS/demeter2/data/CCLE_name_corrections.csv', check.names = F, stringsAsFactors = F)
CCLE_name_map <- CCLE_name_correction$new_name %>% set_names(CCLE_name_correction$old_name)

lineage_info <- load.from.taiga(data.name='lineage-f95f', data.version=5) %>% 
  dplyr::rename(CCLE_ID = CCLE_name,
                disease = primary_tissue,
                disease_subtype = secondary_tissue,
                disease_sub_subtype = tertiary_tissue)
Marc_CL_data <- load.from.taiga(data.name='demeter2-marcotte-a703', data.version=5, data.file='CL_data_comb')
DRIVE_CL_data <- load.from.taiga(data.name='demeter2-drive-0591', data.version=7, data.file='CL_data_comb')
Ach_CL_data <- load.from.taiga(data.name='demeter2-achilles-5386', data.version=9, data.file='CL_data_comb')
Comb_CL_data <- load.from.taiga(data.name='demeter2-combined-dc9c', data.version=6, data.file='CL_data_comb')

sample_info <- data.frame(CCLE_ID = rownames(Comb_CL_data)) %>% 
  left_join(lineage_info, by = 'CCLE_ID') %>% 
  mutate(in_DRIVE = CCLE_ID %in% rownames(DRIVE_CL_data),
         in_Achilles = CCLE_ID %in% rownames(Ach_CL_data),
         in_Marcotte = CCLE_ID %in% rownames(Marc_CL_data))

load('~/CPDS/data/CCLE/Annotations.RData')
Novartis_sample_info <- read.csv('~/CPDS/data/Novartis/Novartis_Sample_info.csv', 
                                 check.names = FALSE, 
                                 stringsAsFactors = FALSE) %>% 
  dplyr::mutate(Novartis_name = CELLLINE,
                CELLLINE = toupper(CELLLINE),
                CLEANNAME_PRIMARYSITE = toupper(CLEANNAME_PRIMARYSITE)) %>% 
  dplyr::mutate(CCLE_ID = CleanCellLineName(CELLLINE),
                CCLE_ID = ifelse(is.na(CCLE_ID), CLEANNAME_PRIMARYSITE, CCLE_ID)) %>% 
  dplyr::select(CCLE_ID, Novartis_name,
                Novartis_Primary_site = PRIMARY_SITE, Novartis_Pathologist_Annotation = PATHOLOGIST_ANNOTATION) %>% 
  mutate(CCLE_ID = revalue(CCLE_ID, c(`GISTT1_GASTROINTESTINAL_TRACT_(SITE_INDETERMINATE)` = 'GISTT1_GASTROINTESTINAL_TRACT')))

sample_info %<>% left_join(Novartis_sample_info, by = 'CCLE_ID')

Marcotte_annot <- read_tsv('~/CPDS/data/Marcotte/cell_line_subtypes.txt') %>% 
  as.data.frame() %>% 
  mutate(Marcotte_name = cell_line,
         cell_line = toupper(cell_line),
         CCLE_ID = CleanCellLineName(cell_line),
         CCLE_ID = ifelse(is.na(CCLE_ID), paste0(cell_line, '_BREAST'), CCLE_ID)) %>% 
  dplyr::select(CCLE_ID, 
                Marcotte_name,
                Marcotte_subtype_three_receptor = subtype_three_receptor,
                Marcotte_subtype_neve = subtype_neve,
                Marcotte_subtype_intrinsic = subtype_intrinsic)

sample_info %<>% left_join(Marcotte_annot, by = 'CCLE_ID')
sample_info %<>% mutate(disease = ifelse(in_Marcotte & is.na(disease), 'breast', disease))

sample_info %<>% mutate(CCLE_ID = plyr::revalue(CCLE_ID, CCLE_name_map))

write.csv(sample_info, '~/CPDS/data/D2_figshare/sample_info.csv', row.names = FALSE)