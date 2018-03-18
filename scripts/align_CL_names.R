library(plyr)
library(dplyr)
library(magrittr)
library(readr)
library(taigr)

load('~/CPDS/data/CCLE/Annotations.RData')

lineage_info <- load.from.taiga(data.name='lineage-f95f', data.version=5) %>% 
  dplyr::rename(CCLE_ID = CCLE_name)
Marc_CL_data <- load.from.taiga(data.name='demeter2-marcotte-a703', data.version=5, data.file='CL_data_comb')
DRIVE_CL_data <- load.from.taiga(data.name='demeter2-drive-0591', data.version=6, data.file='CL_data_comb')
Ach_CL_data <- load.from.taiga(data.name='demeter2-achilles-5386', data.version=9, data.file='CL_data_comb')
Comb_CL_data <- load.from.taiga(data.name='demeter2-combined-dc9c', data.version=5, data.file='CL_data_comb')

name_map <- data.frame(old_name = setdiff(rownames(DRIVE_CL_data), lineage_info$CCLE_ID)) %>% 
  mutate(match_name = toupper(old_name),
         CCLE_ID = CleanCellLineName(match_name))

Novartis_sample_info <- read.csv('~/CPDS/data/Novartis/Novartis_Sample_info.csv', 
                                 check.names = FALSE, 
                                 stringsAsFactors = FALSE) %>% 
  dplyr::mutate(CELLLINE = toupper(CELLLINE), CLEANNAME_PRIMARYSITE = toupper(CLEANNAME_PRIMARYSITE)) %>% 
  dplyr::mutate(match_name = ifelse(CELLLINE %in% name_map$match_name, CELLLINE, CLEANNAME_PRIMARYSITE))

name_map %<>% left_join(
  Novartis_sample_info, by = 'match_name'
) %>% dplyr::select(old_name, new_name = CLEANNAME_PRIMARYSITE) %>% 
  mutate(new_name = revalue(new_name, c(`GISTT1_GASTROINTESTINAL_TRACT_(SITE_INDETERMINATE)` = 'GISTT1_GASTROINTESTINAL_TRACT')))

write.csv(name_map, '~/CPDS/demeter2/results/name_change_map.csv', row.names = FALSE)
