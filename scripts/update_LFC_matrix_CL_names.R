library(plyr)
library(dplyr)
library(magrittr)
library(readr)
library(taigr)

new_name_table <- read.csv('~/CPDS/demeter2/results/name_change_map.csv', check.names = FALSE, stringsAsFactors = FALSE)
new_name_map <- new_name_table$new_name %>% set_names(new_name_table$old_name)

CCLE_name_correction = read.csv('~/CPDS/demeter2/data/CCLE_name_corrections.csv', check.names = F, stringsAsFactors = F)
CCLE_name_map <- CCLE_name_correction$new_name %>% set_names(CCLE_name_correction$old_name)

#DRIVE BFPD
cur_data <- load.from.taiga(data.name = 'drive-lfc-matrices-3867', data.version = 9, data.file = 'BGPD_LFC_mat')
colnames(cur_data) %<>% plyr::revalue(new_name_map) %>% plyr::revalue(CCLE_name_map)
write.csv(cur_data, '~/Desktop/BGPD_LFC_mat.csv')

#DRIVE poolA
cur_data <- load.from.taiga(data.name = 'drive-lfc-matrices-3867', data.version = 9, data.file = 'poolA_LFC_mat')
colnames(cur_data) %<>% plyr::revalue(new_name_map) %>% plyr::revalue(CCLE_name_map)
write.csv(cur_data, '~/Desktop/poolA_LFC_mat.csv')

#DRIVE poolB
cur_data <- load.from.taiga(data.name = 'drive-lfc-matrices-3867', data.version = 9, data.file = 'poolB_LFC_mat')
colnames(cur_data) %<>% plyr::revalue(new_name_map) %>% plyr::revalue(CCLE_name_map)
write.csv(cur_data, '~/Desktop/poolB_LFC_mat.csv')

#Marcotte data
cur_data <- load.from.taiga(data.name = 'marcotte-lfc-data-ccle-id-d459', data.version =1)
colnames(cur_data) %<>% plyr::revalue(new_name_map) %>% plyr::revalue(CCLE_name_map)
write.csv(cur_data, '~/Desktop/Marcotte_LFC_matrix.csv')
