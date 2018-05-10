library(plyr)
library(dplyr)
library(magrittr)
library(readr)
library(taigr)

new_name_table <- read.csv('~/CPDS/demeter2/results/name_change_map.csv', check.names = FALSE, stringsAsFactors = FALSE)
new_name_map <- new_name_table$new_name %>% set_names(new_name_table$old_name)

CCLE_name_correction = read.csv('~/CPDS/demeter2/data/CCLE_name_corrections.csv', check.names = F, stringsAsFactors = F)
CCLE_name_map <- CCLE_name_correction$new_name %>% set_names(CCLE_name_correction$old_name)

#---------------------------------------------------------------------------------------------------------
#Achilles LFC data
#---------------------------------------------------------------------------------------------------------
#Achilles 98k
cur_data <- load.from.taiga(data.name = 'achilles-98k-repcollapsed-lfc-19ce', data.version = 1)
write.csv(cur_data, '~/CPDS/data/D2_figshare/achilles-98k-repcollapsed-lfc.csv')

#55k batch 1
cur_data <- load.from.taiga(data.name = 'achilles-55k-batch1-repcollapsed-lfc-d708', data.version = 1)
write.csv(cur_data, '~/CPDS/data/D2_figshare/achilles-55k-batch1-repcollapsed-lfc.csv')

#55k batch 2
cur_data <- load.from.taiga(data.name = 'achilles-55k-batch2-repcollapsed-lfc-bd7f', data.version = 1)
write.csv(cur_data, '~/CPDS/data/D2_figshare/achilles-55k-batch2-repcollapsed-lfc.csv')


#---------------------------------------------------------------------------------------------------------
#DRIVE LFC data
#---------------------------------------------------------------------------------------------------------
#DRIVE BFPD
cur_data <- load.from.taiga(data.name = 'drive-lfc-matrices-3867', data.version = 11, data.file = 'BGPD_LFC_mat')
write.csv(cur_data, '~/CPDS/data/D2_figshare/drive-bgpd-lfc-mat.csv')

#DRIVE poolA
cur_data <- load.from.taiga(data.name = 'drive-lfc-matrices-3867', data.version = 11, data.file = 'poolA_LFC_mat')
write.csv(cur_data, '~/CPDS/data/D2_figshare/drive-poola-lfc-mat.csv')

#DRIVE poolB
cur_data <- load.from.taiga(data.name = 'drive-lfc-matrices-3867', data.version = 11, data.file = 'poolB_LFC_mat')
write.csv(cur_data, '~/CPDS/data/D2_figshare/drive-poolb-lfc-mat.csv')


#---------------------------------------------------------------------------------------------------------
#Marcotte LFC data
#---------------------------------------------------------------------------------------------------------
cur_data <- load.from.taiga(data.name = 'marcotte-lfc-data-ccle-id-d459', data.version =2)
write.csv(cur_data, '~/CPDS/data/D2_figshare/Marcotte_LFC_matrix.csv')


#---------------------------------------------------------------------------------------------------------
#shRNA mapping
#---------------------------------------------------------------------------------------------------------
cur_data <- load.from.taiga(data.name='gpp-shrna-mapping-8759', data.version=6, data.file='shmap_19mer_noXLOC_Entrezonly')
write.csv(cur_data, '~/CPDS/data/D2_figshare/shRNA-mapping.csv', row.names = FALSE)


#---------------------------------------------------------------------------------------------------------
#Pos neg controls
#---------------------------------------------------------------------------------------------------------
#pos cons
cur_data <- load.from.taiga(data.name='demeter2-pos-neg-controls-a5c6', data.version=1, data.file='hart_pos_controls')
write.csv(cur_data, '~/CPDS/data/D2_figshare/Hart-pos-controls.csv', row.names = FALSE)

#neg cons
cur_data <- load.from.taiga(data.name='demeter2-pos-neg-controls-a5c6', data.version=1, data.file='hart_neg_controls')
write.csv(cur_data, '~/CPDS/data/D2_figshare/Hart-neg-controls.csv', row.names = FALSE)


#---------------------------------------------------------------------------------------------------------
#D2 Achilles model fit
#---------------------------------------------------------------------------------------------------------
#CL data
cur_data <- load.from.taiga(data.name='demeter2-achilles-5386', data.version=10, data.file='CL_data_comb')
write.csv(cur_data, '~/CPDS/data/D2_figshare/D2_Achilles_CL_data.csv')

#hp data
cur_data <- load.from.taiga(data.name='demeter2-achilles-5386', data.version=10, data.file='hp_data_comb')
write.csv(cur_data, '~/CPDS/data/D2_figshare/D2_Achilles_hp_data.csv')

#gene means
cur_data <- load.from.taiga(data.name='demeter2-achilles-5386', data.version=10, data.file='gene_means_proc')
write.csv(cur_data, '~/CPDS/data/D2_figshare/D2_Achilles_gene_dep_scores.csv')

#gene SDs
cur_data <- load.from.taiga(data.name='demeter2-achilles-5386', data.version=10, data.file='gene_SDs_proc')
write.csv(cur_data, '~/CPDS/data/D2_figshare/D2_Achilles_gene_dep_score_SDs.csv')



#---------------------------------------------------------------------------------------------------------
#D2 DRIVE model fit
#---------------------------------------------------------------------------------------------------------
#CL data
cur_data <- load.from.taiga(data.name='demeter2-drive-0591', data.version=8, data.file='CL_data_comb')
write.csv(cur_data, '~/CPDS/data/D2_figshare/D2_DRIVE_CL_data.csv')

#hp data
cur_data <- load.from.taiga(data.name='demeter2-drive-0591', data.version=8, data.file='hp_data_comb')
write.csv(cur_data, '~/CPDS/data/D2_figshare/D2_DRIVE_hp_data.csv')

#gene means
cur_data <- load.from.taiga(data.name='demeter2-drive-0591', data.version=8, data.file='gene_means_proc')
write.csv(cur_data, '~/CPDS/data/D2_figshare/D2_DRIVE_gene_dep_scores.csv')

#gene SDs
cur_data <- load.from.taiga(data.name='demeter2-drive-0591', data.version=8, data.file='gene_SDs_proc')
write.csv(cur_data, '~/CPDS/data/D2_figshare/D2_DRIVE_gene_dep_score_SDs.csv')



#---------------------------------------------------------------------------------------------------------
#D2 Combined model fit
#---------------------------------------------------------------------------------------------------------
#CL data
cur_data <- load.from.taiga(data.name='demeter2-combined-dc9c', data.version=8, data.file='CL_data_comb')
write.csv(cur_data, '~/CPDS/data/D2_figshare/D2_combined_CL_data.csv')

#hp data
cur_data <- load.from.taiga(data.name='demeter2-combined-dc9c', data.version=8, data.file='hp_data_comb')
write.csv(cur_data, '~/CPDS/data/D2_figshare/D2_combined_hp_data.csv')

#gene means
cur_data <- load.from.taiga(data.name='demeter2-combined-dc9c', data.version=8, data.file='gene_means_proc')
write.csv(cur_data, '~/CPDS/data/D2_figshare/D2_combined_gene_dep_scores.csv')

#gene SDs
cur_data <- load.from.taiga(data.name='demeter2-combined-dc9c', data.version=8, data.file='gene_SDs_proc')
write.csv(cur_data, '~/CPDS/data/D2_figshare/D2_combined_gene_dep_score_SDs.csv')


#---------------------------------------------------------------------------------------------------------
#WES_SNP_CN_data
#---------------------------------------------------------------------------------------------------------
cur_data <- load.from.taiga(data.name = 'depmap-wes-cn-data-97cc', data.version = 5, data.file = 'public_18q1_gene_cn', transpose = TRUE)
D2_comb_CLs <- load.from.taiga(data.name='demeter2-combined-dc9c', data.version=6, data.file='CL_data_comb')
cur_data <- cur_data[, rownames(D2_comb_CLs)[rownames(D2_comb_CLs) %in% colnames(cur_data)]]
write.csv(cur_data, '~/CPDS/data/D2_figshare/WES_SNP_CN_data.csv')


