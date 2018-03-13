# Separate file for kubeque
source('~/libraries.R')
source('~/Documents/RNAi Pipeline Scripts/Pull_Data.R')
source('~/Documents/RNAi Pipeline Scripts/fit_models.R')

gene_avg_combination(all_LFC_mats, c(rep('Achilles', 3), rep('DRIVE', 3), 'Marcotte'), Achilles_DRIVE_Marcotte_RDS_fname)
