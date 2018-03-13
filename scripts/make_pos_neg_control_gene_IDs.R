library(plyr); library(magrittr); 
library(reshape2)
library(tibble)
library(dplyr)
library(reshape2)
library(taigr);
library(ggplot2)
library(ggrepel)
library(dependr)
library(readr)
library(stringr)
library(cdsr)
library(data.table)
library(biomaRt)

#load curated pos and neg gene sets
# hart_ess <- read.csv('~/CPDS/data/gene_lists/core-essentials-291-genes-24-of-top48-screens-217constitutive', 
#                      header = FALSE, stringsAsFactors = F)$V1
# hart_non_ess <- read_tsv('~/CPDS/data/gene_lists/hart_training_nonessential.txt')$Gene

hart_info <- read.csv('~/CPDS/demeter2/data/hart_supplemental_gene_table.csv', stringsAsFactors = F, check.names=F)

sh_targets <- load.from.taiga(data.name = 'gpp-shrna-mapping-8759', data.version = 1, data.file = 'CP1175_20171102_19mer')
sh_target_gene_info <- distinct(sh_targets, `Gene Symbol`, `Gene ID`, .keep_all=T) %>% 
    dplyr::select(Gene_Symbol = `Gene Symbol`, Gene_ID = `Gene ID`)

hart_neg_controls <- hart_info$`Nonessential Genes (NE)_updated_symbols`
hart_pos_controls <- hart_info$`ConstitutiveCoreEssentials(CCE)_updated_symbols`

hart_neg_controls <- hart_neg_controls[hart_neg_controls != '']
hart_pos_controls <- hart_pos_controls[hart_pos_controls != '']

hart_pos_controls <- data.frame(Gene_Symbol = hart_pos_controls) %>% 
    left_join(sh_target_gene_info, by = 'Gene_Symbol')
hart_neg_controls <- data.frame(Gene_Symbol = hart_neg_controls) %>% 
    left_join(sh_target_gene_info, by = 'Gene_Symbol')

manual_maps <- c(DUX4L7 = '653543', LINC01599 = '196913', C8orf17 = '100507249')
hart_neg_controls[hart_neg_controls$Gene_Symbol %in% names(manual_maps), 'Gene_ID'] <- manual_maps[hart_neg_controls[hart_neg_controls$Gene_Symbol %in% names(manual_maps), 'Gene_Symbol']]

write.csv(hart_pos_controls, '~/CPDS/demeter2/data/hart_pos_controls.csv', row.names=F)
write.csv(hart_neg_controls, '~/CPDS/demeter2/data/hart_neg_controls.csv', row.names=F)
