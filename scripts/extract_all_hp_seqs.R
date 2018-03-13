library(plyr); library(magrittr); 
library(reshape2)
library(tibble)
library(dplyr)
library(taigr);
library(readr)
library(stringr)

hairpin_map_DRIVE <- read.csv('~/CPDS/data/Novartis/shRNAID_seq_map.csv', stringsAsFactors = F)
hairpin_map_Ach <- read.csv('~/CPDS/data/Achilles/shRNA/full_shrna_mapping_ACH.csv', stringsAsFactors = F)
hairpin_map_Marc <- read_tsv('~/CPDS/data/Marcotte/GPL21133_GMAP-Uts520601_probe_sequences.tab.txt')

load('~/CPDS/data/Marcotte/breast_screens.eset')
hairpin_map_Marc %<>% filter(probeset_name %in% unique(fData(breast_screens)$trps_id))

DRIVE_hairpins <- unique(hairpin_map_DRIVE$SEQ)
Ach_hairpins <- unique(hairpin_map_Ach$SEQ)
Marc_hairpins <- unique(sapply(hairpin_map_Marc$probe_sequence, function(x) {
    if (nchar(x) > 21) {
        x <- str_sub(x, 1, 21)
    }
    return(x)
}))
library(Biostrings)
dna = DNAStringSet(Marc_hairpins)
Marc_hairpins_a <- complement(dna)

all_hairpins <- data.frame(dataset = 'Achilles', hairpin_seq = Ach_hairpins) %>% 
    rbind(data.frame(dataset = 'DRIVE', hairpin_seq = DRIVE_hairpins)) %>% 
    rbind(data.frame(dataset = 'Marcotte', hairpin_seq = Marc_hairpins_a))

write.csv(all_hairpins, '~/CPDS/demeter2/data/all_hairpin_table.csv', row.names=F)
