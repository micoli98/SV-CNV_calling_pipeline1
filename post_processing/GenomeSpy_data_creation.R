###################### Create GenomeSpy segmentation data
library(dplyr)
library(readr)

segments_new <- read.table("~/mnt/storageBig8/work/micoli/SCNA_Purple/results/231006/segmentation_info_final.tsv", sep="\t", header = T)
meta_new <- segments_new %>% select(sample, patient, ploidy, purity, non_curated_purity) %>% unique()

#Add tissue 
tissue <- gsub("^[A-Z]+[0-9]+_", "", meta_new$sample)
tissue <- gsub('_.*', "", tissue)
tissue <- substring(tissue, 2)
tissue <- gsub('[0-9]+', '', tissue)
tissue <- gsub("(L|R)$", "", tissue)
meta_new$tissue <- tissue

#sub("^(\\w+\\d+)_([piro]\\d?)(\\w+)(\\d?)_(.*)$", "\\3", meta_new$sample) work in progress

#Add phase
meta_new$phase <- sub( '^(\\w+\\d+)_([piro])(.*)', '\\2', meta_new$sample)

#Write tables in GenomeSpy folder
write_tsv(segments_new, "~/mnt/storageBig8/web-intra/home/micoli/GenomeSpy/segmentation/data/segments_sets1-15.tsv")
write_tsv(meta_new, "~/mnt/storageBig8/web-intra/home/micoli/GenomeSpy/segmentation/data/meta_sets1-15.tsv")
