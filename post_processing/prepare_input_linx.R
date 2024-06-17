library(dplyr)

# Input: sample_info_extended file 
arg <- commandArgs(trailingOnly = TRUE)
#sample_data <- read_tsv(arg[1])

segment_path <- arg[1]
segment_path <- "~/mnt/storageBig8/work/micoli/SCNA_Purple/results/231006/segmentation_info_final.tsv"
segments <- read.table(segment_path, sep="\t", header = T)

ids <- segments %>% select(sample, patient) %>% unique()
# # Filtering and selecting
# samples_filtered <- sample_data %>%
#   filter(usable==T) %>%
#   filter(contamFilter == F) %>%
#   filter(sampleType=="fresh frozen" | sampleType=="blood") %>%
#   filter(normal==F) %>%
#   select(sample, variantCallingNormalSample, patient, bamFile, variantCallingNormalBamFile) %>%
#   rename("normalSample" = "variantCallingNormalSample",
#          "normalBamFile"="variantCallingNormalBamFile")
# 
# samples_filtered <- samples_filtered[!is.na(samples_filtered$normalSample),] 

write_tsv(samples_filtered, "ids.tsv")
