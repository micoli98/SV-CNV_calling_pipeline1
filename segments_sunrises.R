# Produces the final segmentation results for HMF pipeline
# Written by Giulia 
# Sunrise plots and subclonality plots are created using https://github.com/hartwigmedical/hmftools/blob/master/purple/src/main/resources/r/somaticVariantPlots.R functions

####### INPUT:
####### 1. Sample_info_extended
####### 2. publishDir

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(magrittr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggfx)
  library(tidyverse)
})
source("/mnt/storageBig8/work/micoli/SCNA_Purple/src/functions_plots.R")
theme_set(theme_bw())

#### Candidates table
arg <- commandArgs(trailingOnly = TRUE)
sample_info<-read.table(arg[1], sep="\t", header=T)
publishDir <- arg[2]

#Filtering for usable samples
samples_filtered <- sample_info %>%
  filter(usable==T) %>%
  filter(contamFilter == F) %>%
  filter(sampleType=="fresh frozen" | sampleType=="blood") %>% 
  filter(normal == F)
samples_filtered <- samples_filtered[!is.na(samples_filtered$normalSample),] 

# Selection of sample ids and patients of tumor samples 
real_normals <- samples_filtered[, "variantCallingNormalSample"] %>% unique()

samples_filtered <- samples_filtered %>%
  filter(!(sample %in% real_normals)) %>%
  dplyr::select(sample, patient)

# Retrieve purity and ploidy values from purple results
for (i in 1:nrow(samples_filtered)) {
  file <- paste0(publishDir, "/", samples_filtered[i, "patient"], "/", samples_filtered[i, "sample"], ".purple.purity.tsv")
  if (file.exists(file)) {
    pur_pl <-read.table(file, sep="\t", header=T)
    samples_filtered[i, "ploidy"] <- pur_pl[1, "ploidy"]
    samples_filtered[i, "non_curated_purity"] <- pur_pl[1, "purity"]
    
    #Add checking of aberrancy using the purity-range method with a threshold of 0.5
    samples_filtered[i, "aberrant"] <- ifelse(pur_pl$maxPurity - pur_pl$minPurity <= 0.5, T, F)
  }
}
# Check for completed samples (not necessary)
# samples_filtered <- samples_filtered %>% mutate(sum = rowSums(.[3:4]))
# samplesReady <- subset(samples_filtered, sum!=0)

# Join all segments retrieved from PURPLE results
# Create a function to read and add sample and patient columns
read_and_annotate <- function(sample, patient, publishDir) {
  file <- paste0(publishDir, "/", patient, "/", sample, ".purple.cnv.somatic.tsv")
  if (file.exists(file)) {
    df <- read_tsv(file, show_col_types = FALSE) %>%
      mutate(sample = sample, patient = patient)
    return(df)
  } else {
    return(NULL)
  }
}

# Iterate over samples and patients
all_segs <- map2_df(samples_filtered[["sample"]], samples_filtered[["patient"]], read_and_annotate, publishDir = publishDir) %>%
  as.data.frame()

# Calculare logR and LOH 
all_segs <- all_segs %>%
  inner_join(samples_filtered) %>%  
  mutate(logR = ifelse(copyNumber<=0, -10, log(copyNumber/ploidy)),
         Loh = abs(baf - 0.5) * 2,
         purity = ifelse(aberrant == T, non_curated_purity, 0.08)) #purity correction for non-aberrant samples

all_segs <- all_segs %>% dplyr::select(sample, patient, 1:16, "ploidy", "purity", 20:23)

write_tsv(all_segs, paste0(publishDir, "/segmentation_info_final.tsv"))

##### Sunrise and subclonal plots
for (i in 1:nrow(samples_filtered)) {
  bestFitFile <- paste0(publishDir, "/", samples_filtered[i, "patient"], "/", samples_filtered[i, "sample"], ".purple.purity.tsv")
  rangeFile <- paste0(publishDir, "/", samples_filtered[i, "patient"], "/", samples_filtered[i, "sample"], ".purple.purity.range.tsv")
  if (file.exists(bestFitFile) && file.exists(rangeFile))
  {
    bestFitDF = read.table(file = bestFitFile, sep = "\t", header = T, comment.char = "!") %>% dplyr::select(purity, ploidy, score)
    rangeDF = read.table(file = rangeFile, sep = "\t", header = T, comment.char = "!") %>%
      dplyr::select(purity, ploidy, score)
    rangePlot = purity_ploidy_range_plot(bestFitDF, rangeDF)
    ggsave(filename = paste0(samples_filtered[i, "sample"], ".purity.range.png"), rangePlot, path = paste0(publishDir, "/sunrises/"), units = "in", height = 4, width = 4.8, scale = 1)
  }
  
  file <- paste0(publishDir, "/", samples_filtered[i, "patient"], "/", samples_filtered[i, "sample"], ".purple.somatic.hist.tsv")
  if (file.exists(file) == TRUE)
  {
    somaticBuckets = read.table(file, sep = "\t", header = T, numerals = "no.loss", skipNul = T)
    clonalityModel = read.table(paste0(publishDir, "/", samples_filtered[i, "patient"], "/", samples_filtered[i, "sample"], ".purple.somatic.clonality.tsv"), sep = "\t", header = T, numerals = "no.loss", skipNul = T) %>%
      mutate(isSubclonal = isSubclonal == "true", isValid = isValid == "true", peak = as.character(peak), bucketWeight = as.numeric(as.character(bucketWeight))) %>% filter(isValid)
    clonalityModelPlot = clonality_plot(somaticBuckets, clonalityModel)
    ggsave(filename = paste0(samples_filtered[i, "sample"], ".somatic.clonality.png"), clonalityModelPlot, path = paste0("~/Desktop/surnrises"), units = "in", height = 6, width = 8, scale = 1)
  }
}

