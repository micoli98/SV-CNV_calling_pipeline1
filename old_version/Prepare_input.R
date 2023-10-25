library(dplyr)
library(stringr)
library(tidyverse)
library(readr)

arg <- commandArgs(trailingOnly = TRUE)
#arg <- "~/mnt/storageBig8/work/micoli/pipeline/Common_data/sample_info_test.csv"
#proj_dir <- "/Users/micoli/mnt/storageBig8/work/micoli/051122/"

samples <- read.table(arg[1], sep="\t", header=T)
proj_dir <- arg[2]

#Filtering step
#Filtering out: 
#-non usable sample
#-contaminated samples
#-not fresh-frozen samples
#-normalSample is not NA

samples_filtered <- samples %>%
  filter(usable==T) %>%
  filter(contamFilter == F) %>%
  filter(sampleType=="fresh frozen" | sampleType=="blood") 

samples_filtered <- samples_filtered[!is.na(samples_filtered$normalSample),] 

# File for all tools
ids_table <- samples_filtered %>%
  filter(normal!=TRUE) %>%
  select(sample, variantCallingNormalSample, bamFile, variantCallingNormalBamFile, patient)

colnames(ids_table) <- c("sample", "normalSample", "bamFile", "normalBamFile", "patient")
write_tsv(ids_table, "ids_table.tsv")

# Prepare file for GRIDSS: normal and all tumor bams related and put in line
# Select the non-processed samples
normal_samples <- samples_filtered[samples_filtered$normal==TRUE,]
gridss_input <- normal_samples[,c(17, 2, 18)] %>% unique()

gridss_input$tumorBam <- ""
for (i in gridss_input$variantCallingNormalSample) {
  j <- NULL
  #select the bam files of the samples to process for each variantCallingNormalSample (the second conditions allows to include also normals with different platform)
  j <- as.data.frame(cbind(j, ifelse((samples_filtered$variantCallingNormalSample == i & samples_filtered$sample!=samples_filtered$variantCallingNormalSample), samples_filtered$bamFile, NA))) %>%
    na.omit(j)
  j<- toString(c(j$V1[1:nrow(j)]))
  j<- str_remove_all(j, ",")
  
  gridss_input$tumorBam <- ifelse(gridss_input$variantCallingNormalSample == i, j, gridss_input$tumorBam)
}

colnames(gridss_input) <- c("normalSample", "patient", "bamFile", "tumorBam")
gridss_input <- gridss_input[order(gridss_input$normalSample),] %>% as.data.frame()
gridss_input$tumorBam <- str_replace(gridss_input$tumorBam, "^NA$", "")

#collapse bams in a single column "Bams"
gridss_input <- gridss_input %>% unite("bams", bamFile:tumorBam, sep=" ", remove = FALSE) %>% select(normalSample, patient, bams)

write_tsv(gridss_input, "gridss_input.tsv", na = "")


