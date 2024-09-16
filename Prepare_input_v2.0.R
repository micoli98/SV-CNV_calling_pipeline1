## TODO: there is a problem in the sample selection for the new run. If a sample is added to an already existing patient, this patient is not selected.
##       Recover the missing samples and rerun the analysis. Moreover, the moving file is also wrong. 

suppressMessages({
  library(dplyr)
  library(stringr)
  library(tidyverse)
  library(readr)})

arg <- commandArgs(trailingOnly = TRUE)
# arg <- "~/mnt/storageBig8/work/micoli/pipeline/Common_data/sample_info_test.csv"
# proj_dir <- "/Users/micoli/mnt/storageBig8/work/micoli/051122/"
# old_sample_info <- "~/mnt/storageBig8/work/micoli/SCNA_Purple/resources/sample_info_extended_231004.csv"
# new_sample_info <- "~/mnt/storageBig8/work/micoli/SCNA_Purple/resources/sample_info_extended_240617.csv"
# patients_to_add <- "~/mnt/storageBig8/work/micoli/SCNA_Purple//trial_patients.txt"

new_samples <- read.table(arg[1], sep="\t", header=T)
#old_samples <- read.table(arg[2], sep="\t", header=T)
old_segmentation <- read.table(arg[2], sep="\t", header=T)

sample_info <- new_samples %>%  
  filter(usable==T) %>%
  filter(contamFilter == F) %>%
  filter(sampleType=="fresh frozen" | sampleType=="blood") %>%
  filter(! is.na(variantCallingNormalSample))

patients_to_do <- sample_info %>% 
  filter(!normal) %>% 
  filter(!sample %in% old_segmentation$sample | patient %in% additional_patients) 
# This should do the job... try to clean all


# 
# #Patient selection
# #Selection of patients to do comparing the sample_info_extended with the previous version
# #Additional patients are included in the list to do according to external directions
# new_patients <- new_samples[!new_samples$patient %in% old_samples$patient, "patient"] %>% unique()
# 
# additional_patients <- readLines(arg[3])
# sample_info <- new_samples %>% filter(patient %in% new_patients | patient %in% additional_patients)
# 
# #patients already done should have their folder copied to the new location
# old_samples <- old_samples %>%  
#   filter(usable==T) %>%
#   filter(contamFilter == F) %>%
#   filter(sampleType=="fresh frozen" | sampleType=="blood") 
# patients_to_move <- data.frame(patient = unique(old_samples[!old_samples$patient %in% sample_info$patient, "patient"]))
# write_tsv(patients_to_move, "patients_to_move.tsv")
# 
# #Filtering step
# #Filtering out: 
# #-non usable sample
# #-contaminated samples
# #-not fresh-frozen samples
# #-normalSample is not NA
# 
# samples_filtered <- sample_info %>%
#   filter(usable==T) %>%
#   filter(contamFilter == F) %>%
#   filter(sampleType=="fresh frozen" | sampleType=="blood") 
# 
# samples_filtered <- samples_filtered[!is.na(samples_filtered$normalSample),] 
# 
# # File for all tools
# ids_table <- samples_filtered %>%
#   filter(normal!=TRUE) %>%
#   select(sample, variantCallingNormalSample, bamFile, variantCallingNormalBamFile, patient)
# 
# colnames(ids_table) <- c("sample", "normalSample", "bamFile", "normalBamFile", "patient")
# write_tsv(ids_table, "ids_table.tsv")
# 
# # Prepare file for GRIDSS: normal and all tumor bams related and put in line
# # Select the non-processed samples
# normal_samples <- samples_filtered[samples_filtered$normal==TRUE,]
# gridss_input <- normal_samples[,c(17, 2, 18)] %>% unique()
# 
# gridss_input$tumorBam <- ""
# for (i in gridss_input$variantCallingNormalSample) {
#   j <- NULL
#   #select the bam files of the samples to process for each variantCallingNormalSample (the second conditions allows to include also normals with different platform)
#   j <- as.data.frame(cbind(j, ifelse((samples_filtered$variantCallingNormalSample == i & samples_filtered$sample!=samples_filtered$variantCallingNormalSample), samples_filtered$bamFile, NA))) %>%
#     na.omit(j)
#   j<- toString(c(j$V1[1:nrow(j)]))
#   j<- str_remove_all(j, ",")
#   
#   gridss_input$tumorBam <- ifelse(gridss_input$variantCallingNormalSample == i, j, gridss_input$tumorBam)
# }
# 
# colnames(gridss_input) <- c("normalSample", "patient", "bamFile", "tumorBam")
# gridss_input <- gridss_input[order(gridss_input$normalSample),] %>% as.data.frame()
# gridss_input$tumorBam <- str_replace(gridss_input$tumorBam, "^NA$", "")
# 
# #collapse bams in a single column "Bams"
# gridss_input <- gridss_input %>% unite("bams", bamFile:tumorBam, sep=" ", remove = FALSE) %>% select(normalSample, patient, bams)
# 
# write_tsv(gridss_input, "gridss_input.tsv", na = "")
# 
# 
