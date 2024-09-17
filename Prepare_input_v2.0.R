################################################# Preprocessing #######################################################
### From the file sample_info_extended, select patients to process and create the input files for the various tools

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
# patients_to_add <- "~/mnt/storageBig8/work/micoli/SCNA_Purple/resources/trial_patients.txt"
# new_samples <- read.table(new_sample_info, sep="\t", header=T)
# old_segmentation <- read.table("~/mnt/storageBig8/work/micoli/SCNA_Purple/results/231006/segmentation_info_final.tsv", sep="\t", header=T)

new_samples <- read.table(arg[1], sep="\t", header=T)
#old_samples <- read.table(arg[2], sep="\t", header=T)
old_segmentation <- read.table(arg[2], sep="\t", header=T)

#Patient selection
#Selection of patients to do comparing the sample_info_extended with the segmentation from the previous run
#Additional patients are included in the list to do according to external directions

# Filter non-usable/contaminated samples
#Filtering out:
#-non usable sample
#-contaminated samples
#-not fresh-frozen samples
#-normalSample is not NA
sample_info <- new_samples %>%  
  filter(usable==T) %>%
  filter(contamFilter == F) %>%
  filter(sampleType=="fresh frozen" | sampleType=="blood") %>%
  filter(! is.na(variantCallingNormalSample))

# Select the samples that have been added to already-processes patients and completely new patients
added_samples <- sample_info %>% 
  filter(!normal) %>% 
  filter(!sample %in% old_segmentation$sample) 

# Select the patients from previous selection and force to redo patients in the file
additional_patients <- readLines(patients_to_add)
final_sample_info <- sample_info %>% filter(patient %in% added_samples$patient | patient %in% additional_patients)

# File for all tools
ids_table <- final_sample_info %>%
  filter(normal!=TRUE) %>%
  select(sample, variantCallingNormalSample, bamFile, variantCallingNormalBamFile, patient)
colnames(ids_table) <- c("sample", "normalSample", "bamFile", "normalBamFile", "patient")
write_tsv(ids_table, "ids_table.tsv")

# Prepare file for GRIDSS: normal and all tumor bams related and put in line
# Select the non-processed samples
gridss_input <- final_sample_info %>%
  filter(normal == TRUE) %>%
  select(variantCallingNormalSample, patient, variantCallingNormalBamFile) %>% 
  unique() %>%
  mutate(tumorBam = sapply(variantCallingNormalSample, function(i) {
    final_sample_info %>%
      filter(variantCallingNormalSample == i & sample != variantCallingNormalSample) %>% #excluyde the normal samples and include the duplicated from diff platforms
      pull(bamFile) %>%
      toString() %>% #condense the lsit in a single string
      str_remove_all(",") #remove the commas
  })) %>% 
  setNames(c("normalSample", "patient", "bamFile", "tumorBam")) %>%  # Rename columns
  arrange(normalSample) %>%  # Sort by normalSample
  mutate(tumorBam = str_replace(tumorBam, "^NA$", "")) %>%  # Replace "NA" with empty string
  unite("bams", bamFile:tumorBam, sep = " ", remove = FALSE) %>%  # Combine bamFile and tumorBam
  select(normalSample, patient, bams) 

write_tsv(gridss_input, "gridss_input.tsv", na = "")

#Patients already done should have their folder copied to the new location
patients_to_move <- old_segmentation %>% 
  filter(!patient %in% gridss_input$patient) %>%
  pull(patient) %>%
  unique()

write_tsv(patients_to_move, "patients_to_move.tsv")


