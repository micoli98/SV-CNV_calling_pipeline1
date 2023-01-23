library(dplyr)
library(magrittr)
library(readr)
library(tibble)
library(tidyverse)

arg <- commandArgs(trailingOnly = TRUE)
#arg <- "/Users/micoli/mnt/storageBig8/resources/processed_data/HERCULES/WGSbams/sample_info_extended.csv"
#proj_dir <- "/Users/micoli/mnt/storageBig8/work/micoli/051122/"

samples <- read_tsv(arg[1], show_col_types = FALSE)
samples$normalSample<-replace_na(samples$normalSample, "unknown")
proj_dir <- arg[2]

# File for all tools
ids_table <- samples %>%
  filter(normal!=TRUE) %>%
  select(sample, normalSample, bamFile, normalBamFile, patient) %>%
  drop_na()
write_tsv(ids_table, "ids_table.tsv")

# Prepare file for GRIDSS: normal and all tumor bams related and put in line
# Select the non-processed samples
normal_samples <- samples[samples$normal==TRUE,]
gridss_input <- normal_samples[,c(1, 2, 13)]

gridss_input$tumorBam <- ""
for (i in gridss_input$sample) {
  j <- NULL
  j <- as.data.frame(cbind(j, ifelse(samples$normalSample == i & samples$normal == "FALSE", samples$bamFile, NA))) %>%
    na.omit(j)
  j<- toString(c(j$V1[1:nrow(j)]))
  j<- str_remove_all(j, ",")
  
  gridss_input$tumorBam <- ifelse(gridss_input$sample == i, j, gridss_input$tumorBam)
}

gridss_input <- gridss_input[order(gridss_input$sample),] %>% as.data.frame()
gridss_input$tumorBam <- str_replace(gridss_input$tumorBam, "^NA$", "")

#The samples which have a different normal (different platforms) but belong to the same patient are joined together to 
#increase the power of the joint calling
duplicato<- gridss_input[duplicated(gridss_input$patient),] #Other platform samples
true_normal <- gridss_input[!duplicated(gridss_input$patient),]  #Normal BDNA 
duplicato <- duplicato %>% unite("Bams", bamFile:tumorBam, sep=" ", remove = FALSE) %>% select(patient, Bams)

#If there are patients with multiple normal Bams they are joined to the true normal BDNA

dup_merged <- duplicato %>% 
  group_by(patient) %>% 
  summarise(Bams_all=paste(Bams, collapse=" "))

joined <- full_join(true_normal, dup_merged)
joined$Bams_all <- replace_na(joined$Bams_all, "")
joined <- joined %>% unite("all_bams", bamFile:Bams_all, sep=" ", remove = FALSE) %>% select(sample, patient, all_bams)
colnames(joined)<- c("normalSample", "patient", "Bams")
write_tsv(joined, "gridss_input.tsv", na = "")


