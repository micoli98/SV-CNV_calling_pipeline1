###################### Solving the matching problem in Gripss when there are more BDNA #######################

library(dplyr)
library(readr)
library(tidyverse)

samples<-read_tsv("~/mnt/storageBig8/work/micoli/pipeline/Common_data/sample_info_extended_061222.csv", show_col_types = FALSE)

## Issue with normals
#Checking step already done but redo it for the future (identification of defective samples)
samples_notnorm<- read_tsv("~/mnt/storageBig8/work/micoli/temporary_stuff/gripss_notnorm.tsv")
gripss <- samples_notnorm %>% filter(gripss==0) %>% drop_na()
patients_gripss <- gripss$patient %>% unique()
gripss_samples <- samples %>% filter(patient %in% patients_gripss)
normals_gripss <- gripss$normalSample %>% unique()

normals_gripss <- as.data.frame(normals_gripss)
normals_gripss$patient <- ""
normals_gripss$samples <- ""

for (i in 1:nrow(normals_gripss))
{
  normal <- normals_gripss[[i, 1]]
  patient <- gripss_samples[gripss_samples$sample==normal, "patient"]
  
  normals_gripss[normals_gripss$normals_gripss==normal, "patient"] <- patient
    
  file <- paste0("~/mnt/storageBig8/work/micoli/051122/", patient, "/", normal, "_calls.vcf")
  vcf <- read_tsv(file, comment="##", n_max=2)
  
  samples_patient <- colnames(vcf)[-c(1:9)]
  normal_problematic <- samples_patient[grep("BDNA[0-9]*_(bgiseq|novaseq|hiseq)", samples_patient)]
  samples_patient <- samples_patient[!(samples_patient %in% normal_problematic)]
  
  list <- paste(samples_patient, collapse = ",")
  normals_gripss[normals_gripss$normals_gripss==normal, "samples"] <- list
}

write_tsv(normals_gripss, "~/mnt/storageBig8/work/micoli/temporary_stuff/gripss_to_correct.tsv")

## Issue with samples
samples_notnorm<- read_tsv("~/mnt/storageBig8/work/micoli/temporary_stuff/gripss_notnorm.tsv")
gripss <- samples_notnorm %>% filter(gripss==0) %>% drop_na()

#lista di quelli che devo rifare
gripss_samples <- samples %>% filter(sample %in% gripss$sample)

#file per correggere i vcf (uno per paziente)
#rimuovo dai pazienti incasinati tutti i bgiseq/novaseq/hiseq e festa finita

patients_gripss <- gripss %>% dplyr::select(patient, normalSample) %>% unique()
patients_gripss$samples <- ""

for (i in 1:nrow(patients_gripss))
{
  normal <- patients_gripss[[i, 2]]
  patient <- patients_gripss[[i, 1]]
  
  file <- paste0("~/mnt/storageBig8/work/micoli/051122/", patient, "/", normal, "_calls.vcf")
  vcf <- read_tsv(file, comment="##", n_max=2)
  
  samples_patient <- colnames(vcf)[-c(1:9)]
  sample_problematic <- samples_patient[grep("_(bgiseq|novaseq|hiseq)", samples_patient)]
  samples_patient <- samples_patient[!(samples_patient %in% sample_problematic)]
  
  list <- paste(samples_patient, collapse = ",")
  patients_gripss[i, "samples"] <- list
}

write_tsv(patients_gripss, "~/mnt/storageBig8/work/micoli/temporary_stuff/gripss_to_correct.tsv")
write_tsv(gripss_samples, "~/mnt/storageBig8/work/micoli/temporary_stuff/gripss_redo.tsv")
