library(dplyr)
library(magrittr)
library(readr)
library(tibble)
library(tidyr)
library(tidyverse)

arg <- commandArgs(trailingOnly = TRUE)
#arg <- "/Users/micoli/mnt/storageBig8/resources/processed_data/HERCULES/WGSbams/sample_info_extended.csv"
#proj_dir <- "/Users/micoli/mnt/storageBig8/work/micoli"

samples <- read_tsv(arg[1], show_col_types = FALSE)
samples$normalSample<-replace_na(samples$normalSample, "unknown")
proj_dir <- arg[2]

# File for all tools
ids_table <- samples %>%
  filter(normal!=TRUE) %>%
  select(sample, normalSample, bamFile, normalBamFile, patient) %>%
  drop_na()

## File for Cobalt
#ids_table_cobalt <- data.frame(matrix(ncol=ncol(ids_table), nrow = 0))
#colnames(ids_table_cobalt) <- colnames(ids_table)
#for (i in 1:nrow(ids_table))
#{
#  sample <- ids_table[[i, 1]]
#  patient <- ids_table[[i, 5]]
#  path <- paste0(proj_dir,"/", patient, "/", sample, ".cobalt.ratio.pcf")
#  if (file.exists(path)==F)
#  {
#    ids_table_cobalt <- rbind(ids_table[i, ], ids_table_cobalt)
#  }
#}
#write_tsv(ids_table_cobalt, "ids_table_cobalt.tsv", na = "")

## File for Amber
#ids_table_amber <- data.frame(matrix(ncol=ncol(ids_table), nrow = 0))
#colnames(ids_table_amber) <- colnames(ids_table)
#for (i in 1:nrow(ids_table))
#{
#  sample <- ids_table[[i, 1]]
#  patient <- ids_table[[i, 5]]
#  path <- paste0(proj_dir,"/", patient, "/", sample, ".amber.baf.pcf")
#  if (file.exists(path)==F)
#  {
#    ids_table_amber <- rbind(ids_table[i, ], ids_table_amber)
#  }
#}
#write_tsv(ids_table_amber, "ids_table_amber.tsv", na = "")

## File for Gripss
ids_table_gripss <- data.frame(matrix(ncol=ncol(ids_table), nrow = 0))
colnames(ids_table_gripss) <- colnames(ids_table)
for (i in 1:nrow(other_prob))
{
  sample <- other_prob[[i, 1]]
  patient <- other_prob[[i, 5]]
  path <- paste0(proj_dir,"/", patient, "/", sample, ".gripss.vcf.gz")
  if (file.exists(path)==F)
  {
    ids_table_gripss <- rbind(other_prob[i, ], ids_table_gripss)
  }
}
write_tsv(ids_table_gripss, "ids_table_gripss.tsv", na = "")

## File for Purple
#ids_table_purple <- ids_table#data.frame(matrix(ncol=ncol(ids_table), nrow = 0))
#colnames(ids_table_purple) <- colnames(ids_table)
#for (i in 1:nrow(ids_table))
#{
#  sample <- ids_table[[i, 1]]
#  patient <- ids_table[[i, 5]]
#  path <- paste0(proj_dir,"/", patient, "/", sample, ".purple.cnv.somatic.tsv")
#  if (file.exists(path)==F)
#  {
#    ids_table_purple <- rbind(ids_table[i, ], ids_table_purple)
#  }
#}


## File for GRIDSS
#normal_samples <- samples[samples$normal==TRUE,]
#gridss_input <- normal_samples[,c(1, 2, 13)]

#gridss_input$tumor <- ""
#for (i in gridss_input$sample) {
#  j <- NULL
#  j <- as.data.frame(cbind(j, ifelse(samples$normalSample == i & samples$normal == "FALSE", samples$sample, NA))) %>%
#    na.omit(j)
#  j<- toString(c(j$V1[1:nrow(j)]))
#  j<- str_remove_all(j, ",")
  
#  gridss_input$tumor <- ifelse(gridss_input$sample == i, j, gridss_input$tumor)
#}

#gridss_input <- gridss_input[order(gridss_input$sample),] 
#gridss_input <- as.data.frame(gridss_input)
#gridss_input$tumor <- str_replace(gridss_input$tumor, "^NA$", "")

#gridss_input_2 <- gridss_input#data.frame(matrix(ncol=4))
#colnames(gridss_input_2) <- colnames(gridss_input)

#for (i in 1:nrow(gridss_input))
#{
#  sample <- gridss_input[[i, 1]]
#  patient <- gridss_input[[i, 2]]
#  path <- paste0(proj_dir,"/", patient, "/", sample, "_calls.vcf")
#  if (file.exists(path)==F)
#  {
#    gridss_input_2 <- rbind(gridss_input[i, ], gridss_input_2)
#  }
#}

#The samples which have a different normal (different platforms) but belong to the same patient are joined together to 
#increase the power of the joint calling
#duplicato<- gridss_input[duplicated(gridss_input$patient),] #Other platform samples
#true_normal <- gridss_input[!duplicated(gridss_input$patient),]  #Normal BDNA 
#duplicato <- duplicato %>% unite("Bams", bamFile:tumorBam, sep=" ", remove = FALSE) %>% select(patient, Bams)

#If there are patients with multiple normal Bams they are joined to the true normal BDNA
#if (nrow(duplicato > 0)) {
# true_normal$other_plat<-""
#  for (line in 1:nrow(duplicato)) {
#    pat <- duplicato[line, 1]
#    bams <- duplicato[line, 2]
#    row<-grep(pat, true_normal$patient)
#    true_normal[row, "other_plat"]<-bams
#  }
  
#  true_normal$other_plat <- replace_na(true_normal$other_plat, "")
#  true_normal <- true_normal %>% unite("all_bams", bamFile:other_plat, sep=" ", remove = FALSE) %>% select(sample, patient, all_bams)
#  colnames(true_normal)<- c("normalSample", "patient", "Bams")
#  write_tsv(true_normal, "ids_table_gridss.tsv", na = "")
  
#} else {
#  gridss_input_2 <- gridss_input_2 %>% unite("all_bams", bamFile:tumorBam, sep=" ", remove = FALSE) %>% select(sample, patient, all_bams)
#  colnames(gridss_input_2)<- c("normalSample", "patient", "Bams")
#  write_tsv(gridss_input_2, "ids_table_gridss.tsv", na = "")
#}

#true_normal <- gsub("_(bgiseq|novaseq|hiseq)$", "", ids_table_purple$normalSample)
#ids_table_purple$normalSample <- true_normal
#ids_table_purple <- ids_table_purple %>% select(sample, normalSample, patient)
#write_tsv(ids_table_purple, "ids_table_purple.tsv", na = "")


