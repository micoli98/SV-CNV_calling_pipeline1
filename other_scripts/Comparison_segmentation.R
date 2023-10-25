library(rtracklayer)
library(StructuralVariantAnnotation)
library(VariantAnnotation)
library(magrittr)
library(readr)
library(tibble)
library(tidyverse)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)

a <- read_file("/Users/micoli/mnt/storageBig8/work/micoli/signatures/Interesting_samples.txt")
a <- a %>% str_replace_all("- ", "") %>% strsplit("\n")

sample_info <- read_tsv("/Users/micoli/mnt/storageBig8/work/micoli/pipeline/Common_data/samples_done.tsv", show_col_types = F)
sample_info <- filter(sample_info, sampleType != "CL", sampleType != "organoid")
a <- unlist(a)
b <- list()
for (i in 1:length(a)) {
  name <- paste0(a[i], ".")
  real <- sample_info[str_detect(sample_info$sample, name), 1]
  b <- append(b, real[[1]])
}
b <- unlist(b)
selected <- sample_info[sample_info$sample %in% b, ]

purifyR <- function(purity, ploidy, R) (purity*ploidy*R + 2*(1-purity)*(R-1)) / (purity*ploidy)
purifyLogR <- function(purity, ploidy, logR) log2(pmax(2^-10, purifyR(purity, ploidy, 2^logR)))
home<- "/Users/micoli/mnt/storageBig8/work/micoli/pipeline/Scripts_Dependencies"

diff <- read_tsv("/Users/micoli/mnt/storageBig8/work/micoli/pipeline/output/big_diff.tsv")
diff_names <- unlist(diff$sample)
selected <- sample_info[sample_info$sample %in% diff_names, ]

############################### TABLES CREATION

for (i in 1:nrow(selected)) {
  sample_id <- selected[[i, 1]]
  patient<-selected[[i, 2]]
  path<-paste0(home, '/', patient, '/', patient, '_purple/', sample_id)
  pathCob<-paste0(home, '/', patient, '/', patient,'_cobalt/', sample_id)
  publishDir<-'/Users/micoli/mnt/storageBig8/work/micoli/pipeline/output/segmentation_comparison/'
  
  ############################## COUNTS TABLE
  
  ############### Segments by GRIDSS
  
  if (file.exists(paste0(path, '.purple.cnv.somatic.tsv')) == TRUE) {
    GRIDSS_SOM<- read.table(file = paste0(path, '.purple.cnv.somatic.tsv'), header = TRUE)
    PUR_PL <- read.table(file = paste0(path, '.purple.purity.tsv'), header = TRUE)
    PUR<-PUR_PL[[1]]
    PL<-PUR_PL[[5]]
    
    # Select interesting columns and convert copy number in logR
    
    GRIDSS_SEG <- GRIDSS_SOM %>%
      dplyr::select(chromosome, start, end, minorAlleleCopyNumber, majorAlleleCopyNumber, copyNumber, baf) %>%
      dplyr::rename(
        chr=chromosome,
        nMinor=minorAlleleCopyNumber,
        nMajor=majorAlleleCopyNumber,
      ) 
    GRIDSS_SEG[GRIDSS_SEG$copyNumber<=0, "copyNumber"]<- 0.00001
    GRIDSS_SEG <- GRIDSS_SEG %>%
      mutate(logR=log2(copyNumber/PL)) %>%
      tibble::rowid_to_column("ID")
    
    write_tsv(GRIDSS_SEG, paste0(publishDir, sample_id, "_SEG_HMF.tsv"), na = "")
    
    ############### Segments by GATK
    
    GATK_SEG <- read_tsv("/Users/micoli/mnt/storageBig8/work/klavikka/HERCULES/gatk_ascat/result_gatk_ascat/output/combinedAscatSegments-out.csv", show_col_types = FALSE)
    GATK_SEG <- GATK_SEG %>%
      subset(sample==sample_id) %>%
      dplyr::select(chr, startpos, endpos, purifiedLogR, purifiedBaf)
    
    write_tsv(GATK_SEG, paste0(publishDir, sample_id, "_SEG_GATK.tsv"), na = "")
    
    ############### Raw_counts by GATK
    
    GATK_DEN_COUNTS <- read_tsv(paste0("/Users/micoli/mnt/storageBig8/work/klavikka/HERCULES/gatk_ascat/result_gatk_ascat/", sample_id, "-denoiseCounts-bash/", sample_id, ".denoisedCR.tsv"), show_col_types = FALSE, col_names = TRUE, comment="@")
    colnames(GATK_DEN_COUNTS)<-c("CHR", "start", "end", "logR")
    
    GATK_DEN_COUNTS<- GATK_DEN_COUNTS %>%
      dplyr::select(CHR, start, end, logR) %>%
      mutate(pos=((start+end)/2)) %>%
      mutate(pur_logR=purifyLogR(PUR, PL, logR)) %>%
      dplyr::select(CHR, pos, logR, pur_logR) %>%
      mutate (type="GATK")
    
    cut_GATK_DEN_COUNTS <- GATK_DEN_COUNTS %>%
      filter(pur_logR<= -0.001 | pur_logR>=0.001) %>%
      sample_n(nrow(GATK_DEN_COUNTS)/2)
    Tot_GATK<-sum(GATK_DEN_COUNTS$pur_logR)/nrow(GATK_DEN_COUNTS)
    
    ############### Raw_counts by COBALT
    
    COBALT_COUNTS<-read_tsv(paste0(pathCob, ".cobalt.ratio.tsv"), show_col_types = FALSE, col_names = TRUE)
    COBALT_COUNTS <- COBALT_COUNTS %>%
      mutate(logR=log2(tumorGCRatio/referenceGCDiploidRatio)) %>%
      dplyr::select(chromosome, position, logR) %>%
      mutate(pur_logR=purifyLogR(PUR, PL, logR))
    colnames(COBALT_COUNTS)<- c("CHR", "pos", "logR", "pur_logR")
    
    cut_COBALT_COUNTS <- COBALT_COUNTS %>%
      filter(pur_logR<= -0.001 | pur_logR>=0.001) %>%
      sample_n(nrow(COBALT_COUNTS)/2)
    Tot_COBALT<-sum(COBALT_COUNTS$pur_logR)/nrow(COBALT_COUNTS)
    
    ############### Merging counts in a unique dataframe
    
    COUNTS <- cut_GATK_DEN_COUNTS %>% 
      full_join(cut_COBALT_COUNTS, by=c("CHR" = "CHR", "pos" = "pos", "logR"="logR", "pur_logR"="pur_logR"))
    
    COUNTS$pos <- ceiling(COUNTS$pos)
    COUNTS$type <- COUNTS$type %>% replace_na("HMF")
    COUNTS<-COUNTS[order(COUNTS$CHR, COUNTS$pos),]
    
    write_tsv(COUNTS, paste0(publishDir, sample_id, "_COUNTS.tsv"), na = "")
    
    ############################## STRUCTURAL VARIANTS TABLE
    
    SVs_vcf <- readVcf(paste0(path, ".purple.sv.vcf.gz"), "hg19")
    SVs<-as.data.frame(rowRanges(SVs_vcf))
    SVs<- dplyr::select(SVs, c(1,2,3,10))
    SVs<-mutate(SVs, pos=(start + end)/2) 
    
    # Select only passed and inferred
    final_SVs<-SVs[grep(pattern="PASS|INFERRED", SVs$FILTER),]
    final_SVs<-final_SVs %>% mutate(PASS=ifelse(final_SVs$FILTER=="PASS", TRUE, FALSE))
    
    
    # Add frequency
    id_col<-grep(sample_id, colnames(geno(SVs_vcf)$AF))
    AF<-geno(SVs_vcf)$AF[rownames(final_SVs), id_col, drop=FALSE]
    colnames(AF)<- "AF"
    AF<- as.data.frame(AF)
    SVsAF<-merge(final_SVs, AF, by=0, all=TRUE)
    SVsAF$AF<-as.numeric(SVsAF$AF) %>%
      replace_na(1)
    rownames(SVsAF)<-SVsAF$Row.names
    write_tsv(SVsAF, paste0(publishDir, sample_id, "_SVsAF.tsv"), na = "")
    
  ####################################### CREATE THE JSON PLOT
  
  json <- paste0(
    '{
      "$schema": "https://unpkg.com/@genome-spy/core/dist/schema.json",
      "genome": { "name": "hg38" },
      "resolve": { "scale": { "color": "independent" }},
      "layer": [
        {
          "name": "SVs",
          "data": {"url": "', sample_id, '_SVsAF.tsv"},
          "mark": {"type": "rule"},
          "encoding":{
            "color": {
              "field": "FILTER",
              "type":"nominal",
              "scale": {
                "domain": ["PASS", "INFERRED"],
                "range": ["grey", "#e7298a"]
              } },
            "x": {
              "type": "locus",
              "chrom": "seqnames",
              "pos": "pos"
            },
            "opacity": {"field": "AF", "type": "quantitative"}
          }
        },
        {
          "name": "Counts",
          "data": {"url": "', sample_id, '_COUNTS.tsv"},
          "mark": {"type": "point",
            "geometricZoomBound": 10
          },
          "encoding":{
            "color": {
              "field": "type",
              "type":"nominal",
              "scale": {
                "domain": ["GATK", "HMF"],
                "range": ["#1b9e77", "#d95f02"]
              } },
            "x": {
              "type": "locus",
              "chrom": "CHR",
              "pos": "pos"
            },
            "y": {
              "type": "quantitative",
              "field": "pur_logR", 
              "scale": {"domain": [-4, 4]}
            },
            "size": { "value": 100 },
            "opacity": { "value": 0.3 }
          }
        },
        {
          "name": "HMF",
          "data": {"url": "', sample_id, '_SEG_HMF.tsv"},
          "mark": {
            "type": "rule",
            "size": 3,
            "strokeCap": "butt"
          },
          "encoding":{
            "color": {"value": "#d95f02"},
            "x": {
              "type": "locus",
              "chrom": "chr",
              "pos": "start"
            },
            "y": {
              "type": "quantitative",
              "field": "logR", 
              "scale": {}
            },
            "x2": {
              "chrom": "chr",
              "pos": "end"
            }
          }  
        },
        {
          "name": "GATK",
          "data": {"url": "', sample_id, '_SEG_GATK.tsv"},
          "mark": {
            "type": "rule",
            "size": 3,
            "strokeCap": "butt"
          },
          "encoding":{
            "color": { "value": "#1b9e77" },
            "x": {
              "type": "locus",
              "chrom": "chr",
              "pos": "startpos"
            },
            "y": {
              "type": "quantitative",
              "field": "purifiedLogR", 
              "scale": {}
            },
            "x2": {
              "chrom": "chr",
              "pos": "endpos"
            }
          } 
        }
      ]
    }'
  )
  write(json, file = paste0(publishDir, sample_id, "_comparison.json"))
  }
} 







