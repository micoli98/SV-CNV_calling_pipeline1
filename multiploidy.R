# Plot ploidy statistics of patients that have samples with heterogeneous ploidies
# Written by Giulia
# Function taken from ascat-algorithm.R

library(magrittr)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(readxl)

# Visualization of multiploidy patients if present
# The R script takes as input the csv table with patient information (sample_info_extended.csv), the home and publish directories. 

# It creates a multiploidy plot per patient

conclusionLevels <- c("default", "adjust", "discard")
conclusionColors <- c("black", "#00e050", "#F04030")

######### Create Candidates table

#### Functions from ASCAT

distancesToMatrix <- function(distances) {
  rho <- sort(unique(distances$purity))
  psi <- sort(unique(distances$ploidy))
  
  m <- distances %>%
    arrange(purity, ploidy) %>%
    pull(score) %>%
    matrix(ncol = length(rho),
           nrow = length(psi))
  
  colnames(m) <- rho 
  rownames(m) <- psi
  m <- t(m)
}

findValleys <- function(distanceMatrix, radius = 3) {
  d <- distanceMatrix
  r <- radius
  
  minima <- data.frame(i = numeric(), j = numeric(), value = numeric())
  
  for (i in seq_len(ncol(d))) { # psi
    for (j in seq_len(nrow(d))) { # rho
      m <- d[j, i]
      horiz <- (i-r):(i+r)
      vert <- (j-r):(j+r)
      
      seld <- d[
        ifelse(vert > 0 & vert <= nrow(d), vert, NA),
        ifelse(horiz > 0 & horiz <= ncol(d), horiz, NA)
      ]
      seld[1+r, 1+r] = NA
      
      if (min(seld, na.rm = T) > m) {
        minima <- rbind(minima, data.frame(i, j, value = m))
      }
    }
  }
  
  minima
}

findCandidates <- function(distances) {
  rho <- sort(unique(distances$purity))
  psi <- sort(unique(distances$ploidy))
  
  distancesToMatrix(distances) %>%
    findValleys() %>%
    transmute(psi = psi[i], rho = rho[j], score = value)
}

#### Candidates table 
arg <- commandArgs(trailingOnly = TRUE)
sample_info <- suppressMessages(read_tsv(arg[1], show_col_types = FALSE))
#sample_info <- read_tsv("/Users/micoli/mnt/storageBig8/work/micoli/pipeline/Common_data/H098.csv")
publishDir <- arg[2]
#publishDir <- "/Users/micoli/mnt/storageBig8/work/micoli/pipeline/"

samples_filtered<- sample_info %>%
  filter(normal==FALSE) %>%
  dplyr::select(sample, patient)

all_candidates <- as.data.frame(matrix(ncol=5, nrow=0))
for (i in 1:nrow(samples_filtered)) {
  sample <- samples_filtered[i, 1]
  patient <- samples_filtered[i, 2]
  pur_pl_path <- paste0(publishDir, "/", patient, "/", sample, ".purple.purity.range.tsv")

  if (file.exists(pur_pl_path)) {
    pur_pl <-read_tsv(paste0(publishDir, "/", patient, "/", sample, ".purple.purity.range.tsv"), show_col_types = FALSE)
    pur_pl_decision <- as.data.frame(read_tsv(paste0(publishDir, "/", patient, "/", sample, ".purple.purity.tsv"), show_col_types = FALSE))
    candidates <- findCandidates(pur_pl) %>%
      arrange(score)

    matching_row <- candidates %>%
      filter(psi == pur_pl_decision$ploidy, rho == pur_pl_decision$purity)
    if (nrow(matching_row) == 0) {
      new_row <- data.frame(psi = pur_pl_decision$ploidy, rho = pur_pl_decision$purity, score = pur_pl_decision$score)
      candidates <- bind_rows(candidates, new_row)
    }

    candidates <- candidates %>%
      mutate(
        sample = sample,
        patient = patient,
        rank = ifelse(psi == pur_pl_decision$ploidy & rho == pur_pl_decision$purity, 1, NA)
      ) %>%
      arrange(desc(score)) %>%
      mutate(rank = ifelse(is.na(rank), row_number() + 1, rank)) %>%
      transmute(sample, patient, rank, psi, rho, score)

    all_candidates <- rbind(all_candidates, candidates)
  }
}

# write.table(all_candidates, file.path(publishDir, "opt_temp/all_candidates.tsv"), sep="\t", col.names=T)

######### Aberrant samples
aberrantSamples <- all_candidates %>%
  dplyr::filter(all_candidates$rank == 1)

find_peaks <- function(d) d$x[c(F, diff(diff(d$y) >= 0) < 0, F) & d$y > 0.01]
find_ploidy_peaks <- function(ploidy) find_peaks(density(ploidy, from = 0.08, to = 8, bw = 0.2))
find_ploidy_peak_count <- function(ploidy) length(find_ploidy_peaks(ploidy))

######### Multiple ploidy patients
patientsWithMultiplePloidies <- aberrantSamples %>%
  dplyr::group_by(patient) %>%
  dplyr::summarise(ploidyPeaks = find_ploidy_peak_count(psi),
            n = n()) %>%
  ungroup() %>%
  dplyr::filter(ploidyPeaks > 1)

######## Plots for multiple ploidy patients
for (id in unique(patientsWithMultiplePloidies$patient)) {
  h <- all_candidates %>% 
    filter(patient == id) %>%
    filter(rank %in% 1:3) %>%
    mutate(cut_score = cut(score, seq(0, 1.6, 0.2), include.lowest = TRUE)) %>%
    ggplot() +
    geom_point(aes(x = psi,
                   y = sample,
                   fill = cut_score,
                   size = rank == 1,
                   stroke = rank == 1),
               shape = 21) +
    xlim(c(1, 8)) +
    scale_size_manual(values = c(4, 7), guide = "none") +
    labs(fill = "score") +
    scale_fill_brewer(drop = F) +
    ggtitle(paste0(id, "_multiploidy"))
  h
  ggsave(paste0(id, "_multiploidy.png"), path = paste0(publishDir, "/multiploidy/"), plot = h, width = 9, height = 6, dpi = 100)
}


