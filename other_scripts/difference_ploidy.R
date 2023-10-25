library(magrittr)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(readxl)

#### Candidates table 
#arg <- commandArgs(trailingOnly = TRUE)
arg <- "/Users/micoli/mnt/storageBig8/resources/processed_data/HERCULES/WGSbams/sample_info_extended.csv"
sample_info<-read_tsv(arg, show_col_types = F )
sample_info<- sample_info %>%
  filter(normal==FALSE, sampleType!=c("organoid", "CL")) %>%
  dplyr::select(sample, patient)

publishDir <- "/Users/micoli/mnt/storageBig8/work/micoli/pipeline/output/"
ascat_path <- read_tsv("/Users/micoli/mnt/storageBig8/work/micoli/combinedAscatEstimates.csv")

sample_info$ploidyPurple <- 0
sample_info$ploidyAscat <- 0
sample_info$label <- ""
for (i in 1:nrow(sample_info)) {
  sample <- sample_info[[i, 1]]
  patient <- sample_info[[i, 2]]
  pur_pl_path <- paste0("/Users/micoli/mnt/storageBig8/work/micoli/pipeline/Scripts_Dependencies/", patient, "/", patient, "_purple/", sample, ".purple.purity.tsv")
  
  if (file.exists(pur_pl_path) == TRUE) {
    pur_pl <-read_tsv(paste0("/Users/micoli/mnt/storageBig8/work/micoli/pipeline/Scripts_Dependencies/", patient, "/", patient, "_purple/", sample, ".purple.purity.tsv"), show_col_types = FALSE)
    sample_info[i, 3] <- pur_pl[[1,5]]
  }
  if (sum(sample==ascat_path$sample)!=0) {
    ascat_pl <- ascat_path[sample==ascat_path$sample, 5][[1]]
    sample_info[i, 4] <- ascat_pl
  }
}

sample_info <- mutate(sample_info, diff = abs(ploidyPurple - ploidyAscat))
sample_info <- sample_info[order(sample_info$diff, decreasing = T),]
sample_info <- filter(sample_info, ploidyPurple!=0, ploidyAscat!=0)
analyse <- sample_info[1:10,]

write_tsv(sample_info, paste0(publishDir, "ploidy_difference.tsv"))
write_tsv(analyse, paste0(publishDir, "big_diff.tsv"))

######################### SCATTER PLOT

sample_info$label <- ""

for (i in 1:nrow(sample_info)) {
  if (sample_info[i, 6] > 1) {
    sample_info[i, 5] <- sample_info[i, 1]
  }
}

jitter <- position_jitter(width = 0.1, height = 0.1)

plot <- ggplot(sample_info, aes(x=ploidyPurple, y=ploidyAscat, label=label)) + 
  geom_point(position = jitter, size=1, color = ifelse(sample_info$label == "", "grey50", "red")) +
  geom_text_repel(box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf, size = 2, segment.size = 0.2, segment.alpha = 5.0) +
  geom_abline(slope = 1, intercept = 0, na.rm = FALSE, show.legend = NA, color = "blue") 
plot
ggsave("ploidy_comparison.png", path = publishDir, plot = plot, width = 40, height = 25, units = "cm")



