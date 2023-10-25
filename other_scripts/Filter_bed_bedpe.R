######### Filtering PoN before GRIPSS

library(readr)
library(dplyr)
library(magrittr)
library(GenomicRanges)
library(rtracklayer)

# Upload the PoN_bed created by GRIDSS
bed_gr<- import("/Users/micoli/mnt/storageBig8/work/micoli/051122/PoN/gridss_pon_single_breakend.bed", "bed")

#Filter variants with count>1
filtered_bed <- bed_gr[bed_gr$score != 1]

#Export filtered bed
export(filtered_bed, "/Users/micoli/mnt/storageBig8/work/micoli/051122/PoN/gridss_pon_single_breakend_filtered.bed", "bed")

#Upload the PoN_bedpe created by GRIDSS
bedpe_gr<- import("/Users/micoli/mnt/storageBig8/work/micoli/051122/PoN/gridss_pon_breakpoint.bedpe", "bedpe")

#Filter variants with count>1
filtered_bedpe <- bedpe_gr[bedpe_gr@elementMetadata@listData$score !=1,]

#Export filtered bedpe
rtracklayer::export(filtered_bedpe, con="/Users/micoli/mnt/storageBig8/work/micoli/051122/PoN/gridss_pon_breakpoint_filtered.bedpe")
