##INFO=<ID=ASRP,Number=1,Type=Integer,Description="Count of read pairs incorporated into any breakpoint assembly">                                                                                               │
##INFO=<ID=ASSR,Number=1,Type=Integer,Description="Count of split, soft clipped or indel-containing reads incorporated into any breakpoint assemblies">                                                          │
##INFO=<ID=BASRP,Number=1,Type=Integer,Description="Count of read pairs incorporated into any breakend assembly">                                                                                                │
##INFO=<ID=BASSR,Number=1,Type=Integer,Description="Count of split, soft clipped or indel-containing reads incorporated into any breakend assemblies">                                                           │
##INFO=<ID=BSC,Number=1,Type=Integer,Description="Count of soft clips supporting just local breakend">                                                                                                           │
##INFO=<ID=BUM,Number=1,Type=Integer,Description="Count of read pairs (with one read unmapped) supporting just local breakend">                                                                                  │
##INFO=<ID=IC,Number=1,Type=Integer,Description="Count of read indels supporting breakpoint">                                                                                                                    │
##INFO=<ID=REF,Number=1,Type=Integer,Description="Count of reads mapping across this breakend">                                                                                                                  │
##INFO=<ID=RP,Number=1,Type=Integer,Description="Count of read pairs supporting breakpoint">                                                                                                                     │
##INFO=<ID=SR,Number=1,Type=Integer,Description="Count of split reads supporting breakpoint">                                                                                                                    │

#### Retrieve read support to breakpoints/breakends
library(StructuralVariantAnnotation)
library(VariantAnnotation)
library(dplyr)

#Load data 
base <- "/home/micoli/mnt/storageBig8/work/micoli/051122"
sample_dirs <- list.dirs(base, full.names = F, recursive = F)
unwanted <- c("multiploidy", "linx_trial", "PoN", "sunrises", "surnrises")
sample_dirs <- sample_dirs[!(sample_dirs %in% unwanted)]

bpe_readcounts <- data.frame()
for (p in 179:length(sample_dirs))
{
  sample_files <- list.files(paste0(base, "/", sample_dirs[p]), full.names = T)
  sample_files <- sample_files[grep(".purple.sv.vcf.gz$", sample_files)]
  
  if (length(sample_files) == 0) {next}
  
  for(s in 1:length(sample_files))
  {
    sample_i <- gsub(paste0(base, "/", sample_dirs[p], "/"), "", sample_files[s])
    sample_i <- gsub(".purple.sv.vcf.gz", "", sample_i)
    sv_file <- readVcf(sample_files[s])
    
    #from the info dataframe, select only the fields wanted
    count_fields <- c("ASRP", "ASSR", "BASRP", "BASSR", "BSC", "BUM", "IC", "REF", "RP", "SR")
    count_df <- as.data.frame(info(sv_file))
    count_df <- count_df[, count_fields]
    
    #filter only the breakpoints pairs
    bpgr <- breakpointRanges(sv_file)
    bpbedpe <- breakpointgr2bedpe(bpgr)
    count_df_bp <- count_df[bpbedpe$name, ]
    if(nrow(count_df_bp)>0) {count_df_bp$Type <- "breakpoint"}
    
    #filter only breakpoints 
    begr <- breakendRanges(sv_file)
    count_df_be <- count_df[names(begr),]
    count_df_be <- count_df_be[!is.na(rowSums(count_df_be)),]
    if(nrow(count_df_be)>0) {count_df_be$Type <- "breakend"}
    
    #Join, add columns and reorder
    count_df_all <- rbind(count_df_bp, count_df_be)
    count_df_all$Name <- rownames(count_df_all)
    rownames(count_df_all) <- NULL
    count_df_all$Sample <- sample_i
    count_df_all$Patient <- sample_dirs[p]
    count_df_all <- count_df_all %>% dplyr::select(Sample, Patient, Type, Name, 1:10)
    
    bpe_readcounts <- rbind(bpe_readcounts, count_df_all)
  }
}
write.table(bpe_readcounts, paste0(base, "/breakpoint-breakends_readcounts.tsv"), sep="\t", col.names =T)

#Group gene stuff

genes_df <- data.frame()

for(p in 124:length(sample_dirs))
{
  sample_files <- list.files(paste0(base, "/", sample_dirs[p]), full.names = T)
  sample_files <- sample_files[grep(".cnv.gene.tsv$", sample_files)]
  if(length(sample_files)==0) {next}
  for(s in 1:length(sample_files))
  {
    sample_i <- gsub(paste0(base, "/", sample_dirs[p], "/"), "", sample_files[s])
    sample_i <- gsub(".purple.cnv.gene.tsv", "", sample_i)
    gene_data <- read.table(sample_files[s], sep="\t", header=T)
    gene_data$sample <- sample_i
    gene_data$patient <- sample_dirs[p]
    gene_data <- gene_data %>% dplyr::select(sample, patient, 1:17)

    genes_df <- rbind(genes_df, gene_data)
  }
}
#arrived at p 224 s null
write.table(genes_df, paste0(base, "/genes_data.tsv"), sep="\t", col.names = T)
