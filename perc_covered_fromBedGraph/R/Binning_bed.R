library(binr)
library("parallel")
library("foreach")
library("doParallel")
# Data wrangling
library(dplyr)

workDir <- "perc_covered_fromBedGraph"
outDir <- file.path(workDir, "output")

# Read bed file (which contains the peak list)
bed <- read.delim(file.path(workDir, "data/example.bed"), header=F)
bed <- mutate(bed, mergeID = paste(V1, V2, V3, sep="_"))

#-----------------------------------#
# Seperate each region into bins of 100
#-----------------------------------#

start_tmp <- bed$V2
end_tmp <- bed$V3
width_tmp <- (end_tmp - start_tmp) + 1


# Using multi core [Get the start position]
cl <- makeCluster(detectCores())
registerDoParallel(cl, cores = detectCores())
tmp_list_start <- list()
tmp_list_start <- foreach (i = 1:length(start_tmp), .packages = "binr") %dopar% {
  tmp <- bins(start_tmp[i]:end_tmp[i],target.bins = 100,max.breaks = 100, exact.groups = TRUE)
  s <- gsub(",.*", "", names(tmp$binct))
  s <- gsub("\\[","",s)
  tmp_list_start[[i]] <- as.numeric(s) 
}
stopCluster(cl)


# Using multi core [Get the end position]
cl <- makeCluster(detectCores())
registerDoParallel(cl, cores = detectCores())
tmp_list_end <- list()
tmp_list_end <- foreach (i = 1:length(end_tmp), .packages = "binr") %dopar% {
  tmp <- bins(start_tmp[i]:end_tmp[i],target.bins = 100,max.breaks = 100, exact.groups = TRUE)
  e <- gsub(".*,", "", names(tmp$binct))
  e <- gsub("]","",e)
  tmp_list_end[[i]] <- as.numeric(e)
}
stopCluster(cl)

v1 <- unlist(tmp_list_start)
v2 <- unlist(tmp_list_end)
df_bin <- data.frame(start=v1, end=v2)


# Add back important information into the bined data (An ID that can identify the orginal region)
chr_bin <- rep(as.character(bed$V1), each=100)
ID_bin <- rep(as.character(bed$mergeID), each=100)

df_bin$seqnames <- chr_bin
df_bin$mergeID <- ID_bin

# Reorganise the format
df_bin <- select(df_bin, seqnames, start, end, mergeID)

# Save output to disk in bed format
write.table(df_bin, file.path(outDir, "bined.bed"), quote=F ,sep="\t", row.names=F, col.names=F)


# Command to use on cluster later on 

#  bedtools coverage -a bined.bed -b new_bg.bed -counts > bedtools_cov_output.txt