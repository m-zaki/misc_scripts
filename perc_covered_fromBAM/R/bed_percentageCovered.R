# Given a bed file, we would like to know the percentage of bases covered by more than x sequencing reads

workDir <- "perc_covered_fromBAM"
bedDir <- file.path(workDir, "data")
bamDir <- file.path(workDir, "data")
outDir <- file.path(workDir, "output")
dir.create(outDir)

source(file.path(workDir, "R/function.R"))


library(GenomicRanges)
# Packages to split bed and bin data
library(binr)
library("parallel")
library("foreach")
library("doParallel")
# Packages to read BAM file into R
library(GenomicAlignments)
# Data wrangling
library(dplyr)
library(tidyr)
library(reshape)
library(lazyeval)
# Packages for read counting
library(Rsubread)

# List all the files in the folder
files <- list.files(bedDir, pattern="*.bed",full.names = TRUE) # List all .bed files
names(files) <- sapply(files, function(fn) { gsub('.bed', '', basename(fn)) })


for (gg in seq_along(files)){

# Read bed file into R
bed <- read.delim(files[gg], header=F)
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

# Transform bineed object into G ranges

x <- df_bin
x$rep <- rep(1:100, (nrow(x) / 100) )

G <- GRanges(
  seqnames = x$seqnames,
  ranges = IRanges(start=x$start,
                   end=x$end),
  id = paste(x$mergeID, x$rep, sep="_r"),
  mergeID = x$mergeID
)

#-----------------------------------#
# Read the bam file into R
#-----------------------------------#
sbp <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=F, 
                                     isNotPrimaryRead=F, 
                                     isNotPassingQualityControls=F, 
                                     isDuplicate=NA))

bam_list <- list.files(file.path(bamDir),
                       recursive=FALSE, include.dirs=FALSE,
                       pattern="*.bam$", full=TRUE)


#-----------------------------------#
# Counting using feature Counts
#-----------------------------------#
# Need to create AnnotationFile to use with feature counts
Ganno <- createAnnotationFile(G)

res <- featureCounts(files=bam_list,
                         annot.ext=Ganno,
                         isGTFAnnotationFile=FALSE,
                         allowMultiOverlap=TRUE, 
                         nthreads=8,
                         strandSpecific=0, ### Adjust strand accordingly
                         countMultiMappingReads=FALSE,
                         isPairedEnd=FALSE,
                         requireBothEndsMapped=FALSE,
                         checkFragLength=FALSE)

dat <- as.data.frame(res$counts)
colnames(dat) <- gsub("perc_covered_fromBAM.data.", "", colnames(dat))
col_tmp <- colnames(dat)
names(col_tmp) <- colnames(dat)
dat$mergeID <- gsub("_r.*", "", row.names(dat)) # Remove everything after "_r"



perc_sample_list <- list()



for (r in 1:(length(seq_along(dat)) - 1)) {
c1 <- col_tmp[r]
d1 <- select_(dat, "mergeID", interp(~c1, 
                                  c1=as.name(c1)))
colnames(d1) <- c("mergeID", "cov")
# Split data according to unique identifier
d_l <- split( d1 , f = d1$mergeID )


# Create some temporary objects
perc_cov_list <- list()
gene_name <- vector() # Store the gene_name
val_cov  <- vector() # Temporatly store the percentage of genes covered of individual gene

# Use the function "per.covered"
for (i in seq_along(d_l)) {
  perc_cov_list[[i]] <- perc.covered(input = d_l, min.depth = 5)
  gene_name[i] <- names(perc_cov_list[[i]])
  val_cov[i] <- as.vector(perc_cov_list[[i]])
}

perc_cov_df <- data.frame(mergeID = gene_name,
                          coverage_perc = val_cov)
colnames(perc_cov_df)[2] <- paste0("perc_", names(col_tmp)[r]) # Rename the columns so each sample has a different name
perc_sample_list[[r]] <- perc_cov_df
}

perc_df <- Reduce(function(x,y) {merge(x,y)}, perc_sample_list)

write.table(perc_df, file.path(outDir, paste0(names(files)[gg], "_percentage.txt")),
            sep="\t", col.names=NA, quote=F)
}
