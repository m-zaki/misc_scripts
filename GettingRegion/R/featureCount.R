# Use featureCounts to determine mapping to each region

library(Rsamtools)
library(GenomicRanges)
library(Rsubread)

# Using the same Directory structure
#workDir <- "/Volumes/StemCellBiol1/Personal Folders/zaki/project/misc/misc_scripts/GettingRegion/"
workDir <- "GettingRegion"
dataDir <- file.path(workDir, "data")
outputDir <- file.path(workDir, "output")
outputRegionDir <- file.path(outputDir, "bedRegion")

# Location of bam files
bam.dir <- dataDir
bam.files <- list.files(bam.dir, pattern=".bam$", full.names=TRUE)


# Read the custom annotation file prepared from previous step
exampleSAF <- read.delim(file.path(outputDir, "allRegion.saf"))

res <- featureCounts(files=bam.files,
                     annot.ext=exampleSAF,
                     isGTFAnnotationFile=FALSE,
                     useMetaFeatures=TRUE,
                     nthreads=8,
                     strandSpecific=0,
                     countMultiMappingReads=FALSE,
                     isPairedEnd=FALSE,
                     requireBothEndsMapped=FALSE,
                     checkFragLength=FALSE)

save(res, file=file.path(outputRegionDir, "Rsubread_raw_output.Rdata"))
