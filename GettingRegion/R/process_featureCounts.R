# Generate plot to visualise number of reads assigned to each region

library(dplyr)
library(ggplot2)
library(reshape)
library(scales) 

# Using the same Directory structure
#workDir <- "/Volumes/StemCellBiol1/Personal Folders/zaki/project/misc/misc_scripts/GettingRegion/"
workDir <- "GettingRegion"
dataDir <- file.path(workDir, "data")
outputDir <- file.path(workDir, "output")
outputRegionDir <- file.path(outputDir, "bedRegion")


# load featureCount results
load(file.path(outputRegionDir, "Rsubread_raw_output.Rdata"))
# Number of reads assigned are stored in the "res" object
df <- res$counts
colnames(df) <- gsub("X.scratch.zaki.project.lowInput.20160510_comparison.bam.sub1.", "", colnames(df))
colnames(df) <- gsub(".bam", "", colnames(df))     

# Melt data
m <- melt(df)
colnames(m) <- c("region", "sample", "value")


# Draw bar plot
pdf(file.path(outputRegionDir, "Region_distibution.pdf"), height=7, width=11)
ggplot(m, aes(x=sample, y=value, fill=region)) +
  geom_bar(stat="identity", position="fill")  +
  scale_y_continuous(labels = percent_format(), name = "Percentage of reads") +
  ggtitle("Number of reads assigned to each genomic region") +
  theme(axis.text.x = element_text(angle = 45, hjust= 1))
dev.off()






  
