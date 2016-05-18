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
head(df)

# Melt data
m <- melt(df)
colnames(m) <- c("region", "sample", "value")
# Add column for colouring graph later
m$simple <- gsub("_r.*", "", m$sample)

cols <- c("#b30000", 
          "#d7301f", 
          "#084081", "#0868ac", "#2b8cbe",
          "#004529", "#006837", "#78c679",
          "#084081", "#0868ac", "#2b8cbe",
          "#004529", "#006837", "#78c679")

# Draw boxplot
pdf(file.path(outputRegionDir, "Region_distibution.pdf"), height=7, width=11)
ggplot(m, aes(x=region, y=value, fill=sample)) +
  geom_bar(stat="identity", position="dodge") +
  scale_y_continuous(labels = comma, name = "Number of reads") +
  scale_fill_manual(values=cols) +
  ggtitle("Number of reads assigned to each genomic region")
dev.
  
  
