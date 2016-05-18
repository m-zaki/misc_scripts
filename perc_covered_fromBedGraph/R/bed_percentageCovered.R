# Process the bedtools coverage output

library(dplyr)


source("perc_covered_fromBedGraph/R/function.R")

workDir <- "perc_covered_fromBedGraph"

dat <- read.delim(file.path(workDir, "output/bedtools_cov_output.txt"), header= F)

# Divide data into list based on individual ID
# This step takes the longest on a big data (anyway of doing this with dplyr??)
dat_l <- split( dat , f = dat$V4 ) #column V4 is the peak unique ID

# Create some temporary objects
perc_cov_list <- list()
gene_name <- vector() # Store the gene_name
val_cov  <- vector() # Temporatly store the percentage of genes covered of individual gene

# Use the function "per.covered"
for (i in seq_along(dat_l)) {
  perc_cov_list[[i]] <- perc.covered(input = dat_l, min.depth = 1)
  gene_name[i] <- names(perc_cov_list[[i]])
  val_cov[i] <- as.vector(perc_cov_list[[i]])
}

perc_cov_df <- data.frame(gene_name = gene_name,
                          coverage_perc = val_cov)
# colnames(perc_cov_df)[2] <- paste0("perc_", names(files)[r]) # Rename the columns so each sample has a different name
# perc_sample_list[[r]] <- perc_cov_df
