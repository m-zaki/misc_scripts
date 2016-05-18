# Function to extract the percentage of gene covered
perc.covered <- function(input, min.depth) {
  dat_l <- input[[i]]
  gene_id <- as.character(dat_l$mergeID) # column mergeID has the unique identifier
  cov_value <- dat_l$cov # column which contains the coverage
  
  # Getting the percentage of gene covered more than 1 base
  perc_cov <- round(sum(cov_value >= min.depth) / length(cov_value) * 100, 2)
  names(perc_cov) <- unique(gene_id)
  
  return(perc_cov)
}
