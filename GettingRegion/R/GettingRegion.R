# Determine the proportion of reads mapped to various region
# Eg, exonic, intronic, intergenic


# Generate region that contains the following ;
# 1) Exon region
# 2) Intergenic region 
# 3) Intronic region

library(GenomicAlignments)
library(EnsDb.Hsapiens.v75)
library(GenomicRanges)

#workDir <- "/Volumes/StemCellBiol1/Personal Folders/zaki/project/misc/misc_scripts/GettingRegion/"
workDir <- "GettingRegion"
dataDir <- file.path(workDir, "data")
outputDir <- file.path(workDir, "output")
dir.create(outputDir)
outputRegionDir <- file.path(outputDir, "bedRegion")
dir.create(outputRegionDir)

# Generate region
# 1) Getting exonic region

# Extract the exonic region for all genes
er <- exons(EnsDb.Hsapiens.v75, 
             columns=c("tx_biotype","gene_id", "gene_name", "gene_biotype"),
             filter=list(SeqnameFilter(c(1:22, "X", "Y", "MT"))))

# Remove LRG_gene 
er_sub <- er[!er$gene_biotype == "LRG_gene",]
# Collapse overlaping region
er_reduce <- reduce(er_sub)
exon_df <- as.data.frame(er_reduce)
# Now we have all the exonic region

# 2) Getting intronic region
# Begin by getting the whole region of a particular gene
# This information is available before hand
geneRegion <- read.delim(file.path(dataDir, "fullGeneregion.txt"))
# Convert to Granges
Gr <- GRanges(seqnames = geneRegion$seqnames,
              IRanges(start = geneRegion$start,
                      end = geneRegion$end),
                      strand = geneRegion$strand)
# To obtain intronic region, find the difference between 1) Gene region and 2) Exonic region
# https://support.bioconductor.org/p/69614/
intron <- setdiff(Gr, er_reduce)
intron_df <- as.data.frame(intron)

# 3) Getting intergenic region
# Begin by obtaining the size of a whole chromosome
# Intonic region just be a substraction of the whole crhomosome to Gene region

# Size of hg19 chromosome can be found here - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
chrSize <- read.delim(file.path(dataDir, "hg19_ChrSize.txt"))
GrSize <- GRanges(seqnames = chrSize$seqnames,
              IRanges(start = chrSize$start,
                      end = chrSize$end))
# To obtain intergeinic region, find the difference between 1) Whole Chromosome and 2) Gene Region
intergenic <- setdiff(GrSize, Gr, ignore.strand=TRUE)
intergenic_df <- as.data.frame(intergenic)

# Combine everything into one bed file to be used as input for featureCounts

# Add appropriate identifier
exon_df$GeneID <- "exon"
intron_df$GeneID <- "intron"
intergenic_df$GeneID <- "intergenic"


allRegion <- rbind(exon_df,
                      intron_df,
                      intergenic_df)


# Rearangge the file into SAF format
# http://bioinf.wehi.edu.au/featureCounts/
allRegion_df <- data.frame(GeneID = allRegion$GeneID,
                           Chr = allRegion$seqnames,
                           Start = allRegion$start,
                           End = allRegion$end,
                           Strand = allRegion$strand)

write.table(allRegion_df, file.path(outputRegionDir, "allRegion.saf"), 
            quote=F, sep="\t", row.names = F)





