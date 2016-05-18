# 29/04/16
# Modify the format of bed graph to use with bedtools intersect on luster

workDir <- "perc_covered_fromBedGraph"
outDir <- file.path(workDir, "output")
dir.create(outDir)
# Read bedGraph file
bg <- read.delim(file.path(workDir, "data/example.bg"), header =F)

# Reformat the bedGraph file
bg_multi <- bg$V4 # This column contains how many reads are in a particular region
                  # Need to multiply with this number

# Multiply each bed position to thier corresponding number of reads
chr <- list()
start <- list()
end <- list()

# This is the step which takes a long time (the actual file has about 100 million lines)
for (i in seq_along(bg_multi)){
  chr[[i]] <- as.character(rep(bg$V1[i], each=bg_multi[i]))
  start[[i]] <- rep(bg$V2[i], each=bg_multi[i])
  end[[i]] <- rep(bg$V3[i], each=bg_multi[i])
}

chr_vec <- unlist(chr)
start_vec <- unlist(start)
end_vec <- unlist(end)

bg_new <- data.frame(chr = chr_vec,
                     start = start_vec,
                     end = end_vec)


write.table(bg_new, file.path(outDir,"new_bg.bed"), quote=F, sep="\t", row.names=F, col.names=F)

# ------------------------------------- # 
# Misc comments, please ignore
# ------------------------------------- # 

#write.table(bg_new, "new_bg.bed", quote=F, sep="\t", row.names=F, col.names=F)


# Command to use on cluster later on 


#  bedtools coverage -a example.bed -b new_bg.bed -counts > out.txt

#x <- capture.output(for (i in seq_along(bg_multi))
#  print(as.character(rep(bg$V1[i], each=bg_multi[i]))))


#x <- capture.output(for (i in 1:5)
#  print(as.character(rep(bg$V1[i], each=bg_multi[i]))))

