##############################################
## RNA-seq pipeline [Quantification]
## args: GTF annotation file, Bam file from alignment
##############################################

setRepositories(ind=1:8)
library(stringr)
library(dplyr)
library(foreach)
library(parallel)
library(doParallel)
library(iterators)

print("quantification.R started")

# Define args 
args = commandArgs(trailingOnly=TRUE)
#args[1] <- annotation file                    <Gallus_gallus.GRCg6a.105.gtf.gz> 
#args[2] <- Bam file from alignment             <Gallus_gallus_SRR42_sorted.bam>
#args[3] <- ...

# Set parallel option
ci <- makeCluster(16)
registerDoParallel(ci)

# Set directory
dir=paste0(system("pwd",intern = TRUE),"/")
setwd(dir)

# Set sub directory
bam_path <- paste0(dir,"bam/")
annotation_path <- paste0(dir,"annotation/")
quantification_path <- paste0(dir,"quan/")

# Perform quantification with |featureCounts|
foreach(i=seq(2,length(args)),.packages = "stringr",.combine = c) %dopar% {
  name <- gsub("sorted.bam","",args[i])
  comm <- "featureCounts -a "
  system(paste0(comm,annotation_path,args[1]," -o ",quantification_path,name,"quantification.txt ",bam_path,args[i]))
}
# *[featureCounts -a annotation_file -o output_file input_file]

print("quantification.R finished")
