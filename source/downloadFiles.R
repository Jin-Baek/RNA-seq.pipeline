##############################################
## RNA-seq pipeline [Download files]
##############################################
setRepositories(ind=1:8)
library(stringr)
library(foreach)
library(parallel)
library(doParallel)
library(iterators)

print("downloadFiles.R started")

# Define args 
args = commandArgs(trailingOnly=TRUE)
#args[1] <- reference genome download URL           <http://ftp.ensembl.org/pub/release-105/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz>
#args[2] <- GTF annotation download URL             <http://ftp.ensembl.org/pub/release-105/gtf/gallus_gallus/Gallus_gallus.GRCg6a.105.gtf.gz>

# Set parallel option
ci <- makeCluster(16)
registerDoParallel(ci)

# Set directory
dir=paste0(system("pwd",intern = TRUE),"/")
setwd(dir)

# Set sub directory 
reference_path <- paste0(dir,"reference/")
annotation_path <- paste0(dir,"annotation/")

# Regardless of the extension, download the file to the appropriate folder and pass if it's already exist
foreach(i=1:length(args),.packages = "stringr",.combine = c) %dopar% {
  frag_slash <- unlist(str_split(args[i],"/"))
  frag_dot <- unlist(str_split(gsub(".gz","",frag_slash[length(frag_slash)]),fixed('.')))
  # Check if the extension of args link is fa or gtf
  if(frag_dot[length(frag_dot)]=="fa"){
    # Pass if it is already downloaded
    if(frag_slash[length(frag_slash)] %in% list.files(reference_path) || gsub(".gz","",frag_slash[length(frag_slash)]) %in% list.files(reference_path)){
      print("reference file is already exist!!")
    }else{
      system(paste0("wget -P ",reference_path," ",args[i]))
      print("reference file is downloaded")
    }
  }else if(frag_dot[length(frag_dot)]=="gtf"){
    # Pass if it is already downloaded
    if(frag_slash[length(frag_slash)] %in% list.files(annotation_path) || gsub(".gz","",frag_slash[length(frag_slash)]) %in% list.files(annotation_path)){
      print("gtf file is already exist!!")
    }else{
      system(paste0("wget -P ",annotation_path," ",args[i]))
      print("gtf file is downloaded")
    }
  }
}
# *[wget http://ftp.ensembl.org/pub/release-105/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz]
# *[wget http://ftp.ensembl.org/pub/release-105/gtf/gallus_gallus/Gallus_gallus.GRCg6a.105.gtf.gz]

print("downloadFiles.R finished")