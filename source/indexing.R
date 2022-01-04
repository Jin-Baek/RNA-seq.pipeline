##############################################
## RNA-seq pipeline [Indexing]    
##############################################

setRepositories(ind=1:8)
library(stringr)
library(foreach)
library(parallel)
library(doParallel)
library(iterators)

print("indexing.R started")

# Define args 
args = commandArgs(trailingOnly=TRUE)
# args[1] <- reference genome file     <Gallus_gallus.GRCg6a.dna.toplevel.fa.gz>
# args[2] <- ...

# Set parallel option
ci <- makeCluster(16)
registerDoParallel(ci)

# Set directory
dir=paste0(system("pwd",intern = TRUE),"/")
setwd(dir)

# Set sub directory
reference_path <- paste0(dir,"reference/")
index_path <- paste0(dir,"index/")

# Perform indexing for assigned reference genome with |Hisat2|
foreach(i=1:length(args),.combine = c,.packages = "stringr") %dopar% {
  if(str_detect(args[i],".gz")){
    # Unzip reference genome
    system(paste0("gzip -d ",reference_path,args[i]))
    print("Unzip is performed!!")
  }else{
    print("Unzip is already performed!!")
  }
  comm <- "hisat2-build -p 8 "
  reference.type <- unlist(str_split(args[i],fixed('.')))[1]
  if(length(grep(paste0(reference.type,"_IndexFile"),list.files(path=index_path)))!=0){
    print("index file is already exist!!!")
  }else{
    comm <- paste0(comm,reference_path,gsub(".gz","",args[i])," ",index_path,reference.type,"_IndexFile")
    system(comm)
  }
}
# *[gzip -d Gallus_gallus.GRCg6a.dna.toplevel.fa.gz]
# *[hisat2-build -p 8 filename.fa indexfileName]

print("indexing.R finished")