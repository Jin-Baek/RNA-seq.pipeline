##############################################
## RNA-seq pipeline [Main]
##############################################

setRepositories(ind=1:8)
library(stringr)
library(dplyr)
library(foreach)
library(parallel)
library(doParallel)
library(iterators)

# Define args (args is variable factor)
args = commandArgs(trailingOnly=TRUE)

# Set base directory (all files should be on base path)
dir=paste0(system("pwd",intern = TRUE))
setwd(dir)
samples_path <- paste0(dir,"/samples/") 
trim_path <- paste0(dir,"/trim/")
bam_path <- paste0(dir,"/bam/")

# Make directories 
DirNeed <- c("reference","samples","annotation","trim","index","sam","bam","quan")

for(dirname in paste0(dir,"/",DirNeed)){
  x<- dirname %in% list.dirs(dir)
  if(x){
    print(paste0(dirname," is already exist!"))
  }else{
    system(paste0("mkdir ",dirname))
  }
}


#======= downloadFiles.R 
vf <- unlist(str_split(args," "))

dind <- grep("-df",vf)+1  
url.list <- c()
bol <- TRUE
while(bol){
  if(str_detect(vf[dind],pattern = "-[a-z]")||is.na(vf[dind])){
    bol <- FALSE
  }else{
    url.list <- append(url.list,vf[dind])
    dind <- dind+1
  }
}

system(paste0("Rscript downloadFiles.R ",paste(url.list,collapse = " ")))


#======= download samples
setwd(samples_path)
textname <- vf[grep("-sm",vf)+1]
ndsample <- c()
for(i in readLines(paste0(samples_path,textname))){
  if(length(grep(i,list.files(path = samples_path,pattern = "*.fastq.gz")))!=0){
    print(paste0(i," is already exist!!"))
  }else{
    print(paste0(i," should be download!!"))
    ndsample <- append(ndsample,i)
  }
}

# Make new sample txt file(delete the sample which are already downloaded) overload on original sample txt file
write.table(ndsample,paste0(samples_path,textname),row.names = FALSE,quote = FALSE,col.names = FALSE)
system(paste0(vf[grep("-dump",vf)+1]," ",samples_path,textname))
setwd(dir)


#======= trimming.R
trimmomatic_path <- vf[grep("-Ptrim",vf)+1]
adapters_path <- vf[grep("-Padapt",vf)+1]
endMode <- vf[grep("-em",vf)+1]

# trimmomatic parameter
if(grepl("-tp1",args)){  
  tp_LEADING <- as.integer(vf[grep("-tp1",vf)+1])
}else{
  tp_LEADING <- 3
}

if(grepl("-tp2",args)){  
  tp_TRAILING <- as.integer(vf[grep("-tp2",vf)+1])
}else{
  tp_TRAILING <- 3
}

if(grepl("-tp3",args)){  
  tp_MINLEN <- as.integer(vf[grep("-tp3",vf)+1])
}else{
  tp_MINLEN <- 36
}

# Samples
sample.list <- list.files(path=samples_path,pattern = "*.fastq.gz")

system(paste0("Rscript trimming.R ",trimmomatic_path," ",adapters_path," ",endMode," ",tp_LEADING," ",tp_TRAILING," ",tp_MINLEN," ",
              paste(sample.list,collapse = " ")))



#======= indexing.R 
reference <- vf[grep("-ref",vf)+1]

system(paste0("Rscript indexing.R ",reference))



#======= alignment.R
species <- vf[grep("-species",vf)+1]
indexFile <- vf[grep("-i",vf)+1]

# trimmed samples
if(endMode=="SE"){
  tsample.list <- list.files(path=trim_path,pattern = "_output.fastq.gz")
}else if(endMode=="PE"){
  tsample.list <- list.files(path=trim_path,pattern = "_paired.fastq.gz")
}

system(paste0("Rscript alignment.R ",indexFile," ",endMode," ",paste(tsample.list,collapse = " ")))



#======= quantification.R
annotation <- vf[grep("-an",vf)+1]

# bam files which mapped by each samples
bam.list <- list.files(path = bam_path,pattern = paste0("^",species,".*_sorted.bam"))

system(paste0("Rscript quantification.R ",annotation," ",paste(bam.list,collapse = " ")))

print("main.R finished")

