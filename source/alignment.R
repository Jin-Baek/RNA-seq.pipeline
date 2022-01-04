##############################################
## RNA-seq pipeline [Alignment], [sam to sorted bam]
##############################################

setRepositories(ind=1:8)
library(stringr)
library(dplyr)
library(foreach)
library(parallel)
library(doParallel)
library(iterators)

print("alignment.R started")

# Define args 
args = commandArgs(trailingOnly=TRUE)
# args[1] <- Index reference genome file               <Gallus_gallus_IndexFile>
# args[2] <- End mode                                                 <SE or PE>
# Input according to End Mode 
# args[3] <- Trimmed sample fastq.gz file   <SRR42_output_forward_paired.fastq.gz>
# args[4] <- Trimmed sample fastq.gz file   <SRR42_output_reverse_paired.fastq.gz>
# args[5] <- ...

# Set parallel option
ci <- makeCluster(16)
registerDoParallel(ci)

# Set directory
dir=paste0(system("pwd",intern = TRUE),"/")
setwd(dir)

# Set sub directory
index_path <- paste0(dir,"index/")
sam_path <- paste0(dir,"sam/")
bam_path <- paste0(dir,"bam/")
trim_path <- paste0(dir,"trim/")

# Perform alignment with |Hisat2|
name <- gsub("_IndexFile","",args[1])
if(args[2]=="SE"){
  foreach(i=seq(3,length(args)),.packages = "stringr",.combine = c) %dopar% {
    sample <- unlist(str_split(args[i],"_"))[1]
    if(length(grep(paste0(name,"_",sample,"_"),list.files(path = bam_path)))==0){
      comm <-"hisat2 -x "
      system(paste0(comm,index_path,args[1]," -U ",trim_path,args[i]," -S ",sam_path,name,"_",sample,".sam"))
      
      # Perform sam to bam & bam to sorted bam with |samtools|
      system(paste0("samtools view -bS ",sam_path,name,"_",sample,".sam > ",bam_path,name,"_",sample,".bam"))
      system(paste0("samtools sort ",bam_path,name,"_",sample,".bam -o ",bam_path,name,"_",sample,"_sorted.bam"))
    }else{
      print(paste0(name,"_",sample,"_"," samples bam file is already exist!!"))
    }
  }
}else if(args[2]=="PE"){
  foreach(i=seq(3,length(args),2),.packages = "stringr",.combine = c) %dopar% {
    sample <- unlist(str_split(args[i],"_"))[1]
    if(length(grep(paste0(name,"_",sample,"_"),list.files(path = bam_path)))==0){
      comm <-"hisat2 -x "
      system(paste0(comm,index_path,args[1]," -1 ",trim_path,args[i]," -2 ",trim_path,args[i+1]," -S ",sam_path,name,"_",sample,".sam"))
      
      system(paste0("samtools view -bS ",sam_path,name,"_",sample,".sam > ",bam_path,name,"_",sample,".bam"))
      system(paste0("samtools sort ",bam_path,name,"_",sample,".bam -o ",bam_path,name,"_",sample,"_sorted.bam"))
    }else{
      print(paste0(name,"_",sample,"_"," samples bam file is already exist!!"))
    }
  }
}else{
  print("Need to specify correct [End Mode]")
}
# *[hisat2 -x IndexFile -1 SRR1758114_1.fastq.gz -2 SRR1758114_2.fastq.gz -S chicken.sam 2> log_file &]
# *[samtools view -bS eg2.sam > eg2.bam]
# *[samtools sort eg2.bam -o eg2.sorted.bam]
print("alignment.R finished")
