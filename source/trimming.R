##############################################
## RNA-seq pipeline [Trimming]    
##############################################

setRepositories(ind=1:8)
library(stringr)
library(foreach)
library(parallel)
library(doParallel)
library(iterators)

print("trimming.R started")

# Define args 
args = commandArgs(trailingOnly=TRUE)
# args[1] <- trimmomatic path        </program/Trimmomatic/trimmomatic-0.39.jar>
# args[2] <- adapters path           </program/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2:True>
# args[3] <- End mode                                                 <SE or PE>
# args[4] <- trimmomatic parameter : LEADING                                 <3>
# args[5] <- trimmomatic parameter : TRAILING                                <3>
# args[6] <- trimmomatic parameter : MINLEN                                 <36>
# Input according to sample End Mode 
# args[7] <- sample fastq.gz file                             <SRR42_1.fastq.gz>
# args[8] <- sample fastq.gz file                             <SRR42_2.fastq.gz>
# args[9] <- ... 

# Set parallel option
ci <- makeCluster(16)
registerDoParallel(ci)

# Set directory
dir=paste0(system("pwd",intern = TRUE),"/")
setwd(dir)

# Set sub directory
samples_path <- paste0(dir,"samples/") 
trim_path <- paste0(dir,"trim/")
trimmomatic_path <- args[1] 
adapters_path <- args[2]

# Perform trimming on all samples with |Trimmomatic|
if(args[3]=="SE"){
  foreach(i=seq(7,length(args)),.packages = "stringr",.combine = c) %dopar% {
    sampleCode <- paste0(unlist(str_split(args[i],"_"))[1],"_")
    if(length(grep(sampleCode,list.files(path = trim_path)))==0){
      trimming_cmd <- paste0("java -jar ",trimmomatic_path," SE -trimlog ",trim_path,sampleCode,"log.txt ",samples_path,args[i]," ",
                             trim_path,sampleCode,"output.fastq.gz ","ILLUMINACLIP:",adapters_path,":2:30:10:2:True LEADING:"
                             ,args[4]," TRAILING:",args[5]," MINLEN:",args[6])
      system(trimming_cmd)
    }else{
      print(paste0(sampleCode," trimmed samples is already exist!!"))
    }
  }
}else if(args[3]=="PE"){
  foreach(i=seq(7,length(args),2),.packages = "stringr",.combine = c) %dopar% {
    sampleCode <- paste0(unlist(str_split(args[i],"_"))[1],"_")
    if(length(grep(sampleCode,list.files(path = trim_path)))==0){
      trimming_cmd <- paste0("java -jar ",trimmomatic_path," PE -trimlog ",trim_path,sampleCode,"log.txt ",samples_path,args[i]," ",samples_path,
                             args[i+1]," ",trim_path,sampleCode,"output_forward_paired.fastq.gz ",trim_path,sampleCode,"output_forward_unpaired.fastq.gz "
                             ,trim_path,sampleCode,"output_reverse_paired.fastq.gz ",trim_path,sampleCode,"output_reverse_unpaired.fastq.gz ",
                             "ILLUMINACLIP:",adapters_path," LEADING:",args[4]," TRAILING:",args[5]," MINLEN:",args[6])
      system(trimming_cmd)
    }else{
      print(paste0(sampleCode," trimmed samples is already exist!!"))
    }
  }
}else{
  print("Need to specify correct [End Mode]")
}
# *[java -jar /program/Trimmomatic/trimmomatic-0.39.jar PE -trimlog trimlog.txt SRR42_1.fastq.gz SRR42_2.fastq.gz output_forward_paired.fastq.gz output_forward_unpaired.fastq.gz output_reverse_paired.fastq.gz output_reverse_unpaired.fastq.gz ILLUMINACLIP:/program/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36]

print("trimming.R finished")