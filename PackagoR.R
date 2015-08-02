#script to combine the illumina demultiplexed script into one,
#whihc can then be demultiplexed again using the DemultiplexoR.
#The output names are R1 (forward read); R2 (Illumina index read);
#R3 Reverse read. 

library("ShortRead")
#input adn output paths

input<-"C:/TCR/input/"
output<-"C:/TCR/output/"

#select different Illumina reads
R1<-grep("_R1_",as.character(dir(input)))
R2<-grep("_R2_",as.character(dir(input)))
I1<-grep("_I1_",as.character(dir(input)))
#loop though files and combine all the I1, R1 and R2 reads together; 
# Rename and save them as R2, R1 and R3 respectively

i<-1
for (i in 1:length(I1)){
index<-readFastq(paste(input,dir(input)[I1[i]],sep=""))
writeFastq(index,file=paste(output,"R2.fastq.gz",sep=""),full=TRUE,mode="a",compress=TRUE)
print(dir(input)[I1[i]])
}

for(i in 1:length(R1)){
  file_R1<-readFastq(paste(input,dir(input)[R1[i]],sep=""))
  writeFastq(file_R1,file=paste(output,"R1.fastq.gz",sep=""),full=TRUE,mode="a",compress=TRUE)
  print(dir(input)[R1[i]])
}

for(i in 1:length(R2)){
  file_R2<-readFastq(paste(input,dir(input)[R2[i]],sep=""))
  writeFastq(file_R2,file=paste(output,"R3.fastq.gz",sep=""),full=TRUE,mode="a",compress=TRUE)
  print(dir(input)[R2[i]])
}

