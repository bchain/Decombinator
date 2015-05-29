#this function calculates the lengths of three files called R1.fastq.gz, R2.fastq.gz
#R3.fastq.gz as a vector l of length 3. 
fastq_length(path="pathfile"){

l<-c()
file1<-paste(pathfile,R1.fastq.gz,sep="")
file2<-paste(pathfile,R2.fastq.gz,sep="")
file3<-paste(pathfile,R3.fastq.gz,sep="")
###########################################
l1<-0
seqR1<-FastqStreamer(file1,n=1E6)
while (length(seqsR1_fq<-yield(seqR1))){l1<-l1+length(seqsR1_fq)}
l[1]<-l1
close(seqR1)
############################
l2<-0
n<-1
seqR2<-FastqStreamer(file2,n=1E6)
while (length(seqsR2_fq<-yield(seqR2))){l2<-l2+length(seqsR2_fq)
n<-n+1}
close(seqR2)
l[2]<-l2
#################################
l3<-0
n<-1
seqR2<-FastqStreamer(file3,n=1E6)
while (length(seqsR3_fq<-yield(seqR3))){l3<-l3+length(seqsR3_fq)}
close(seqR3)
l[3]<-l3
return(l)
        }

##############################
seqsR1_fq<-yield(seqR1)
seqsR3_fq<-yield(seqR3)
writeFastq(head(seqsR1_fq),"R1_fq.fastq",compress=FALSE)
writeFastq(head(seqsR3_fq),"R3_fq.fastq",compress=FALSE)
head(seqsR3_fq)
seqsR1_fq[2][[1]]
writeFastq(tail(seqsR1_fq),"R1_fq_tail.fastq",compress=FALSE)
writeFastq(tail(seqsR3_fq),"R3_fq_tail.fastq",compress=FALSE)

writeFastq(seqsR3_fq[17539550:17539560],"R3_fq_double.fastq",compress=FALSE)

