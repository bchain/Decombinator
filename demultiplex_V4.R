#If not already set, set path which ontains the index files; 
#and will contain the output folder. 
#setwd("31_12_2014\\")
#set output directory within this working directory
output<-"output"
dir.create(output)
library(parallel)
library("ShortRead")
library("stringdist")

cores<-detectCores()
cl <- makeCluster(cores) 
#######################################################################################################################
#number of missmatches allowed 
mm1<-1
mm2<-2

#expt number
expt<-"38"

#position of barcode in read 3
start_bc<-1
end_bc<-6

#position of index in read 1 (at the moment usually positions 7-12 immediately after barcode)
indexR1_b<-7
indexR1_e<-12
#path to data
path_data<-paste("c:\\New folder\\expt_",expt,"\\",sep="")

#data files names
name1<-"R1.fastq.gz"
name2<-"R2.fastq.gz"
name3<-"R3.fastq.gz"
########################################################################################################################
#read in indeces for SP1
SP1_i<-read.table("SP1_indeces.txt",sep="\t",stringsAsFactors=FALSE)
SP2_i<-read.table("SP2_indeces.txt",sep="\t",stringsAsFactors=FALSE)
#calcualte hamming distnace between indeces.
#a<-stringdistmatrix(SP1_i[2,],SP1_i[2,],method="hamming")
#note that 5 and 14 only differ by two
#read in indeces for SP2


#function which does fuzzy match to SP1 and SP2
fmatch_SP1<-function(x){agrep(x,SP1_i[2,],max.distance=list(ins=0,del=0,sub=mm1))}
fmatch_SP2<-function(x){agrep(x,SP2_i[2,],max.distance=list(ins=0,del=0,sub=mm2))}

#fields with for data
file1<-paste(path_data,name1,sep="")
file2<-paste(path_data,name2,sep="")
file3<-paste(path_data,name3,sep="")

#open file connections
seqR1<-FastqStreamer(file1,n=1E6)
seqR2<-FastqStreamer(file2,n=1E6)
seqR3<-FastqStreamer(file3,n=1E6)


#iterate over file
ptm <- proc.time()
while (length(seqsR1_fq<-yield(seqR1))){
  #seqsR1_fq<-yield(seqR1)
#extract barcode sequence  from start of read 3 and paste to end of read 1, and create a new fastQ object
seqsR3_fq<-yield(seqR3)
if(length(seqsR3_fq)==length(seqsR1_fq)){bar_codes<-substr(sread(seqsR3_fq),start_bc,end_bc)} else {bar_codes<-substr(sread(seqsR3_fq[1:length(seqsR1_fq)]),start_bc,end_bc)}         
#seqsR1_bc<-DNAStringSet(paste(sread(seqsR1_fq),bar_codes,sep=""))
R1_reads<-sread(seqsR1_fq)
seqsR1_bc<-DNAStringSet(paste(R1_reads,bar_codes,sep=""))
#adjust quality reads to same length 
q<-BStringSet(paste(quality(quality(seqsR1_fq)),"XXXXXX",sep=""))
#create a new FastQ file for output
seqsR1_fq_bc<-ShortReadQ(sread=seqsR1_bc,quality=q,id=id(seqsR1_fq))

#a<-BStringSet(quality(seqsR1_fq))

#extract index from read 1 (usually found between position 7 and 12; set above)
#index_seq<-substr(R1_reads,7,12)
index_seq<-substr(R1_reads,indexR1_b,indexR1_e)
#index using agrep
#first move the necessary data files and objects to the parallel cores
clusterExport(cl,c("fmatch_SP1","mm1","SP1_i"))
#then run parallel version of Lapply
index_1<-parLapply(cl,index_seq,fmatch_SP1)
#index_1<-parSapply(cl,index_seq,fmatch_SP1)

#check for ambiguous elements
l_list<-sapply(index_1,length)
#all ambiguous indices will be indexed 99
index_1[which(l_list>1)]<-99
index_1[which(is.na(index_1>0))]<-0


#Repeat for SP2 index
#extract index sequences 
index_seq2<-substr(sread(yield(seqR2)),1,6)

#index using agrep
clusterExport(cl,c("fmatch_SP2","mm2","SP2_i"))
index_2<-parLapply(cl,index_seq2,fmatch_SP2)


#check for ambiguous indices
l_list<-sapply(index_2,length)
index_2[which(l_list>1)]<-99
index_2[which(is.na(index_2>0))]<-0

#combine indeces
comb_seq_index<-paste(unlist(index_1,use.names=FALSE), unlist(index_2,use.names=FALSE),sep="_")
#factor all possible indeces
index_com<-levels(factor(comb_seq_index))

#save each index set separately as fast_q file with index name, and expt name 

for(fn in index_com){
  writeFastq(seqsR1_fq_bc[which(comb_seq_index==fn)],file=paste(output,"/",fn,"_",expt,".fastq.gz",sep=""),full=TRUE,mode="a",compress=TRUE)
          }
    }   #end of while loop cycling through sequences
close(seqR1)
close(seqR2)
close(seqR3)
#stop the clock
proc.time() - ptm
#close clusters
stopCluster(cl)

########################################################################################################################
#count number of seqeunces in each file and output a summary Table
dir_output<-dir("output\\")
seq_names<-strsplit(dir_output,"_")
output_summary<-rep(0,length(dir_output)*3)
dim(output_summary)<-c(length(dir_output),3)
output_summary[,1]<-sapply(seq_names,"[[", 1,simplify=TRUE )
output_summary[,2]<-sapply(seq_names,"[[", 2, simplify =TRUE )
#cycle thyough files and count the reads
f<-1
for (file in dir_output){
  file_n<-paste("output\\",file,sep="")
  
Seq<-FastqStreamer(file_n)
l<-0
while (length(count<-yield(Seq))){
  l<-length(count)+l
    }#end of fastq_file
close(Seq)
output_summary[f,3]<-l
f<-f+1
      }#end of directory
write.table(output_summary,"output_summary.txt",sep="\t",row.names=FALSE,col.names=FALSE)
#plot all columns whihc contain more than 0.1% reads
plot_sum<-sort(as.numeric(output_summary[,3]),decreasing=TRUE)
p<-plot_sum/sum(plot_sum)
plot_sum_f<-plot_sum[which(p>0.001)]
x<-(1:length(plot_sum_f))
barplot(plot_sum_f,names.arg=x,cex.names=0.8,cex.axis=0.8)
imageSave(paste("frequency_plot_",expt,".png",sep=""))


#rm(list=ls(all=TRUE))
