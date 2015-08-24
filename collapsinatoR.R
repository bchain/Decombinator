#script to collapse seqeunce data to single barcodes
#and collect numbers of reads, TCRs, unique TCRs

#set expt name
expt<-47

#Define barcode sequence
#I8_6N_I8_6N
index<-"GTCGTGAT......GTCGTGAT......"

#Data input 
expt_path<-"C:/TCR/input/"
#Data output
output_path<-"C:/TCR/output/"

#from my work computer
#dropbox<-"D:/dropbox/"
#at home
#dropbox<-"C:/Users/Benny Chain/dropbox/"

####################################################################################################################
####################################################################################################################
#cycle through directory analysing each file in turn and 
#extracting unique TCR RNAs (i.e. collaping on bar codes) and 
#then unique TCRs and storing the 6 part identifier (including CDR3)
#in list file.
i<-1
output_sum<-c()
output_seq<-list()
#path to expt 
directory<-dir(expt_path)
output_sum<-matrix(rep(NA,7*length(directory)),ncol=7)

for( i in 1:length(directory)){
sample_namef<-directory[i]

#read in data 
data_path<-paste(expt_path,sample_namef,sep="")
TCRs<-read.table(data_path,stringsAsFactors=FALSE,header=FALSE,sep=",")
sample_name<-substr(directory[i],1,(nchar(directory[i])-4))
print(sample_name)
colnames(TCRs)<-c("V","J","Vdel","Jdel","insert","seqid","barcode","X")

#number of reads
n<-dim(TCRs)[1]
n

#look for correct barcode pattern when suitable index present
#use fuzzy matching to allow two missmatches
isf<-0
isf<-agrep(index,TCRs[,7],fixed=FALSE,max.distance=list(ins=0,del=0,sub=2))
n_u<-length(isf)
length(isf)/n
#select only those  sequences with barcodes
if(length(isf)<100){cat ("Very small number of sequences in sample ", sample_name)}
if (length(isf) != 0) {TCRs_bc<-TCRs[isf,]} else {TCRs_bc<-TCRs}

#look at distribution of barcode freqeuncies
bc_f<-tabulate(factor(TCRs_bc$barcode))
#mean and variance of barcode copy number
mean(bc_f)
var(bc_f)
max(bc_f)
#plot and save histogram of barcode family size
#hist(bc_f,breaks=c(0,1,10,100,max(bc_f)))$counts
#output<-paste(output_path,"barcode_f",sample_name,".png",sep="")
#imageSave(output)

#check for unique TCR RNA sequences (i.e. unique combination of barcode and TCR)
TCRs_M<-paste(TCRs_bc[,1],TCRs_bc[,2],TCRs_bc[,3],TCRs_bc[,4],TCRs_bc[,5],TCRs_bc[,7],sep=",")
TCRS_U<-levels(factor(TCRs_M))

#number of TCR RNA molecules (i.e. collapse multiple TCRs with same barcode)
RNA_u<-length(TCRS_U)
#split up TCRS_U into components again
RNA_U<-as.data.frame(matrix(unlist(strsplit(TCRS_U,split=",")),byrow=TRUE,ncol=6))
colnames(RNA_U)<-c("V","J","Vdel","Jdel","insert","bar_code")

########################################################
#collect sequences for later analysis in list called e.g. output_seq
output_seq[[i]]<-list(RNA_U)
names(output_seq)[i]<-sample_name


########################################################
#look for unique T sequences
RNA_U_paste<-paste(RNA_U[,1],RNA_U[,2],RNA_U[,3],RNA_U[,4],RNA_U[,5],sep=",")
TCRs_unique<-levels(factor(RNA_U_paste))
#number of unique T seqeunces
T_u<-length(TCRs_unique)
freq<-tabulate(factor(RNA_U_paste))
TCRs_unique<-TCRs_unique[order(freq,decreasing=TRUE)]
freq<-freq[order(freq,decreasing=TRUE)]
#save a CSV file with all unqiue TCRs and their frequency in order of frequency
freq_table<-paste(TCRs_unique,freq,sep=",")
write.table(freq_table,paste(output_path,sample_name,"_TCR",expt,".csv",sep=""),sep="\t",row.names=FALSE)


###############################################################
#collect summary file data
output_sum[i,1]<-sample_name
output_sum[i,2]<-n
output_sum[i,3]<-n_u
output_sum[i,4]<-mean(bc_f)
output_sum[i,5]<-var(bc_f)
output_sum[i,6]<-RNA_u
output_sum[i,7]<-T_u


#end of file loop
}
##############################################################################
##############################################################################
#OUTPUT

#save tab delimited spreadsheet  file with summary of all seqeunce data for all files
colnames(output_sum)<-c("Sample name", "Total reads", "Reads_with_correct_barcode","Mean_BC_freq","Var_BC_freq","TCR_RNAs","Unique_TCRs")
write.table(output_sum,paste(output_path,"summary",expt,".txt",sep=""),sep="\t",row.names=FALSE)
#save sequence list file in relevant directory
output_all<-paste(output_path,"TCR_seq",sep="")
save(output_seq,file=output_all)

#if necessayr can reload file
load(file="FILE")
