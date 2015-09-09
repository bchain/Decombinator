#script to collapse seqeunce data to single barcodes and collect numbers 
#this script is specifically for ligation to SP26N
 
#initialise file to contain output

output_sum<-c()
output_seq<-list()
#from my work computer
dropbox<-"D:/dropbox/"
#at home
dropbox<-"C:/Users/Benny Chain/dropbox/"

#Data input 
expt_path<-"C:/TCR/input/"
#output
output_path<-"C:/TCR/output/"
#expt number
expt<-45

####################################################################################################################
####################################################################################################################

directory<-dir(expt_path)
output_sum<-matrix(rep(NA,4*length(directory)),ncol=4)
#cycle through directory analysing each file in turn and 
#extracting unique TCR RNAs (i.e. collaping on bar codes) and 
#then unique TCRs and storing the 6 part identifier (including CDR3)
#in list file.
i<-1
for( i in 1:length(directory)){
sample_name<-directory[i]
sample_name
#read in data 
data_path<-paste(expt_path,sample_name,sep="")
TCRs<-read.table(data_path,stringsAsFactors=FALSE,header=FALSE,sep=",")
colnames(TCRs)<-c("V","J","Vdel","Jdel","insert","seqid","barcode","X")

#look for correct SP2 barcode when suitbale index present
#index<-"GTCGTGAT......GTCGTGAT......"
#number of reads
n<-dim(TCRs)[1]
n

TCRs_bc<-TCRs
TCRs_bc[,7]<-substr(TCRs[,7],1,7)

#look for distribtuion of barcode freqeuncies
bc_f<-tabulate(factor(TCRs_bc$barcode))
#mean and variance of barcode copy number
mean(bc_f)
var(bc_f)
#plot and save histogram of barcode family size
#hist(bc_f,breaks=c(0,1,10,100,max(bc_f)))$counts
#output<-paste(output,"barcode_f",sample_name,".png",sep="")
#imageSave(output)

#check for unique TCR RNA sequences (i.e. unique combination of barcode and TCR)
TCRs_M<-paste(TCRs_bc[,1],TCRs_bc[,2],TCRs_bc[,3],TCRs_bc[,4],TCRs_bc[,5],TCRs_bc[,7],sep=",")
TCRS_U<-levels(factor(TCRs_M))
#number of unique TCR RNA molecules
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

###############################################################
#collect summary file data
output_sum[i,1]<-sample_name
output_sum[i,2]<-n
output_sum[i,3]<-RNA_u
output_sum[i,4]<-T_u
#end of file loop
}
##############################################################################
##############################################################################
#OUTPUT

#save tab delimited spreadsheet  file with summary of all seqeunce data for all files
colnames(output_sum)<-c("Sample name", "Total reads", "Unique_RNAs","Unique_TCRs")
write.table(output_sum,paste(output_path,"summary",expt,".txt",sep=""),sep="\t",row.names=FALSE)
#save sequence list file in relevant directory
output_all<-paste(output_path,"TCR_seq",sep="")
save(output_seq,file=output_all)

#if necessayr can reload file
load(paste(dropbox,"R/26_05_2015/","processed_data_115",sep=""))
#save each set of results as a separate tab delimited file.
