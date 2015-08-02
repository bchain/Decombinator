#script to collapse seqeunce data to single barcodes and collect numbers
#initialise file to contain output
expt<-45
output<-list()
#path to expt 
#from my work computer
dropbox<-"D:/dropbox/"
#at home
dropbox<-"C:/Users/Benny Chain/dropbox/"

#Data input 
expt_path<-"C:/TCR/input/"
#output
output_path<-"c:/output/"

####################################################################################################################
####################################################################################################################
#cycle through directory analysing each file in turn and 
#extracting unique TCR RNAs (i.e. collaping on bar codes) and 
#then unique TCRs and storing the 6 part identifier (including CDR3)
#in list file.
i<-1
directory<-dir(expt_path)
sample_name<-directory[i]
sample_name
#read in data 

data_path<-paste(expt_path,sample_name,sep="")
TCRs<-read.table(data_path,stringsAsFactors=FALSE,header=FALSE,sep=",")

colnames(TCRs)<-c("V","J","Vdel","Jdel","insert","seqid","barcode","X")

#look for correct SP2 barcode when suitbale index present
index<-"GTCGTGAT......GTCGTGAT......"
#number of reads
n<-dim(TCRs)[1]
n
#number with correct index seqeunce; 
#use fuzzy matching to allow two missmatches
isf<-0
isf<-agrep(index,TCRs[,7],fixed=FALSE,max.distance=list(ins=0,del=0,sub=2))
length(isf)
length(isf)/n
#select only those  sequences with barcodes

if (isf != 0) {TCRs_bc<-TCRs[isf,]} else {TCRs_bc<-TCRs}

#look for distribtuion of barcode freqeuncies
bc_f<-tabulate(factor(TCRs_bc$barcode))
#mean and variance of barcode copy number
mean(bc_f)
var(bc_f)
#plot and save histogram of barcode family size
hist(bc_f,breaks=c(0,1,10,100,max(bc_f)))$counts
output<-paste(output,"barcode_f",sample_name,".png",sep="")
imageSave(output)

#check for unique TCR RNA sequences (i.e. unique combination of barcode and TCR)
TCRs_M<-paste(TCRs_bc[,1],TCRs_bc[,2],TCRs_bc[,3],TCRs_bc[,4],TCRs_bc[,5],TCRs_bc[,7],TCRs_bc[,9],sep=",")
TCRS_U<-levels(factor(TCRs_M))
#number of unique TCR RNA molecules
length(TCRS_U)

#split up TCRS_U into components again
RNA_U<-as.data.frame(matrix(unlist(strsplit(TCRS_U,split=",")),byrow=TRUE,ncol=7))
colnames(RNA_U)<-c("V","J","Vdel","Jdel","insert","bar_code","CDR3")
######################################################
#collect data for later analysis in list called e.g. TCR_115 (for patient 115)
TCR_115[[i]]<-list(RNA_U)
names(TCR_115)[i]<-sample_name
#save list file in relevant directory 
output_all<-paste(dropbox,"R/26_05_2015/","processed_data_115",sep="")
save(TCR_115,file=output_all)
i<-i+1
#}
########################################################
#look for unique T sequences
RNA_U_paste<-paste(RNA_U[,1],RNA_U[,2],RNA_U[,3],RNA_U[,4],RNA_U[,5],RNA_U[,7],sep=",")
TCRs_unique<-levels(factor(RNA_U_paste))
#number of unique T seqeunces
length(TCRs_unique)

#hist(tabulate(factor(RNA_U_paste)),breaks=c(0,1,10,100,1000,10000,100000))$counts
#zipf plot
plot(log(1:length(TCRs_unique)),log(sort(tabulate(factor(RNA_U_paste)),decreasing=TRUE)))
#split back into Decombinated elements
TCRs_u<-as.data.frame(matrix(unlist(strsplit(TCRs_unique,split=",")),byrow=TRUE,ncol=6))
colnames(TCRs_u)<-c("V","J","Vdel","Jdel","insert","CDR3")


#increment file counter
i<-i+1

#if necessayr can reload file
load(paste(dropbox,"R/26_05_2015/","processed_data_115",sep=""))
#save each set of results as a separate tab delimited file.
l<-length(TCR_115)
j<-1
file_name<-substring(names(TCR_115)[j],1, (nchar(names(TCR_115)[j])-6))
output_each<-paste(dropbox,"R/26_05_2015/",file_name,".txt",sep="")
for (j in 1:l){
  file_name<-substring(names(TCR_115)[j],1, (nchar(names(TCR_115)[j])-6))
  output_each<-paste(dropbox,"R/26_05_2015/",file_name,".txt",sep="")
  write.table(TCR_115[[j]],output_each,sep="\t",row.names=FALSE)
}