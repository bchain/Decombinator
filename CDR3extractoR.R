#scrpt for working out aa translation and CDR3 from
#5 part Decombinator output
library(Biostrings)
#at work Cruciform
dropbox<-"D:/dropbox/"
#at home
dropbox<-"C:/Users/Benny Chain/dropbox/"

#folder with Collapsinator output (input for CDR3 extractor)
folder<-"TcR_seqs/Expt42/results_1_9_42.fastq/"
#Decombinator output file
file<-"1_9_42.fastq_beta.txt"

#output folder
output_folder<-paste(dropbox,"R/15_06_2015/",sep="")
#output file
file_out<-paste(substr(file,1,nchar(file)-4),"_CDR3.txt",sep="")
output<-paste(output_folder,file_out,sep="")
#data_path<-"D:/Dropbox/TcR_seqs/Expt42/results_1_9_42.fastq/"
file_in<-paste(dropbox,folder,file,sep="")
#ids is file with five part identifier (output of Decombinator or Collapsinator)

ids<-read.table(file_in,sep=",",stringsAsFactors=FALSE)
################################################################
################################################################
#The next session loads all the different V and J region files
#from github/decombinator and tehn works out where the C and FGXG motifs are 
#as a distnace from each end.
#load V and J regions
AV_in<-paste(dropbox,"github/decombinator/","human_TRAV_region.fasta",sep="")
AV = readDNAStringSet(AV_in)
AJ_in<-paste(dropbox,"github/decombinator/","human_TRAJ_region.fasta",sep="")
AJ = readDNAStringSet(AJ_in)
BV_in<-paste(dropbox,"/github/decombinator/","human_TRBV_region.fasta",sep="")
BV = readDNAStringSet(BV_in)
BJ_in<-paste(dropbox,"github/decombinator/","human_TRBJ_region.fasta",sep="")
BJ = readDNAStringSet(BJ_in)

#find positions of Cs in V regions
AV_Cs<-gregexpr("C",as.character(translate(AV)))
AV_C<-sapply(AV_Cs,max)
#the C in TRDV2*01 is the penultimate one
AV_C[46]<-92
#same for beta chains
BV_Cs<-gregexpr("C",as.character(translate(BV)))
BV_C<-sapply(BV_Cs,max)

#find  GXG in alpha J regions
lengths<-nchar(AJ)
AJ_sh<-subseq(AJ,start=(lengths - 48),end=(lengths - 1))
AJ_shaa<-as.character(translate(AJ_sh))
AJ_GXGs<-gregexpr("[F,W]G.G",AJ_shaa)
AJ_GXG<-sapply(AJ_GXGs,max)
#anomalous CDR3 motif in TRAJ16*0
AJ_GXG[7]<-6
AJ_GXG_e<-nchar(AJ_shaa)-AJ_GXG-3

#find GXG in beta chains
lengths<-nchar(BJ)
BJ_sh<-subseq(BJ,start=(lengths - 45),end=(lengths))
BJ_shaa<-as.character(translate(BJ_sh))
BJ_GXGs<-gregexpr("[F,W]G.G",BJ_shaa)
BJ_GXG<-sapply(BJ_GXGs,max)
BJ_GXG_e<-nchar(BJ_shaa)-BJ_GXG-3

##################################################################################################
###################################################################################################
#now analyse the data file and obtain teh full amino acid translation
#and also teh CDR3 sequences
if(substr(file_in,start=nchar(file_in)-8,stop=nchar(file_in)-4)=="alpha"){chain<-"alpha"}
if(substr(file_in,start=nchar(file_in)-7,stop=nchar(file_in)-4)=="beta"){chain<-"beta"}

if(chain=="alpha"){
 
    V<-subseq(AV[ids [,1]+1],start=1,end=width(AV[ids [,1]+1])-ids[,3])
    J<-subseq(AJ[ids [,2]+1],start=(ids[,4]+1),end=width(AJ[ids [,2]+1]))
    V1<-paste(V,gsub("\\s","",ids[,5]),J,sep="")
    VDJ<-suppressWarnings(translate(DNAStringSet(V1)))
    CDR3<-subseq(VDJ,start=AV_C[ids [,1]+1],end=width(VDJ)-AJ_GXG_e[(ids [,2]+1)])
                }

if(chain=="beta"){

V<-subseq(BV[ids [,1]+1],start=1,end=width(BV[ids [,1]+1])-ids[,3])
J<-subseq(BJ[ids [,2]+1],start=(ids[,4]+1),end=width(BJ[ids [,2]+1]))
V1<-paste(V,gsub("\\s","",ids[,5]),J,sep="")
VDJ<-suppressWarnings(translate(DNAStringSet(V1)))
CDR3<-subseq(VDJ,start=BV_C[ids [,1]+1],end=width(VDJ)-BJ_GXG_e[(ids[,2]+1)])
                }


ids_CDR3<-cbind(ids,as.vector(VDJ),as.vector(CDR3))
write.table(ids_CDR3, output,sep="\t")


