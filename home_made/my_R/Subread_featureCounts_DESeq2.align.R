
library(Rsubread)


rm(list=ls(all=TRUE));


args <- commandArgs(TRUE);


if (args[1]=="Paired-End") {

readfile1_name<-paste(args[3],"/Processed_FASTQ/",args[2],"_1.fastq",sep="");
readfile2_name<-paste(args[3],"/Processed_FASTQ/",args[2],"_2.fastq",sep="");
output_file_name<-paste(args[3],"/Subread_featureCounts_DESeq2/Processed_BAM/",args[2],".unique.bam",sep="");

align(index=args[4],readfile1=readfile1_name,readfile2=readfile2_name,input_format="gzFASTQ",output_format="BAM",output_file=output_file_name,unique=TRUE,indels=5,minFragLength=args[5],maxFragLength=args[6],PE_orientation=args[7]);

} else {


readfile1_name<-paste(args[3],"/Processed_FASTQ/",args[2],".adapt.qc.fastq",sep="");
output_file_name<-paste(args[3],"/Subread_featureCounts_DESeq2/Processed_BAM/",args[2],".unique.bam",sep="");

align(index=args[4],readfile1=readfile1_name,input_format="gzFASTQ",output_format="BAM",output_file=output_file_name,unique=TRUE,indels=5);

}


