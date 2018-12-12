library(Rsubread)
library(limma)
library(edgeR)



rm(list=ls(all=TRUE));



args <- commandArgs(TRUE);



bam_file_folder<-paste(args[7],"/Subread_featureCounts_DESeq2/Processed_BAM",sep="");


setwd(bam_file_folder);



full_1<-strsplit(args[4],"\\;");
full_2<-strsplit(args[5],"\\;");

full<-c(full_1[[1]],full_2[[1]]);


bam_file_name<-paste(full,".sorted.unique.bam",sep="");



if (args[1]=="Paired-End") {

fc<-featureCounts(files=bam_file_name,annot.ext=args[6],isGTFAnnotationFile=TRUE,GTF.featureType="exon",isPairedEnd=TRUE,nthreads=args[9],strandSpecific=args[10],checkFragLength=TRUE,minFragLength=args[11],maxFragLength=args[12],GTF.attrType="gene_id",useMetaFeatures=TRUE)

} else {

fc<-featureCounts(files=bam_file_name,annot.ext=args[6],isGTFAnnotationFile=TRUE,GTF.featureType="exon",isPairedEnd=FALSE,nthreads=args[9],strandSpecific=args[10],GTF.attrType="gene_id",useMetaFeatures=TRUE)

}



x <- DGEList(counts=fc$counts, genes=fc$annotation)
x_rpkm <- rpkm(x,x$genes$Length,log=FALSE)



output_file_folder<-paste(args[7],"/Subread_featureCounts_DESeq2/DEG",sep="");


setwd(output_file_folder);





a_vs_b<-paste(args[2],args[3],sep="_vs_");



file_name_x_csv<-paste(args[8],".Subread_featureCounts_DESeq2.",a_vs_b,".gene.counts.csv",sep="");

file_name_x_txt<-paste(args[8],".Subread_featureCounts_DESeq2.",a_vs_b,".gene.counts.txt",sep="");

file_name_x_rpkm_csv<-paste(args[8],".Subread_featureCounts_DESeq2.",a_vs_b,".gene.RPKM.csv",sep="");

file_name_x_rpkm_txt<-paste(args[8],".Subread_featureCounts_DESeq2.",a_vs_b,".gene.RPKM.txt",sep="");



write.csv(as.matrix(x),file=file_name_x_csv)
write.table(as.matrix(x),file=file_name_x_txt,sep="\t",quote = FALSE,row.names = TRUE)

write.csv(as.matrix(x_rpkm),file=file_name_x_rpkm_csv)
write.table(as.matrix(x_rpkm),file=file_name_x_rpkm_txt,sep="\t",quote = FALSE,row.names = TRUE)



library(DESeq2)



output_file_folder<-paste(args[7],"/Subread_featureCounts_DESeq2/DEG",sep="");


setwd(output_file_folder);



file_name_col<-paste(args[8],".Subread_featureCounts_DESeq2.",a_vs_b,".col.csv",sep="");


countData<-read.csv(file=file_name_x_csv,sep=",",header=TRUE,row.names=1)
colData<-read.csv(file=file_name_col,sep=",",header=TRUE,row.names=1)

dds<-DESeqDataSetFromMatrix(countData= countData,colData= colData,design=~condition)



dds$condition<-relevel(dds$condition, ref = args[2])



dds<-DESeq(dds)
res<-results(dds)
resOrdered<-res[order(res$padj),]



file_name_DESeq2_csv<-paste(args[8],".Subread_featureCounts_DESeq2.",a_vs_b,".DESeq2.csv",sep="");
file_name_DESeq2_txt<-paste(args[8],".Subread_featureCounts_DESeq2.",a_vs_b,".DESeq2.txt",sep="");



write.csv(as.data.frame(resOrdered),file=file_name_DESeq2_csv)
write.table(as.data.frame(resOrdered),file=file_name_DESeq2_txt,sep="\t",quote = FALSE,row.names = TRUE)

head(resOrdered)


