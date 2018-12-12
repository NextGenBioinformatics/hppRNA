
library(Rsubread)
library(limma)
library(edgeR)


rm(list=ls(all=TRUE));


args <- commandArgs(TRUE);



bam_file_folder<-paste(args[2],"/Subread_featureCounts_DESeq2/Processed_BAM",sep="");


setwd(bam_file_folder);


full<-strsplit(args[4],"\\;");



bam_file_name<-paste(full[[1]],".sorted.unique.bam",sep="");



if (args[1]=="Paired-End") {

fc<-featureCounts(files=bam_file_name,annot.ext=args[3],isGTFAnnotationFile=TRUE,GTF.featureType="exon",isPairedEnd=TRUE,nthreads=args[6],strandSpecific=args[7],checkFragLength=TRUE,minFragLength=args[8],maxFragLength=args[9],GTF.attrType="gene_id",useMetaFeatures=TRUE)

} else {

fc<-featureCounts(files=bam_file_name,annot.ext=args[3],isGTFAnnotationFile=TRUE,GTF.featureType="exon",isPairedEnd=FALSE,nthreads=args[6],strandSpecific=args[7],GTF.attrType="gene_id",useMetaFeatures=TRUE)

}



x <- DGEList(counts=fc$counts, genes=fc$annotation)
x_rpkm <- rpkm(x,x$genes$Length,log=FALSE)



output_file_folder<-paste(args[2],"/Subread_featureCounts_DESeq2/Gene_expression_matrix",sep="");


setwd(output_file_folder);



file_name_x_csv<-paste(args[5],".Subread_featureCounts_DESeq2.gene.counts.csv",sep="");

file_name_x_txt<-paste(args[5],".Subread_featureCounts_DESeq2.gene.counts.txt",sep="");

file_name_x_rpkm_csv<-paste(args[5],".Subread_featureCounts_DESeq2.gene.RPKM.csv",sep="");

file_name_x_rpkm_txt<-paste(args[5],".Subread_featureCounts_DESeq2.gene.RPKM.txt",sep="");



write.csv(as.matrix(x),file=file_name_x_csv)
write.table(as.matrix(x),file=file_name_x_txt,sep="\t",quote = FALSE,row.names = TRUE)



write.csv(as.matrix(x_rpkm),file=file_name_x_rpkm_csv)
write.table(as.matrix(x_rpkm),file=file_name_x_rpkm_txt,sep="\t",quote = FALSE,row.names = TRUE)


