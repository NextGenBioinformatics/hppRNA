


library(edgeR)



rm(list=ls(all=TRUE));



args <- commandArgs(TRUE);



counts_file_name<-paste(args[5],"/Bowtie_eXpress_edgeR/DEG/",args[6],".Bowtie_eXpress_edgeR.",args[1],"_vs_",args[2],".eff_counts.txt",sep="");



counts = read.table(counts_file_name,header=TRUE,stringsAsFactors=FALSE,sep="\t",quote="");




full_1<-strsplit(args[3],"\\;");
full_2<-strsplit(args[4],"\\;");

full<-c(full_1[[1]],full_2[[1]]);



a_1<-replicate(length(full_1[[1]]),args[1]);
a_2<-replicate(length(full_2[[1]]),args[2]);

a<-c(a_1,a_2);


names(a)<-full;



group<-unname(a[c(colnames(counts))]);



d <- DGEList(counts=counts, group=group)

d <- calcNormFactors(d)

d <- estimateCommonDisp(d)

d <- estimateTagwiseDisp(d)

et <- exactTest(d,pair=c(args[1],args[2]))

final<-topTags(et,n=100000)



file_name_x_csv<-paste(args[5],"/Bowtie_eXpress_edgeR/DEG/",args[6],".Bowtie_eXpress_edgeR.",args[1],"_vs_",args[2],".transcript.DEG.csv",sep="");

file_name_x_txt<-paste(args[5],"/Bowtie_eXpress_edgeR/DEG/",args[6],".Bowtie_eXpress_edgeR.",args[1],"_vs_",args[2],".transcript.DEG.txt",sep="");



write.csv(as.data.frame(final),file=file_name_x_csv)

write.table(as.data.frame(final),file=file_name_x_txt,sep="\t",quote = FALSE,row.names = TRUE,col.names=NA)


