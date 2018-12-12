library(CircTest)


rm(list=ls(all=TRUE));


args <- commandArgs(TRUE);



CircRNACount <- read.delim(paste(args[1],"/CircRNACount",sep=""),header=T)
LinearCount <- read.delim(paste(args[1],"/LinearCount",sep=""),header=T)
CircCoordinates <- read.delim(paste(args[1],"/CircCoordinates",sep=""),header=T)

minimum_sample_number <- min(as.numeric(args[2]), as.numeric(args[3]))



CircRNACount_filtered <- Circ.filter(circ = CircRNACount, linear = LinearCount, Nreplicates = minimum_sample_number, filter.sample = minimum_sample_number, filter.count = 5, percentage = 0.1)

CircCoordinates_filtered <- CircCoordinates[rownames(CircRNACount_filtered),]
LinearCount_filtered <- LinearCount[rownames(CircRNACount_filtered),]

test = Circ.test(CircRNACount_filtered, LinearCount_filtered, CircCoordinates_filtered, group=c(rep(1,args[2]),rep(2,args[3])))



output_file_folder<-args[1];



setwd(output_file_folder);



a_vs_b<-paste(args[4],args[5],sep="_vs_");



file_name_circTest_csv<-paste(args[6],".STAR_DCC_circTest.",a_vs_b,".circTest.csv",sep="");
file_name_circTest_txt<-paste(args[6],".STAR_DCC_circTest.",a_vs_b,".circTest.txt",sep="");



write.csv(as.data.frame(test$summary_table),file=file_name_circTest_csv)
write.table(as.data.frame(test$summary_table),file=file_name_circTest_txt,sep="\t",quote = FALSE,row.names = TRUE)



