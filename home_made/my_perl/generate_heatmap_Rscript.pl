#!/usr/bin/perl -w
# Program name:		generate_heatmap_Rscript.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-03-30
# Last update:		2016-03-30



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"i:s","d:s","o:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-h	:	help and usage
	-i	:	data_file_name
	-d	:	folder_name
	-o	:	R_script


USAGE

die $usage unless $opts{i};
die $usage unless $opts{d};
die $usage unless $opts{o};





my $input_data_full_path=$opts{d}."/".$opts{i};
my $output_figure_full_path=$opts{d}."/".$opts{i}.".1.heatmap.pdf";



open Profile,">$opts{o}";



###########################################################################################################
print Profile '

rm(list=ls(all=TRUE))

pdf(file="';


print Profile "$output_figure_full_path";

print Profile '",width=6,height=5,title="",pointsize=3)



on.exit(par(old.par));

old.par<-par(mar=c(2.5,2.5,1,2.5),mgp=c(1.5,0.5,0),fin=c(3.5,3),tck=-0.02,mfrow=c(1,1),bty="l",cex.lab=1,cex.axis=1,cex.main=1,font.main=3,font.lab=2,font.axis=2);

data = read.table("';



print Profile "$input_data_full_path";


print Profile '",header=TRUE,stringsAsFactors=F,sep="\t",quote="");    # Input the data. 

library(gplots);




data1 = data

my.rows <- apply(data1, 1, max) >=1;

data1 <- data1[my.rows, 1: ncol(data1)]





data2 = log2(data1+1);


library(preprocessCore)
data3 = normalize.quantiles(as.matrix(data2)) 


rownames(data3) = rownames(data2)
colnames(data3) = colnames(data2)


heatmap.2(data3, col=greenred(256), dendrogram="both",scale="row",trace="none",margins = c(15,10));

dev.off();


';



close Profile;


