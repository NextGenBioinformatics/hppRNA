#!/usr/bin/perl -w
# Program name:		generate_PCA_Rscript.pl
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
my $output_matrix_full_path=$opts{d}."/".$opts{i}.".1.matrix.txt";
my $output_rotation_full_path=$opts{d}."/".$opts{i}.".1.rotation.txt";
my $output_figure_full_path=$opts{d}."/".$opts{i}.".1.PCA.pdf";



open Profile,">$opts{o}";

print Profile '



rm(list=ls(all=TRUE))

library(preprocessCore)
library(ggfortify);
library(cluster);


pdf(file="';


print Profile "$output_figure_full_path";


print Profile '",width=10,height=10,title="",pointsize=8)



data = read.table("';


print Profile "$input_data_full_path";

print Profile '",header=TRUE,stringsAsFactors=F,sep="\t",quote="");    # Input the data. 



data1 = data

my.rows <- apply(data1, 1, max) >=1;

data1 <- data1[my.rows, 1: ncol(data1)]



write.table(data1,file="';


print Profile "$output_matrix_full_path";

print Profile '",sep="\t",quote = FALSE,row.names = TRUE)



data2 = log2(data1+1);


data3 = normalize.quantiles(as.matrix(data2)) 


rownames(data3) = rownames(data2)
colnames(data3) = colnames(data2)

data4 = t(data3);



autoplot(prcomp(data4),data=data4,label=TRUE,label.size=3,label.colour="red");



fit <- prcomp(data4);

write.table(fit$rotation, file="';



print Profile "$output_rotation_full_path";


print Profile '",sep="\t",quote = FALSE,row.names = TRUE);



dev.off();


';


close Profile;


