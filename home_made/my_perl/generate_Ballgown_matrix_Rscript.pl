#!/usr/bin/perl -w
# Program name:		generate_Ballgown_matrix_Rscript.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-07-30
# Last update:		2016-07-30



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"p:s","d:s","o:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-h	:	help and usage
	-p	:	project_name
	-d	:	directory
	-o	:	R_script

USAGE

die $usage unless $opts{p};
die $usage unless $opts{d};
die $usage unless $opts{o};



my $project_name=$opts{p};
my $directory=$opts{d};



my %group=();
my %replicate=();



open Profile,">$opts{o}";

{

###################################################################################################
###################################################################################################

print Profile 'library(ballgown)

library(gdata)

data_directory<-"';

print Profile "$directory";


print Profile '/Gene_expression_matrix/Ballgown_input";

bg = ballgown(dataDir=data_directory, samplePattern=\'R_\', meas=\'all\')

FPKM_table_transcript_level = texpr(bg, \'all\')







selected_column_name=matchcols(FPKM_table_transcript_level, with=c("t_name", "gene_name","FPKM"), method="or")

FPKM_table_transcript_level_selected=FPKM_table_transcript_level[,c(selected_column_name[[1]],selected_column_name[[2]],selected_column_name[[3]])]



write.table(as.matrix(FPKM_table_transcript_level_selected),file="';



print Profile "$directory";



print Profile '/Gene_expression_matrix/';


print Profile "$project_name";


print Profile '.HISAT_StringTie_Ballgown.transcript.FPKM.txt",sep="\t",quote = FALSE,row.names = FALSE)

FPKM_table_gene_level = gexpr(bg)


write.table(as.matrix(FPKM_table_gene_level),file="';



print Profile "$directory";


print Profile '/Gene_expression_matrix/';


print Profile "$project_name";


print Profile '.HISAT_StringTie_Ballgown.gene.FPKM.txt",sep="\t",quote = FALSE,row.names = TRUE)



';



}

close Profile;


