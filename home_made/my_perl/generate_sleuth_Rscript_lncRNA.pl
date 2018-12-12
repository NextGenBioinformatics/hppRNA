#!/usr/bin/perl -w
# Program name:		generate_sleuth_Rscript_lncRNA.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-03
# Last update:		2016-06-03



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"p:s","d:s","c:s","o:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-h	:	help and usage
	-p	:	project_name
	-d	:	directory
	-c	:	comparison
	-o	:	R_script

USAGE

die $usage unless $opts{p};
die $usage unless $opts{d};
die $usage unless $opts{c};
die $usage unless $opts{o};



my $project_name=$opts{p};
my $directory=$opts{d};
my $comparison=$opts{c};


my @comparison_t=split(/\;/,$comparison);


open Profile,">$opts{o}";

{

###################################################################################################
###################################################################################################



print Profile 'library("sleuth")



rm(list=ls(all=TRUE))



base_dir <- "';


print Profile "$directory";


print Profile '/DEG"';




for (my $i=0;$i<@comparison_t;$i++)

{


my @t=split(/\:/,$comparison_t[$i]);

my $new_name=$t[0]."_vs_".$t[1];


print Profile '


sample_id <- dir(file.path(base_dir,"';


print Profile "$new_name";



print Profile '"))


sample_id

kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "';


print Profile "$new_name";



print Profile '", id, "kallisto"))

kal_dirs

s2c <- read.table(file.path(base_dir, "';


print Profile "$project_name";



print Profile '.lncRNA_denovo.';



print Profile "$new_name";


print Profile '.info.txt"), header = TRUE, stringsAsFactors=FALSE)

s2c <- dplyr::select(s2c, sample = run_accession, condition)
s2c

s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)

so <- sleuth_prep(s2c, ~ condition)

so <- sleuth_fit(so)

so <- sleuth_wt(so, \'condition';


print Profile "$t[0]";


print Profile '\')


models(so)


results_table <- sleuth_results(so, \'condition';



print Profile "$t[0]";


print Profile '\')


results_table_resOrdered <- results_table[order(results_table$qval),]


write.table(as.matrix(results_table_resOrdered),file=file.path(base_dir, "';


print Profile "$project_name";


print Profile '.lncRNA_denovo.';


print Profile "$new_name";


print Profile '.transcript.DEG.txt"),sep="\t",quote = FALSE,row.names = FALSE)';


}#for (my $i=0;$i<@comparison_t;$i++)

###################################################################################################
###################################################################################################



}

close Profile;


