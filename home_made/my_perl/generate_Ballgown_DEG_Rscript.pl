#!/usr/bin/perl -w
# Program name:		generate_Ballgown_DEG_Rscript.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-07-30
# Last update:		2016-07-30



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"p:s","d:s","s:s","g:s","c:s","o:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-h	:	help and usage
	-p	:	project_name
	-d	:	directory
	-s	:	sample_name_array
	-g	:	group_name_array
	-c	:	comparison
	-o	:	R_script

USAGE

die $usage unless $opts{p};
die $usage unless $opts{d};
die $usage unless $opts{s};
die $usage unless $opts{g};
die $usage unless $opts{c};
die $usage unless $opts{o};



my $project_name=$opts{p};
my $directory=$opts{d};
my $sample_name_array=$opts{s};
my $group_name_array=$opts{g};
my $comparison=$opts{c};

my %group=();
my %replicate=();

my @sample_name_array_t=split(/\;/,$sample_name_array);

my @group_name_array_t=split(/\;/,$group_name_array);

my @comparison_t=split(/\;/,$comparison);

for (my $i=0;$i<@sample_name_array_t;$i++)

{

$group{$sample_name_array_t[$i]}=$group_name_array_t[$i];

push(@{$replicate{$group_name_array_t[$i]}},$sample_name_array_t[$i]);

}



open Profile,">$opts{o}";

{

###################################################################################################
###################################################################################################

print Profile 'library(ballgown)

library(gdata)



';



for (my $i=0;$i<@comparison_t;$i++)

{


my @t=split(/\:/,$comparison_t[$i]);

my $new_name=$t[0]."_vs_".$t[1];







print Profile 'data_directory<-"';


print Profile "$directory";

print Profile '/DEG/';

print Profile "$new_name";

print Profile '";

bg = ballgown(dataDir=data_directory, samplePattern=\'R_\', meas=\'all\')


';





my $group_real_array_1="\"".$replicate{$t[1]}[0]."\"";
my $group_number_array_1="1";

for (my $j=1;$j<@{$replicate{$t[1]}};$j++)

{

$group_real_array_1=$group_real_array_1.",\"".$replicate{$t[1]}[$j]."\"";
$group_number_array_1=$group_number_array_1.",1";

}




my $group_real_array_0="\"".$replicate{$t[0]}[0]."\"";
my $group_number_array_0="0";

for (my $j=1;$j<@{$replicate{$t[0]}};$j++)

{

$group_real_array_0=$group_real_array_0.",\"".$replicate{$t[0]}[$j]."\"";
$group_number_array_0=$group_number_array_0.",0";

}




print Profile 'group_real<-c(';


print Profile "$group_number_array_1";

print Profile ',';

print Profile "$group_number_array_0";


print Profile ')


names(group_real)<-c(';

print Profile "$group_real_array_1";

print Profile ',';

print Profile "$group_real_array_0";

print Profile ')



group_array=c();

for (i in 1:length(sampleNames(bg)))

{

group_array=c(group_array,group_real[sampleNames(bg)[i]][[1]]);

}



pData(bg) = data.frame(id=sampleNames(bg), group=group_array)



stat_gene_results = stattest(bg, feature=\'gene\', meas=\'FPKM\', covariate=\'group\', getFC = TRUE)



stat_gene_results_resOrdered <- stat_gene_results[order(stat_gene_results$qval),]



write.table(as.matrix(stat_gene_results_resOrdered),file="';



print Profile "$directory";


print Profile '/DEG/';

print Profile "$project_name";

print Profile '.HISAT_StringTie_Ballgown.';

print Profile "$new_name";

print Profile '.gene.DEG.txt",sep="\t",quote = FALSE,row.names = FALSE)



stat_transcript_results = stattest(bg, feature=\'transcript\', meas=\'FPKM\', covariate=\'group\', getFC = TRUE)



stat_transcript_results_merge=merge(stat_transcript_results, texpr(bg, \'all\'), by.x = "id", by.y = "t_id", all = TRUE)

stat_transcript_results_merge_selected=stat_transcript_results_merge[,c(2,10,14,3,4,5)]

stat_transcript_results_merge_selected_resOrdered <- stat_transcript_results_merge_selected[order(stat_transcript_results_merge_selected$qval),]

write.table(as.matrix(stat_transcript_results_merge_selected_resOrdered),file="';



print Profile "$directory";

print Profile '/DEG/';

print Profile "$project_name";

print Profile '.HISAT_StringTie_Ballgown.';

print Profile "$new_name";

print Profile '.transcript.DEG.txt",sep="\t",quote = FALSE,row.names = FALSE)



';



}



}

close Profile;


