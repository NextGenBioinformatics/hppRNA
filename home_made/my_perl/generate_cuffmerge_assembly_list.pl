#!/usr/bin/perl -w
# Program name:		generate_cuffmerge_assembly_list.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-02
# Last update:		2016-06-02



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"d:s","s:s","o:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-h	:	help and usage
	-d	:	directory
	-s	:	all_sample_array
	-o	:	assembly_list.txt

USAGE

die $usage unless $opts{d};
die $usage unless $opts{s};
die $usage unless $opts{o};



my $directory=$opts{d};
my $all_sample_array=$opts{s};



my @all_sample_array_t=split(/\;/,$all_sample_array);



open Profile,">$opts{o}";

{



for (my $i=0;$i<@all_sample_array_t;$i++)

{

print Profile $directory."/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/".$all_sample_array_t[$i]."_clout/transcripts.gtf\n";

}



}

close Profile;



