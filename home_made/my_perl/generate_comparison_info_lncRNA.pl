#!/usr/bin/perl -w
# Program name:		generate_comparison_info_lncRNA.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-02
# Last update:		2016-06-02



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"p:s","d:s","s:s","g:s","c:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-h	:	help and usage
	-p	:	project_name
	-d	:	directory
	-s	:	sample_name_array
	-g	:	group_name_array
	-c	:	comparison

USAGE

die $usage unless $opts{p};
die $usage unless $opts{d};
die $usage unless $opts{s};
die $usage unless $opts{g};
die $usage unless $opts{c};



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

###########################################################################################################



for (my $i=0;$i<@comparison_t;$i++)

{

my @t=split(/\:/,$comparison_t[$i]);

my $new_name=$t[0]."_vs_".$t[1];

my $info_file_full_name=$directory."/".$project_name.".lncRNA_denovo.".$new_name.".info.txt";

open Profile,">$info_file_full_name";

{

print Profile "run_accession condition\n";

for (my $j=0;$j<@{$replicate{$t[0]}};$j++)

{

print Profile $replicate{$t[0]}[$j]." ".$group{$replicate{$t[0]}[$j]}."\n";

}



for (my $j=0;$j<@{$replicate{$t[1]}};$j++)

{

print Profile $replicate{$t[1]}[$j]." ".$group{$replicate{$t[1]}[$j]}."\n";

}

}

close Profile;


}


