#!/usr/bin/perl -w
# Program name:		generate_DCC_samplesheet_all.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-02
# Last update:		2016-06-02



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"d:s","s:s","n:s","o:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-h	:	help and usage
	-d	:	directory
	-s	:	all_sample_array
	-n	:	type of mapping
	-o	:	samplesheet

USAGE

die $usage unless $opts{d};
die $usage unless $opts{s};
die $usage unless $opts{n};
die $usage unless $opts{o};



my $directory=$opts{d};
my $all_sample_array=$opts{s};
my $mapping_type=$opts{n};



my @all_sample_array_t=split(/\;/,$all_sample_array);



open Profile,">$opts{o}";

{



for (my $i=0;$i<@all_sample_array_t;$i++)

{

print Profile $directory."/".$all_sample_array_t[$i]."/".$mapping_type."_mapping/".$all_sample_array_t[$i].".Chimeric.out.junction\n";

}



}

close Profile;



