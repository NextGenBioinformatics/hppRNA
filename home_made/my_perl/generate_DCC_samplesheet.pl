#!/usr/bin/perl -w
# Program name:		generate_DCC_samplesheet.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-02
# Last update:		2016-06-02



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"d:s","c:s","t:s","n:s","o:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-h	:	help and usage
	-d	:	directory
	-c	:	control_sample_array
	-t	:	treatment_sample_array
	-n	:	type of mapping
	-o	:	samplesheet

USAGE

die $usage unless $opts{d};
die $usage unless $opts{c};
die $usage unless $opts{t};
die $usage unless $opts{n};
die $usage unless $opts{o};



my $directory=$opts{d};
my $control_sample_array=$opts{c};
my $treatment_sample_array=$opts{t};
my $mapping_type=$opts{n};



my @control_sample_array_t=split(/\;/,$control_sample_array);

my @treatment_sample_array_t=split(/\;/,$treatment_sample_array);



open Profile,">$opts{o}";

{



for (my $i=0;$i<@control_sample_array_t;$i++)

{

print Profile $directory."/".$control_sample_array_t[$i]."/".$mapping_type."_mapping/".$control_sample_array_t[$i].".Chimeric.out.junction\n";

}



for (my $i=0;$i<@treatment_sample_array_t;$i++)

{

print Profile $directory."/".$treatment_sample_array_t[$i]."/".$mapping_type."_mapping/".$treatment_sample_array_t[$i].".Chimeric.out.junction\n";

}



}

close Profile;



