#!/usr/bin/perl -w
# Program name:		get_gene_names.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-03-30
# Last update:		2016-03-30



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"i:s","o:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-h	:	help and usage
	-i	:	R_projectname.cuffdiff_cufflink.combine.class.pc.FPKM
	-o	:	hg19.gene.names/mm9.gene.names


USAGE

die $usage unless $opts{i};
die $usage unless $opts{o};



open Profile,">$opts{o}";



open I,"$opts{i}";

while (<I>)

{

chomp;

my @t=split(/\t/,$_);

print Profile "$t[0]\n";

}

close I;



close Profile;


