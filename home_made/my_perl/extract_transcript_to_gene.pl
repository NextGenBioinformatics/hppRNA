#!/usr/bin/perl -w
# Program name:		extract_transcript_to_gene.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-03
# Last update:		2016-06-03



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
	-i	:	hg19_refGene.txt
	-o	:	hg19.transcript_to_gene.txt


USAGE

die $usage unless $opts{i};
die $usage unless $opts{o};



open Profile,">$opts{o}";



open I,"$opts{i}";

while (<I>)

{

chomp;

my @t=split(/\t/,$_);

if (substr($t[0],0,1) ne "#")

{

print Profile "$t[1]\t$t[12]\n";

}

}

close I;



close Profile;


