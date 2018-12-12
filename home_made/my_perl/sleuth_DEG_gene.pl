#!/usr/bin/perl -w
# Program name:		sleuth_DEG_gene.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-04
# Last update:		2016-06-04



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"i:s","j:s","o:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-h	:	help and usage
	-i	:	hg19.transcript_to_gene.txt
	-j	:	R_hppRNA_data.Kallisto_sleuth.testis_vs_brain.transcript.DEG.txt

	-o	:	R_hppRNA_data.Kallisto_sleuth.testis_vs_brain.gene.DEG.txt


USAGE

die $usage unless $opts{i};
die $usage unless $opts{j};

die $usage unless $opts{o};

my %gene=();


open I,"$opts{i}";

while (<I>)
{

chomp;

my @t=split(/\t/,$_);

$gene{$t[0]}=$t[1];

}

close I;

open Profile,">$opts{o}";

open J,"$opts{j}";

while (<J>)
{

chomp;

my @t=split(/\t/,$_);

if ($t[0] eq "target_id") {print Profile "gene_id\t$_\n";}

else

{

#my @tt=split(/\./,$t[0]);

print Profile "$gene{$t[0]}\t$_\n";

}

}
close J;

close Profile;
