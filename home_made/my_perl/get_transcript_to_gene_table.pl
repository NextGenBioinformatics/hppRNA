#!/usr/bin/perl -w
# Program name:		get_transcript_to_gene_table.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-19
# Last update:		2016-06-19



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
	-i	:	genes.selected.pc.gtf
	-o	:	transcript_to_gene.txt



USAGE

die $usage unless $opts{i};
die $usage unless $opts{o};



my %class=();

#chrY	unknown	exon	59358329	59359508	.	-	.	gene_id "DDX11L16"; gene_name "DDX11L16"; transcript_id "NR_110561"; tss_id "TSS3246";



#NM_032291       SGIP1
#NM_001308203    SGIP1



my %gene=();

open I,"$opts{i}";

while (<I>)
{

chomp;

my @t=split(/\t/,$_);

if (substr($t[0],0,2) ne "##")

{

my @tt=split(/gene_id \"/,$t[8]);
my @ttt=split(/\"/,$tt[1]);

my @tt_1=split(/transcript_id \"/,$t[8]);
my @ttt_1=split(/\"/,$tt_1[1]);

$gene{$ttt_1[0]}=$ttt[0];

}

}
close I;



open Profile,">$opts{o}";
{
foreach my $key (keys %gene) {print Profile "$key\t$gene{$key}\n";}
}
close Profile;


