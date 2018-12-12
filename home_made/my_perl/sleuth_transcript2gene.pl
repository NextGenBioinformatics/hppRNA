#!/usr/bin/perl -w
# Program name:		sleuth_transcript2gene.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-05
# Last update:		2016-06-05



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
	-j	:	R_hppRNA_data.Kallisto_sleuth.transcript.TPM.txt
	-o	:	R_hppRNA_data.Kallisto_sleuth.gene.TPM.txt


USAGE

die $usage unless $opts{i};
die $usage unless $opts{j};
die $usage unless $opts{o};


#NM_032291       SGIP1
#NM_001308203    SGIP1
#NM_001080397    SLC45A1


my %gene=();
my %gene_name=();


open I,"$opts{i}";

while (<I>)
{

chomp;

my @t=split(/\t/,$_);

$gene{$t[0]}=$t[1];

$gene_name{$t[1]}=1;

}
close I;


#transcript_id R_brain_3b      R_brain_3c      R_brain_a       R_testis_7a     R_testis_7b     R_testis_7c
#NM_001287050.1  0       0       0       4.9186  3.01216 2.49596
#NM_001039211.2  1.96607 2.94101 0.782442        2.08198 1.8114  2.14341


my $column_number="";
my %gene_expression=();


open Profile,">$opts{o}";

open J,"$opts{j}";

while (<J>)

{

chomp;

my @t=split(/\t/,$_);

$column_number=@t;

if ($t[0] eq "transcript_id") {

print Profile "gene_id\t";

for (my $i=2;$i<@t-1;$i++) {print Profile "$t[$i]\t";}

print Profile "$t[(@t-1)]\n";

} #if ($t[0] eq "transcript_id") {




else

{

for (my $i=2;$i<@t;$i++) {

#my @tt=split(/\./,$t[0]);

if (defined $gene_expression{$gene{$t[0]}}{$i}) {$gene_expression{$gene{$t[0]}}{$i}=$gene_expression{$gene{$t[0]}}{$i}+$t[$i];}

else {$gene_expression{$gene{$t[0]}}{$i}=$t[$i];}

}#for (my $i=1;$i<@t;$i++) {

}#else

}

close J;



foreach my $key (keys %gene_name) {

print Profile "$key\t";

for (my $i=2;$i<$column_number-1;$i++) {

print Profile "$gene_expression{$key}{$i}\t";

}#for (my $i=1;$i<$column_number-1;$i++) {

print Profile "$gene_expression{$key}{($column_number-1)}\n";

}#foreach my $key (keys %gene_name) {



close Profile;


