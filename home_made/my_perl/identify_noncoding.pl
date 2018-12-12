#!/usr/bin/perl -w
# Program name:		identify_noncoding.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-12-01
# Last update:		2016-12-01



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
	-i	:	cuffcmp.combined.selected.result
	-j	:	cuffcmp.combined.selected.gtf
	-o	:	cuffcmp.combined.selected.noncoding.gtf



USAGE

die $usage unless $opts{i};
die $usage unless $opts{j};
die $usage unless $opts{o};



##id     type    score
#TCONS_00001761  coding  0.00450339

my %type=();

open I,"$opts{i}";

while (<I>)
{

chomp;

my @t=split(/\t/,$_);

if (($t[0] ne "#id")&&($t[1] eq "noncoding")) {

$type{$t[0]}=$t[1];

}

}
close I;



#chr1    Cufflinks       exon    317720  317781  .       +       .       gene_id "XLOC_000004"; transcript_id "TCONS_00000007"; exon_number "1"; gene_name "RP4-669L17.10"; oId "TCONS_00000008"; nearest_ref "ENST00000440038.2"; class_code "j"; tss_id "TSS5";



open Profile,">$opts{o}";

open J,"$opts{j}";

while (<J>)
{

chomp;

my @t=split(/\t/,$_);
my @tt=split(/oId \"/,$t[8]);
my @ttt=split(/\"/,$tt[1]);

if (defined $type{$ttt[0]}) {

print Profile "$_\n";

}

}
close J;

close Profile;


