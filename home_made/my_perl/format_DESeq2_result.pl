#!/usr/bin/perl -w
# Program name:		format_DESeq2_result.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-03-30
# Last update:		2016-03-30



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"i:s","j:s","b:s","c:s","o:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-h	:	help and usage
	-i	:	Species_specific_gene_name_list
	-j	:	DESeq.csv
	-b	:	baseline_condition_name
	-c	:	case_condition_name
	-o	:	DESeq.pc.csv


USAGE

die $usage unless $opts{i};
die $usage unless $opts{j};
die $usage unless $opts{b};
die $usage unless $opts{c};
die $usage unless $opts{o};



my %is=();

open I,"$opts{i}";

while (<I>)

{

chomp;

my @t=split(/\t/,$_);

if ($t[0] ne "Gene") {$is{$t[0]}=1;}

}

close I;


open Profile,">$opts{o}";

open J,"$opts{j}";

while (<J>)

{

chomp;

my @t=split(/,/,$_);

if ($t[0] eq "\"\"") {print Profile "$_,\"direction\",\"baseline\",\"comparison\"\n";}

my @tt=split(/\"/,$t[0]);

if ((defined $tt[1])&&(defined $is{$tt[1]})&&($t[6] ne "NA")) {



if (($t[6]<0.05)&&($t[2]>0)) {print Profile "$_,\"Up\",\"$opts{b}\",\"$opts{c} vs $opts{b}\"\n";}



if (($t[6]<0.05)&&($t[2]<0)) {print Profile "$_,\"Down\",\"$opts{b}\",\"$opts{c} vs $opts{b}\"\n";}



}#if ((defined $tt[1])&&(defined $is{$tt[1]})) {

}

close J;

close Profile;



