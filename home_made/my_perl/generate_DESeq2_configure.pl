#!/usr/bin/perl -w
# Program name:		generate_DESeq2_configure.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-03-30
# Last update:		2016-03-30



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"c1:s","c2:s","b1:s","b2:s","o:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-h	:	help and usage

	-c1	:	condition1
	-c2	:	condition2

	-b1	:	bam1_array
	-b2	:	bam2_array

	-o	:	DESeq2_configure_col_csv_file

USAGE



die $usage unless $opts{c1};
die $usage unless $opts{c2};

die $usage unless $opts{b1};
die $usage unless $opts{b2};

die $usage unless $opts{o};



#,condition
#R_BM87_IL7R.sorted.unique.bam,BoneMarrow_Progenitor_Adult_N
#R_BM91_IL7R.sorted.unique.bam,BoneMarrow_Progenitor_Adult_N
#R_BM92_IL7R.sorted.unique.bam,BoneMarrow_Progenitor_Adult_N
#R_FL_CS17_IL7R.sorted.unique.bam,FetalLiver_Progenitor_Early_Late_N
#R_FL_CS21_IL7R.sorted.unique.bam,FetalLiver_Progenitor_Early_Late_N
#R_FL_CS22_IL7R.sorted.unique.bam,FetalLiver_Progenitor_Early_Late_N



my $condition_1=$opts{c1};
my $condition_2=$opts{c2};



my @bam_1=split(/\;/,$opts{b1});
my @bam_2=split(/\;/,$opts{b2});



open Profile,">$opts{o}";

print Profile ",condition\n";

for (my $i=0;$i<@bam_1;$i++) {

print Profile "$bam_1[$i].sorted.unique.bam,$condition_1\n";

}# end of for



for (my $i=0;$i<@bam_2;$i++) {

print Profile "$bam_2[$i].sorted.unique.bam,$condition_2\n";

}# end of for

close Profile;


