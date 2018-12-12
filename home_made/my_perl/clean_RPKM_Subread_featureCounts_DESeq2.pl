#!/usr/bin/perl -w
# Program name:		clean_RPKM_Subread_featureCounts_DESeq2.pl
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
	-i	:	R_hppRNA_data.Subread_featureCounts_DESeq2.gene.RPKM.txt
	-o	:	R_hppRNA_data.Subread_featureCounts_DESeq2.gene.RPKM.clean.txt


USAGE

die $usage unless $opts{i};
die $usage unless $opts{o};



#R_brain_3b.sorted.unique.bam    R_brain_3c.sorted.unique.bam    R_brain_a.sorted.unique.bam     R_testis_7a.sorted.unique.bam   R_testis_7b.sorted.unique.bam   R_testis_7c.sorted.unique.bam



my $line=1;



open Profile,">$opts{o}";

open I,"$opts{i}";

while (<I>)

{

chomp;

my @t=split(/\t/,$_);

if ($line==1) {

for (my $i=0;$i<@t-1;$i++) {

my @tt=split(/\.sorted.unique.bam/,$t[$i]);

print Profile "$tt[0]\t";

}   #####for end


my @ttt=split(/\.sorted.unique.bam/,$t[(@t-1)]);

print Profile "$ttt[0]\n";


$line++;


}  ###if ($t[0] eq "Gene") {

else {print Profile "$_\n";}

} #end of while

close I;

close Profile;
