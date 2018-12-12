#!/usr/bin/perl -w
# Program name:		clean_FPKM_STAR_RSEM_EBSeq.pl
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
	-i	:	R_hppRNA_data.STAR_RSEM_EBSeq.gene.FPKM.txt
	-o	:	R_hppRNA_data.STAR_RSEM_EBSeq.gene.FPKM.clean.txt


USAGE

die $usage unless $opts{i};
die $usage unless $opts{o};

#Gene	R_M3_D21_45RA.genes.fpkm_tracking

open Profile,">$opts{o}";

open I,"$opts{i}";

while (<I>)

{

chomp;

my @t=split(/\t/,$_);

if ($t[0] eq "gene_id") {

for (my $i=1;$i<@t-1;$i++) {

my @tt=split(/\t/,$t[$i]);

print Profile "$tt[0]\t";

}   #####for end


my @ttt=split(/\t/,$t[(@t-1)]);

print Profile "$ttt[0]\n";

}  ###if ($t[0] eq "gene_id") {

else {print Profile "$_\n";}

} #end of while

close I;

close Profile;
