#!/usr/bin/perl -w
# Program name:		Ballgown_fix_FPKM_title.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-05-31
# Last update:		2016-05-31



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
	-i	:	R_hppRNA_data.HISAT_StringTie_Ballgown.gene.FPKM.txt
	-o	:	R_hppRNA_data.HISAT_StringTie_Ballgown.gene.FPKM.clean.txt


USAGE

die $usage unless $opts{i};
die $usage unless $opts{o};

#FPKM.R_brain_3b	FPKM.R_brain_3c	FPKM.R_brain_a	FPKM.R_testis_7a	FPKM.R_testis_7b	FPKM.R_testis_7c

open Profile,">$opts{o}";

my $line=1;

open I,"$opts{i}";

while (<I>)

{

chomp;

my @t=split(/\t/,$_);

if ($line==1)

{

for (my $i=0;$i<@t-1;$i++)

{

my @tt=split(/FPKM./,$t[$i]);

print Profile "$tt[1]\t";

}

my @ttt=split(/FPKM./,$t[@t-1]);

print Profile "$ttt[1]\n";

$line++;

}

else

{

print Profile "$_\n";

}


}

close I;



close Profile;


