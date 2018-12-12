#!/usr/bin/perl -w
# Program name:		select_cuffcompare_class_code.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-12-01
# Last update:		2016-12-01



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
	-i	:	cuffcmp.combined.gtf
	-o	:	cuffcmp.combined.selected.gtf



USAGE

die $usage unless $opts{i};
die $usage unless $opts{o};



#chr1	Cufflinks	exon	29554	30039	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "1"; gene_name "MIR1302-11"; oId "TCONS_00000001"; nearest_ref "ENST00000473358.1"; class_code "="; tss_id "TSS1";
#chr1	Cufflinks	exon	30564	30667	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "2"; gene_name "MIR1302-11"; oId "TCONS_00000001"; nearest_ref "ENST00000473358.1"; class_code "="; tss_id "TSS1";



open Profile,">$opts{o}";

open I,"$opts{i}";

while (<I>)
{

chomp;

my @t=split(/\t/,$_);
my @tt=split(/class_code \"/,$t[8]);
my @ttt=split(/\"/,$tt[1]);

if (($ttt[0] eq "i")||($ttt[0] eq "j")||($ttt[0] eq "o")||($ttt[0] eq "u")||($ttt[0] eq "x"))

{

print Profile "$_\n";

}

}
close I;

close Profile;


