#!/usr/bin/perl -w
# Program name:		only_pc.pl
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
	-i	:	cuffdiff_cufflink.combine.class/signature
	-o	:	cuffdiff_cufflink.combine.class/signature.pc


USAGE

die $usage unless $opts{i};

die $usage unless $opts{o};



open Profile,">$opts{o}";

open I,"$opts{i}";

while (<I>)

{

chomp;

my @t=split(/\t/,$_);

print Profile "$_\n";

#if ($t[0] ne "ncRNA") {print Profile "$_\n";}

} #end of while

close I;

close Profile;
