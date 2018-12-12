#!/usr/bin/perl -w
# Program name:		only_FPKM.pl
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
	-o	:	cuffdiff_cufflink.combine.class/signature.FPKM


USAGE

die $usage unless $opts{i};

die $usage unless $opts{o};

my @index=();
my @name=();

open Profile,">$opts{o}";

open I,"$opts{i}";

while (<I>)

{

chomp;

my @t=split(/\t/,$_);

if ($t[0] eq "Class") {

for (my $i=0;$i<@t;$i++) {

my @tt=split(/\:/,$t[$i]);

if ((defined $tt[1])&&($tt[1] eq "FPKM")) {push(@index,$i);push(@name,$tt[0]);}

} #for (my $i=0;$i<@t;$i++) {

print Profile "Gene\t";

for (my $i=0;$i<@index-1;$i++) {print Profile "$name[$i]\t";}

print Profile "$name[@index-1]\n";

}# if ($t[0] eq "Class") {

else 

{
print Profile "$t[1]\t";

for (my $i=0;$i<@index-1;$i++) {print Profile "$t[$index[$i]]\t";}

print Profile "$t[$index[@index-1]]\n";

}

} #end of while

close I;

close Profile;
