#!/usr/bin/perl -w
# Program name:		combine_cuffdiff_cufflink_new.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-03-30
# Last update:		2016-03-30



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
	-i	:	cuffdiff.combine
	-j	:	cufflink.combine

	-o	:	cuffdiff_cufflink.combine


USAGE

die $usage unless $opts{i};
die $usage unless $opts{j};

die $usage unless $opts{o};

my %content=();
my $title="";


open I,"$opts{i}";

while (<I>)
{

chomp;

my @t=split(/\t/,$_);

my @tt=split(/\:/,$t[0]);

if ((defined $tt[1])&&($tt[1] eq "test_id")) {$title=$_;}

else {$content{$t[0]}=$_;}

}
close I;

open Profile,">$opts{o}";

open J,"$opts{j}";

while (<J>)
{

chomp;

my @t=split(/\t/,$_);

my @tt=split(/\:/,$t[0]);

if ((defined $tt[1])&&($tt[1] eq "isoform")) {my $r1=chop($_);my $r2=chop($title);print Profile "$_\t$title\n";}

else {my $r3=chop($_);my $r4=chop($content{$t[0]});print Profile "$_\t$content{$t[0]}\n";}

}
close J;

close Profile;
