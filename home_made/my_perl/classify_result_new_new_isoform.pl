#!/usr/bin/perl -w
# Program name:		classify_result_new_new.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-03-30
# Last update:		2016-03-30



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"i:s","j:s","o1:s","o2:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-h	:	help and usage
	-i	:	Refgenes.gtf
	-j	:	cuffdiff_cufflink.combine
	-o1	:	cuffdiff_cufflink.combine.class
	-o2	:	cuffdiff_cufflink.combine.signature



USAGE

die $usage unless $opts{i};
die $usage unless $opts{j};


my $sign="";

my @content=();
my %class=();

#chrY	unknown	exon	59358329	59359508	.	-	.	gene_id "DDX11L16"; gene_name "DDX11L16"; transcript_id "NR_110561"; tss_id "TSS3246";

open I,"$opts{i}";

while (<I>)
{

chomp;

my @t=split(/\t/,$_);
my @tt=split(/gene_id \"/,$t[8]);
my @ttt=split(/\"/,$tt[1]);

my @tt_1=split(/transcript_id \"/,$t[8]);
my @ttt_1=split(/\"/,$tt_1[1]);



if (substr($ttt_1[0],0,2) eq "NM") {$class{$ttt[0]}=1;}

}
close I;


#my $output_1=$opts{j}.".class";

open Profile,">$opts{o1}";

open J,"$opts{j}";

while (<J>)
{

chomp;

my @t=split(/\t/,$_);

my @tt=split(/\:/,$t[0]);

if ((defined $tt[1])&&($tt[1] eq "isoform")) {my $lin0="Class\t".$_;print Profile "$lin0\n";push(@content,$lin0);}

else

{

if (defined $class{$t[1]}) {my $lin1="Protein-coding\t".$_;print Profile "$lin1\n";push(@content,$lin1);}

else {my $lin2="ncRNA\t".$_;print Profile "$lin2\n";push(@content,$lin2);}

}

}
close J;

close Profile;

#my $output_2=$opts{j}.".signature";

open Profile2,">$opts{o2}";

{

for (my $i=0;$i<@content;$i++) {

my @t=split(/\t/,$content[$i]);
my @tt=split(/\:/,$t[0]);

$sign=0;

if ($t[0] eq "Class") {print Profile2 "$content[$i]\n";}

else {

for (my $i=0;$i<@t;$i++) {

if ($t[$i] eq "yes") {$sign=1;}

}

if ($sign==1) {print Profile2 "$content[$i]\n";}


}  # end of else

}  #  end of for

}

close Profile2;
