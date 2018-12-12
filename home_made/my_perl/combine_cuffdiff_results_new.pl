#!/usr/bin/perl -w
# Program name:		combine_cuffdiff_results.pl
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
	-i	:	folder
	-o	:	raw_result


USAGE

die $usage unless $opts{i};
die $usage unless $opts{o};

my $dir_to_process="$opts{i}";

my @file_list=();
my @gene_list=();

my %gene=();
my %FC=();
my %q_value=();
my %sig=();


##############################################################################
###############################begin to go through each file################

opendir DH, $dir_to_process or die "can't open $dir_to_process:$!";

while (my $name=readdir DH) {	

next if $name=~/^\./;

my $full_name ="$dir_to_process/$name";

# read file.

#test_id	gene_id	gene	locus	sample_1	sample_2	status	value_1	value_2	log2(fold_change)	test_stat	p_value	q_value	significant
#A1BG	A1BG	A1BG	chr19:58858171-58874214	TAR139G	Csi	OK	7.2083	7.617	0.0795633	0.238474	0.7319	0.997867	no

push(@file_list,$name);

open I,"$full_name";

#############################################################################################
while (<I>)
{

chomp;

my @t=split(/\t/,$_);

if ($t[0] ne "test_id") {

my $test_id=$t[0];
my $gene_id=$t[1];
my $gene=$t[2];
my $locus=$t[3];
my $sample_1=$t[4];
my $sample_2=$t[5];
my $status=$t[6];
my $value_1=$t[7];
my $value_2=$t[8];
my $log2_fold_change=$t[9];
my $test_stat=$t[10];
my $p_value=$t[11];
my $q_value=$t[12];
my $significant=$t[13];

my $core=$name.":".$test_id;

$gene{$t[0]}=1;

$FC{$core}=$t[9];
$q_value{$core}=$t[12];
$sig{$core}=$t[13];

}

}
close I;

}#while (my $name=readdir DH) {	

closedir DH;


foreach my $key (keys %gene) {

push(@gene_list,$key);

}


open Profile,">$opts{o}";

for (my $i=0;$i<@file_list;$i++) {

my @tt=split(/_gene_exp.diff/,$file_list[$i]);


my $test_id2=$tt[0].":test_id";
my $gene_id2=$tt[0].":gene_id";
my $gene2=$tt[0].":gene";
my $locus2=$tt[0].":locus";
my $sample_12=$tt[0].":sample_1";
my $sample_22=$tt[0].":sample_2";
my $status2=$tt[0].":status";
my $value_12=$tt[0].":value_1";
my $value_22=$tt[0].":value_2";
my $log2_fold_change2=$tt[0].":log2_fold_change";
my $test_stat2=$tt[0].":test_stat";
my $p_value2=$tt[0].":p_value";
my $q_value2=$tt[0].":q_value";
my $significant2=$tt[0].":significant";


print Profile "$test_id2\t$log2_fold_change2\t$q_value2\t$significant2\t";

}

print Profile "\n";


for (my $i=0;$i<@gene_list;$i++) {

for (my $j=0;$j<@file_list;$j++) {

my $core_1=$file_list[$j].":".$gene_list[$i];

print Profile "$gene_list[$i]\t$FC{$core_1}\t$q_value{$core_1}\t$sig{$core_1}\t";

}

print Profile "\n";

}

close Profile;

###############################################################################################

