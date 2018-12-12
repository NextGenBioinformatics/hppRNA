#!/usr/bin/perl -w
# Program name:		combine_cufflinks_isoform_results_new.pl
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

my %FPKM=();
my %gene=();
my %gene_symbol=();

##############################################################################
###############################begin to go through each file################

opendir DH, $dir_to_process or die "can't open $dir_to_process:$!";

while (my $name=readdir DH) {

next if $name=~/^\./;

my $full_name ="$dir_to_process/$name";

# read file.

#tracking_id	class_code	nearest_ref_id	gene_id	gene_short_name	tss_id	locus	length	coverage	FPKM	FPKM_conf_lo	FPKM_conf_hi	FPKM_status
#OR4F5	-	-	OR4F5	OR4F5	TSS15648	chr1:69090-70008	-	-	0	0	0	OK

push(@file_list,$name);


open I,"$full_name";

#############################################################################################

while (<I>)
{

chomp;

my @t=split(/\t/,$_);

if ($t[0] eq "tracking_id") {

my $tracking_id_title=$t[0];
my $class_code_title=$t[1];
my $nearest_ref_id_title=$t[2];
my $gene_id_title=$t[3];
my $gene_short_name_title=$t[4];
my $tss_id_title=$t[5];
my $locus_title=$t[6];
my $len_title=$t[7];
my $coverage_title=$t[8];
my $FPKM_title=$t[9];
my $FPKM_conf_lo_title=$t[10];
my $FPKM_conf_hi_title=$t[11];
my $FPKM_status_title=$t[12];

}

else {

my $tracking_id=$t[0];
my $class_code=$t[1];
my $nearest_ref_id=$t[2];
my $gene_id=$t[3];
my $gene_short_name=$t[4];
my $tss_id=$t[5];
my $locus=$t[6];
my $len=$t[7];
my $coverage=$t[8];
my $FPKM=$t[9];
my $FPKM_conf_lo=$t[10];
my $FPKM_conf_hi=$t[11];
my $FPKM_status=$t[12];

my $core=$name.":".$tracking_id;

if (!(defined $FPKM{$core})) {$FPKM{$core}=$t[9];}

else {$FPKM{$core}=$FPKM{$core}+$t[9];}

$gene{$t[0]}=1;

$gene_symbol{$core}=$t[3];

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

my @tt=split(/_genes.fpkm_tracking/,$file_list[$i]);


my $gene_symbol_id2=$tt[0].":gene";
my $tracking_id2=$tt[0].":isoform";
my $class_code2=$tt[0].":class_code";
my $nearest_ref_id2=$tt[0].":nearest_ref_id";
my $gene_id2=$tt[0].":gene_id";
my $gene_short_name2=$tt[0].":gene_short_name";
my $tss_id2=$tt[0].":tss_id";
my $locus2=$tt[0].":locus";
my $len2=$tt[0].":length";
my $coverage2=$tt[0].":coverage";
my $FPKM2=$tt[0].":FPKM";
my $FPKM_conf_lo2=$tt[0].":FPKM_conf_lo";
my $FPKM_conf_hi2=$tt[0].":FPKM_conf_hi";
my $FPKM_status2=$tt[0].":FPKM_status";

print Profile "$tracking_id2\t$gene_symbol_id2\t$FPKM2\t";

}

print Profile "\n";


for (my $i=0;$i<@gene_list;$i++) {

for (my $j=0;$j<@file_list;$j++) {

my $core_1=$file_list[$j].":".$gene_list[$i];

print Profile "$gene_list[$i]\t$gene_symbol{$core_1}\t$FPKM{$core_1}\t";

}

print Profile "\n";

}

close Profile;

###############################################################################################

