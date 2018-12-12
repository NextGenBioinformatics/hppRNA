#!/usr/bin/perl -w
# Program name:		combine_RSEM_gene_results.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-17
# Last update:		2016-06-17



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

##############################################################################
###############################begin to go through each file################

opendir DH, $dir_to_process or die "can't open $dir_to_process:$!";

while (my $name=readdir DH) {

next if $name=~/^\./;

my $full_name ="$dir_to_process/$name";

# read file.

#gene_id	transcript_id(s)	length	effective_length	expected_count	TPM	FPKM
#A1BG	NM_130786	1766	1571.34	18	0.92	0.42


#0:gene_id


#6:FPKM


push(@file_list,$name);


open I,"$full_name";

#############################################################################################

while (<I>)
{

chomp;

my @t=split(/\t/,$_);

if ($t[0] eq "gene_id") {

my $gene_id_title=$t[0];
my $FPKM_title=$t[6];

}

else {

my $gene_id_content=$t[0];
my $FPKM_content=$t[6];



my $core=$name.":".$t[0];

if (!(defined $FPKM{$core})) {$FPKM{$core}=$t[6];}

else {$FPKM{$core}=$FPKM{$core}+$t[6];}

$gene{$t[0]}=1;

}

}
close I;

}#while (my $name=readdir DH) {	

closedir DH;


foreach my $key (keys %gene) {

push(@gene_list,$key);

}



open Profile,">$opts{o}";

print Profile "gene_id\t";

for (my $i=0;$i<@file_list-1;$i++) {

my @tt=split(/.genes.results/,$file_list[$i]);



print Profile "$tt[0]\t";

}

#R_brain_3b.genes.results

my @ttt=split(/.genes.results/,$file_list[(@file_list-1)]);

print Profile "$ttt[0]\n";


for (my $i=0;$i<@gene_list;$i++) {

my @tttt=split(/\./,$gene_list[$i]);

print Profile "$gene_list[$i]\t";

for (my $j=0;$j<@file_list-1;$j++) {

my $core_1=$file_list[$j].":".$gene_list[$i];

print Profile "$FPKM{$core_1}\t";

}

my $core_11=$file_list[(@file_list-1)].":".$gene_list[$i];

print Profile "$FPKM{$core_11}\n";

}

close Profile;

###############################################################################################


