#!/usr/bin/perl -w
# Program name:		combine_eXpress_results.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-24
# Last update:		2016-06-24



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
	-i	:	folder
	-j  :   transcript_to_gene.txt
	-o	:	raw_result


USAGE

die $usage unless $opts{i};
die $usage unless $opts{j};
die $usage unless $opts{o};

my $dir_to_process="$opts{i}";

my @file_list=();
my @gene_list=();

my %fpkm=();
my %gene=();

##############################################################################
###############################begin to go through each file################

opendir DH, $dir_to_process or die "can't open $dir_to_process:$!";

while (my $name=readdir DH) {

next if $name=~/^\./;

my $full_name ="$dir_to_process/$name";

# read file.

#1: target_id

#10: fpkm



#bundle_id	target_id	length	eff_length	tot_counts	uniq_counts	est_counts	eff_counts	ambig_distr_alpha	ambig_distr_beta	fpkm	fpkm_conf_low	fpkm_conf_high	solvable	tpm
#1	NM_014619	5715	5557.795832	1062	4	577.250745	593.578481	1.950081e+01	1.649017e+01	5.370447e+00	4.208136e+00	6.532759e+00	T	1.378008e+01



push(@file_list,$name);


open I,"$full_name";

#############################################################################################

while (<I>)
{

chomp;

my @t=split(/\t/,$_);

if ($t[0] eq "bundle_id") {

my $target_id_title=$t[1];
my $fpkm_title=$t[10];

}

else {

my $target_id_content=$t[1];
my $fpkm_content=$t[10];



my $core=$name.":".$t[1];

if (!(defined $fpkm{$core})) {$fpkm{$core}=$t[10];}

else {$fpkm{$core}=$fpkm{$core}+$t[10];}

$gene{$t[1]}=1;
}

}
close I;

}#while (my $name=readdir DH) {	

closedir DH;


foreach my $key (keys %gene) {

push(@gene_list,$key);

}



#NM_032291       SGIP1
#NM_001308203    SGIP1
#NM_001080397    SLC45A1


my %gene_definition=();


open J,"$opts{j}";

while (<J>)
{

chomp;

my @t=split(/\t/,$_);

$gene_definition{$t[0]}=$t[1];

}
close J;



#R_testis_7c.results.xprs



open Profile,">$opts{o}";

print Profile "transcript_id\tgene_id\t";

for (my $i=0;$i<@file_list-1;$i++) {

my @tt=split(/.results.xprs/,$file_list[$i]);



print Profile "$tt[0]\t";

}

my @ttt=split(/.results.xprs/,$file_list[(@file_list-1)]);

print Profile "$ttt[0]\n";


for (my $i=0;$i<@gene_list;$i++) {

#my @tttt=split(/\./,$gene_list[$i]);

#print Profile "$gene_list[$i]\t$gene_definition{$tttt[0]}\t";


print Profile "$gene_list[$i]\t$gene_definition{$gene_list[$i]}\t";


for (my $j=0;$j<@file_list-1;$j++) {

my $core_1=$file_list[$j].":".$gene_list[$i];

print Profile "$fpkm{$core_1}\t";

}

my $core_11=$file_list[(@file_list-1)].":".$gene_list[$i];

print Profile "$fpkm{$core_11}\n";

}

close Profile;

###############################################################################################


