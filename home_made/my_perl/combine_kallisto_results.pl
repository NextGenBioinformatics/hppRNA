#!/usr/bin/perl -w
# Program name:		combine_kallisto_results.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-01
# Last update:		2016-06-01



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
	-j  :   hg19.transcript_to_gene.txt
	-o	:	raw_result


USAGE

die $usage unless $opts{i};
die $usage unless $opts{j};
die $usage unless $opts{o};

my $dir_to_process="$opts{i}";

my @file_list=();
my @gene_list=();

my %tpm=();
my %gene=();

##############################################################################
###############################begin to go through each file################

opendir DH, $dir_to_process or die "can't open $dir_to_process:$!";

while (my $name=readdir DH) {

next if $name=~/^\./;

my $full_name ="$dir_to_process/$name";

# read file.

#target_id	length	eff_length	est_counts	tpm
#NM_032291.3	10951	10764.4	1003.43	6.67146

push(@file_list,$name);


open I,"$full_name";

#############################################################################################

while (<I>)
{

chomp;

my @t=split(/\t/,$_);

if ($t[0] eq "target_id") {

my $target_id_title=$t[0];
my $length_title=$t[1];
my $eff_length_title=$t[2];
my $est_counts_title=$t[3];
my $tpm_title=$t[4];

}

else {

my $target_id_content=$t[0];
my $length_content=$t[1];
my $eff_length_content=$t[2];
my $est_counts_content=$t[3];
my $tpm_content=$t[4];



my $core=$name.":".$t[0];

if (!(defined $tpm{$core})) {$tpm{$core}=$t[4];}

else {$tpm{$core}=$tpm{$core}+$t[4];}

$gene{$t[0]}=1;
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






open Profile,">$opts{o}";

print Profile "transcript_id\tgene_id\t";

for (my $i=0;$i<@file_list-1;$i++) {

my @tt=split(/.abundance.tsv/,$file_list[$i]);



print Profile "$tt[0]\t";

}

my @ttt=split(/.abundance.tsv/,$file_list[(@file_list-1)]);

print Profile "$ttt[0]\n";


for (my $i=0;$i<@gene_list;$i++) {

#my @tttt=split(/\./,$gene_list[$i]);

#print Profile "$gene_list[$i]\t$gene_definition{$tttt[0]}\t";


print Profile "$gene_list[$i]\t$gene_definition{$gene_list[$i]}\t";


for (my $j=0;$j<@file_list-1;$j++) {

my $core_1=$file_list[$j].":".$gene_list[$i];

print Profile "$tpm{$core_1}\t";

}

my $core_11=$file_list[(@file_list-1)].":".$gene_list[$i];

print Profile "$tpm{$core_11}\n";

}

close Profile;

###############################################################################################


