#!/usr/bin/perl -w
# Program name:		generate_Subread_Rscript.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-06
# Last update:		2016-06-06



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"t:s","b1:s","b2:s","d:s","b:s","r:s","o:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-t	:	Paired-End or Single-End
	-h	:	help and usage
	-b1	:	bam1_array
	-b2	:	bam2_array
	-d	:	folder
	-b	:	basename
	-r	:	reference
	-o	:	R_script


USAGE

die $usage unless $opts{t};

die $usage unless $opts{b1};
die $usage unless $opts{b2};

die $usage unless $opts{d};

die $usage unless $opts{b};
die $usage unless $opts{r};

die $usage unless $opts{o};



my $bam_array="";



my @bam_1=split(/\;/,$opts{b1});
my @bam_2=split(/\;/,$opts{b2});



my $basename_real=$opts{b};
my $reference_real=$opts{r};



#"R_BM87_IL7R.sorted.unique.bam",



for (my $i=0;$i<@bam_1;$i++) {$bam_array=$bam_array."\"R_".$bam_1[$i].".sorted.unique.bam\",";}



for (my $i=0;$i<@bam_2-1;$i++) {$bam_array=$bam_array."\"R_".$bam_2[$i].".sorted.unique.bam\",";}



$bam_array=$bam_array."\"R_".$bam_2[@bam_2-1].".sorted.unique.bam\"";



open Profile,">$opts{o}";



print Profile "library(Rsubread)\n\n\n";



print Profile "rm(list=ls(all=TRUE));\n\n\n";



########################################
######################subread starts####
########################################



print Profile "buildindex(basename=\"".$opts{d}."/Subread_featureCounts_DESeq2/Code/".$basename_real."\",reference=\"".$reference_real."\")\n\n\n";





if ($opts{t} eq "Paired-End") {



my $readfile1_array="";

for (my $i=0;$i<@bam_1;$i++) {$readfile1_array=$readfile1_array."\"$opts{d}/Processed_FASTQ/R_".$bam_1[$i]."_1.fastq\",";}

for (my $i=0;$i<@bam_2-1;$i++) {$readfile1_array=$readfile1_array."\"$opts{d}/Processed_FASTQ/R_".$bam_2[$i]."_1.fastq\",";}

$readfile1_array=$readfile1_array."\"$opts{d}/Processed_FASTQ/R_".$bam_2[(@bam_2-1)]."_1.fastq\"";



my $readfile2_array="";

for (my $i=0;$i<@bam_1;$i++) {$readfile2_array=$readfile2_array."\"$opts{d}/Processed_FASTQ/R_".$bam_1[$i]."_2.fastq\",";}

for (my $i=0;$i<@bam_2-1;$i++) {$readfile2_array=$readfile2_array."\"$opts{d}/Processed_FASTQ/R_".$bam_2[$i]."_2.fastq\",";}

$readfile2_array=$readfile2_array."\"$opts{d}/Processed_FASTQ/R_".$bam_2[(@bam_2-1)]."_2.fastq\"";



my $primary_bam_array="";

for (my $i=0;$i<@bam_1;$i++) {$primary_bam_array=$primary_bam_array."\"$opts{d}/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_1[$i].".unique.bam\",";}

for (my $i=0;$i<@bam_2-1;$i++) {$primary_bam_array=$primary_bam_array."\"$opts{d}/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_2[$i].".unique.bam\",";}

$primary_bam_array=$primary_bam_array."\"$opts{d}/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_2[(@bam_2-1)].".unique.bam\"";



print Profile "align(index=\"".$opts{d}."/Subread_featureCounts_DESeq2/Code/".$basename_real."\",readfile1=c(".$readfile1_array."),readfile2=c(".$readfile2_array."),input_format=\"gzFASTQ\",output_format=\"BAM\",output_file=c(".$primary_bam_array."),unique=TRUE,indels=5,PE_orientation=\"fr\")\n\n\n";

}

else {



my $readfile1_array="";

for (my $i=0;$i<@bam_1;$i++) {$readfile1_array=$readfile1_array."\"$opts{d}/Processed_FASTQ/R_".$bam_1[$i].".adapt.qc.fastq\",";}

for (my $i=0;$i<@bam_2-1;$i++) {$readfile1_array=$readfile1_array."\"$opts{d}/Processed_FASTQ/R_".$bam_2[$i].".adapt.qc.fastq\",";}

$readfile1_array=$readfile1_array."\"$opts{d}/Processed_FASTQ/R_".$bam_2[(@bam_2-1)].".adapt.qc.fastq\"";



my $primary_bam_array="";

for (my $i=0;$i<@bam_1;$i++) {$primary_bam_array=$primary_bam_array."\"$opts{d}/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_1[$i].".unique.bam\",";}

for (my $i=0;$i<@bam_2-1;$i++) {$primary_bam_array=$primary_bam_array."\"$opts{d}/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_2[$i].".unique.bam\",";}

$primary_bam_array=$primary_bam_array."\"$opts{d}/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_2[(@bam_2-1)].".unique.bam\"";



print Profile "align(index=\"$basename_real\",readfile1=c(".$readfile1_array."),input_format=\"gzFASTQ\",output_format=\"BAM\",output_file=c(".$primary_bam_array."),unique=TRUE,indels=5)\n\n\n";

}



print Profile "library(Rsamtools)\n\n\n";



for (my $i=0;$i<@bam_1;$i++) {

print Profile "sortBam(\"".$opts{d}."/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_1[$i].".unique.bam\",\"".$opts{d}."/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_1[$i].".sorted.unique.bam\")\n\n\n";
print Profile "indexBam(\"".$opts{d}."/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_1[$i].".sorted.unique.bam\")\n\n\n";

}



for (my $i=0;$i<@bam_2;$i++) {

print Profile "sortBam(\"".$opts{d}."/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_2[$i].".unique.bam\",\"".$opts{d}."/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_2[$i].".sorted.unique.bam\")\n\n\n";
print Profile "indexBam(\"".$opts{d}."/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_2[$i].".sorted.unique.bam\")\n\n\n";

}



###########################################################################################
close Profile;


