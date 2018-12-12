#!/usr/bin/perl -w
# Program name:		generate_featureCounts_DESeq2_Rscript.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-06-06
# Last update:		2016-06-06



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"t:s","c1:s","c2:s","b1:s","b2:s","f:s","d:s","b:s","r:s","p:s","o:s");
my $usage= <<"USAGE";
	Program : $0
	Usage   : $0 [options]
	-t	:	Paired-End or Single-End
	-h	:	help and usage
	-c1	:	condition1
	-c2	:	condition2
	-b1	:	bam1_array
	-b2	:	bam2_array
	-f	:	featurecount_gtf_file
	-d	:	folder
	-b	:	basename
	-r	:	reference
	-p	:	project
	-o	:	R_script


USAGE

die $usage unless $opts{t};

die $usage unless $opts{c1};
die $usage unless $opts{c2};

die $usage unless $opts{b1};
die $usage unless $opts{b2};

die $usage unless $opts{f};

die $usage unless $opts{d};

die $usage unless $opts{b};
die $usage unless $opts{r};
die $usage unless $opts{p};

die $usage unless $opts{o};



my $bam_array="";



my $file_name_prefix="R_".$opts{p}.".Subread_featureCounts_DESeq2";

my $condition_1=$opts{c1};
my $condition_2=$opts{c2};



my @bam_1=split(/\;/,$opts{b1});
my @bam_2=split(/\;/,$opts{b2});



my $basename_real=$opts{b};
my $reference_real=$opts{r};

my $featurecount_gtf_file=$opts{f};



#"R_BM87_IL7R.sorted.unique.bam",



for (my $i=0;$i<@bam_1;$i++) {$bam_array=$bam_array."\"".$opts{d}."/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_1[$i].".sorted.unique.bam\",";}



for (my $i=0;$i<@bam_2-1;$i++) {$bam_array=$bam_array."\"".$opts{d}."/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_2[$i].".sorted.unique.bam\",";}



$bam_array=$bam_array."\"".$opts{d}."/Subread_featureCounts_DESeq2/Processed_BAM/R_".$bam_2[@bam_2-1].".sorted.unique.bam\"";



open Profile,">$opts{o}";



print Profile "library(Rsubread)\n";
print Profile "library(limma)\n";
print Profile "library(edgeR)\n";



print Profile "rm(list=ls(all=TRUE));\n";



########################################
######################subread starts####
########################################



my $isPairedEnd_value="";



if ($opts{t} eq "Paired-End") {$isPairedEnd_value="TRUE";}

else {$isPairedEnd_value="FALSE";}



print Profile "fc<-featureCounts(files=c(".$bam_array."),annot.ext=\"".$featurecount_gtf_file."\",isGTFAnnotationFile=TRUE,GTF.featureType=\"exon\",isPairedEnd=$isPairedEnd_value,nthreads=8,GTF.attrType=\"gene_id\",useMetaFeatures=TRUE)\n\n\n";



print Profile "x <- DGEList(counts=fc\$counts, genes=fc\$annotation)\n";
print Profile "x_rpkm <- rpkm(x,x\$genes\$Length,log=FALSE)\n";



print Profile "write.csv(as.matrix(x),file=\"".$opts{d}."/Subread_featureCounts_DESeq2/DEG/".$file_name_prefix.".".$condition_1."_vs_".$condition_2."_counts.csv\")\n";
print Profile "write.csv(as.matrix(x_rpkm),file=\"".$opts{d}."/Subread_featureCounts_DESeq2/DEG/".$file_name_prefix.".".$condition_1."_vs_".$condition_2."_rpkm.csv\")\n";
print Profile "write.table(as.matrix(x_rpkm),file=\"".$opts{d}."/Subread_featureCounts_DESeq2/DEG/".$file_name_prefix.".".$condition_1."_vs_".$condition_2."_rpkm.txt\",sep=\"\t\",quote = FALSE,row.names = FALSE)\n";



print Profile "library(DESeq2)\n";



print Profile "rm(list=ls(all=TRUE));\n";



print Profile "countData<-read.csv(file=\"".$opts{d}."/Subread_featureCounts_DESeq2/DEG/".$file_name_prefix.".".$condition_1."_vs_".$condition_2."_counts.csv\",sep=\",\",header=TRUE,row.names=1)\n";
print Profile "colData<-read.csv(file=\"".$opts{d}."/Subread_featureCounts_DESeq2/DEG/".$file_name_prefix.".".$condition_1."_vs_".$condition_2."_col.csv\",sep=\",\",header=TRUE,row.names=1)\n";

print Profile "dds<-DESeqDataSetFromMatrix(countData= countData,colData= colData,design=~condition)\n";



print Profile "dds\$condition<-relevel(dds\$condition, \"".$condition_1."\")\n";



print Profile "dds<-DESeq(dds)\n";
print Profile "res<-results(dds)\n";
print Profile "resOrdered<-res[order(res\$padj),]\n";
print Profile "write.csv(as.data.frame(resOrdered),file=\"".$opts{d}."/Subread_featureCounts_DESeq2/DEG/".$file_name_prefix.".".$condition_1."_vs_".$condition_2.".DESeq2.csv\")\n";
print Profile "write.table(as.data.frame(resOrdered),file=\"".$opts{d}."/Subread_featureCounts_DESeq2/DEG/".$file_name_prefix.".".$condition_1."_vs_".$condition_2.".DESeq2.txt\",sep=\"\t\",quote = FALSE,row.names = FALSE)\n";

print Profile "head(resOrdered)\n";



###########################################################################################
close Profile;


