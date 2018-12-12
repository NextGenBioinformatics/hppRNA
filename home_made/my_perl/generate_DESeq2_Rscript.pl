#!/usr/bin/perl -w
# Program name:		generate_DESeq2_Rscript.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-03-30
# Last update:		2016-03-30



use warnings;
use strict;
use Getopt::Long;
use List::Util qw(shuffle);

my %opts;
GetOptions(\%opts,"t:s","c1:s","c2:s","b1:s","b2:s","f:s","d:s","o:s");
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
	-o	:	R_script


USAGE

die $usage unless $opts{t};

die $usage unless $opts{c1};
die $usage unless $opts{c2};

die $usage unless $opts{b1};
die $usage unless $opts{b2};

die $usage unless $opts{f};

die $usage unless $opts{d};

die $usage unless $opts{o};



my $bam_array="";



my $condition_1=$opts{c1};
my $condition_2=$opts{c2};



my @bam_1=split(/\;/,$opts{b1});
my @bam_2=split(/\;/,$opts{b2});



my $featurecount_gtf_file=$opts{f};



#"R_BM87_IL7R.sorted.unique.bam",



for (my $i=0;$i<@bam_1;$i++) {$bam_array=$bam_array."\"R_".$bam_1[$i].".sorted.unique.bam\",";}



for (my $i=0;$i<@bam_2-1;$i++) {$bam_array=$bam_array."\"R_".$bam_2[$i].".sorted.unique.bam\",";}



$bam_array=$bam_array."\"R_".$bam_2[@bam_2-1].".sorted.unique.bam\"";



open Profile,">$opts{o}";

print Profile '



library(Rsubread)
library(limma)
library(edgeR)



rm(list=ls(all=TRUE));

setwd("';


print Profile "$opts{d}";


print Profile '");



fc<-featureCounts(files=c(';



print Profile "$bam_array";



print Profile '),annot.ext="';


print Profile "$featurecount_gtf_file";



print Profile '",isGTFAnnotationFile=TRUE,GTF.featureType="exon",isPairedEnd=';

if ($opts{t} eq "Paired-End") {print Profile "TRUE";}

else {print Profile "FALSE";}

print Profile ',nthreads=8,GTF.attrType="gene_id",useMetaFeatures=TRUE)



x <- DGEList(counts=fc$counts, genes=fc$annotation)
x_rpkm <- rpkm(x,x$genes$Length,log=FALSE)



write.csv(as.matrix(x),file=paste("';


print Profile "$condition_1";


print Profile '_vs_';


print Profile "$condition_2";



print Profile '_counts.csv",sep=""))
write.csv(as.matrix(x_rpkm),file=paste("';



print Profile "$condition_1";



print Profile '_vs_';



print Profile "$condition_2";



print Profile '_rpkm.csv",sep=""))



library(DESeq2)



rm(list=ls(all=TRUE));

setwd("';


print Profile "$opts{d}";



print Profile '");

countData<-read.csv(file="';



print Profile "$condition_1";


print Profile '_vs_';



print Profile "$condition_2";



print Profile '_counts.csv",sep=",",header=TRUE,row.names=1)
colData<-read.csv(file="';



print Profile "$condition_1";



print Profile '_vs_';



print Profile "$condition_2";



print Profile '_col.csv",sep=",",header=TRUE,row.names=1)

dds<-DESeqDataSetFromMatrix(countData= countData,colData= colData,design=~condition)




dds$condition<-relevel(dds$condition, "';


print Profile "$condition_1";


print Profile '")




dds<-DESeq(dds)
res<-results(dds)
resOrdered<-res[order(res$padj),]
write.csv(as.data.frame(resOrdered),file="';


print Profile "$condition_1";



print Profile '_vs_';



print Profile "$condition_2";


print Profile '.DESeq2.csv")

head(resOrdered)

';

###########################################################################################
close Profile;


