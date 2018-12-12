#!/usr/bin/perl -w
# Program name:		generate_path.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-08-01
# Last update:		2016-08-01



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
	-i	:	The location of software installed
	-o	:	The paths to executable files


USAGE

die $usage unless $opts{i};
die $usage unless $opts{o};



open Profile,">$opts{o}";



{



print Profile "export PATH=$opts{i}/snakemake/Miniconda3/envs/snakemake-tutorial/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/perl-5.22.2/perl-5.22.2/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/jdk1.8.0_91/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/Python-2.7.11:\$PATH\n\n";
print Profile "export PATH=$opts{i}/R-3.3.0/R-3.3.0/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/FastQC:\$PATH\n\n";
print Profile "export PATH=$opts{i}/fastx_toolkit_0.0.13:\$PATH\n\n";
print Profile "export PATH=$opts{i}/Python-2.7.11/Python-2.7.11/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/bowtie2-2.2.9:\$PATH\n\n";
print Profile "export PATH=$opts{i}/tophat-2.1.1.Linux_x86_64:\$PATH\n\n";
print Profile "export PATH=$opts{i}/cufflinks-2.2.1.Linux_x86_64:\$PATH\n\n";
print Profile "export PATH=$opts{i}/samtools-1.3.1/samtools-1.3.1/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/ngsplot-develop/bin:\$PATH\n\n";
print Profile "export NGSPLOT=$opts{i}/ngsplot-develop;\n\n";
print Profile "export PATH=$opts{i}/hisat2-2.0.3-beta:\$PATH\n\n";
print Profile "export PATH=$opts{i}/stringtie-1.2.2.Linux_x86_64:\$PATH\n\n";
print Profile "export PATH=$opts{i}/kallisto_linux-v0.42.5:\$PATH\n\n";
print Profile "export PATH=$opts{i}/STAR-2.5.1b/bin/Linux_x86_64_static:\$PATH\n\n";
print Profile "export PATH=$opts{i}/RSEM-1.2.30:\$PATH\n\n";
print Profile "export PATH=$opts{i}/RSEM-1.2.30/EBSeq:\$PATH\n\n";
print Profile "export PATH=$opts{i}/express-1.5.1-linux_x86_64:\$PATH\n\n";
print Profile "export PATH=$opts{i}/FusionCatcher/fusioncatcher/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/bowtie-1.1.2:\$PATH\n\n";
print Profile "export PATH=$opts{i}/iSeeRNA-1.2.2:\$PATH\n\n";



print Profile "\n\n\n";



}



close Profile;



