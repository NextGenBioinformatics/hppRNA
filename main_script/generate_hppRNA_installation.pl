#!/usr/bin/perl -w
# Program name:		generate_hppRNA_installation.pl
# Programmer:		Dapeng Wang
# Email:			wangdp123@gmail.com
# Date:				2016-08-22
# Last update:		2016-08-22



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
	-o	:	The shell for installation of hppRNA pipeline


USAGE

die $usage unless $opts{i};
die $usage unless $opts{o};



open Profile,">$opts{o}";



{



#########################################################################################################################
#########################################################################################################################
#########################################################################################################################



###hppRNA software list:

###Fundemental dependancies

###Perl: https://www.perl.org/
###Java SE Development Kit: http://www.oracle.com/technetwork/java/javase/downloads/index.html
###Python: https://www.python.org/downloads/
###R: https://www.r-project.org/



###Pre-mapping

###FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
###FASTX-Toolkit: http://hannonlab.cshl.edu/fastx_toolkit/
###PRINSEQ-lite: https://sourceforge.net/projects/prinseq/files/standalone/
###Cutadapt: https://pypi.python.org/pypi/cutadapt



###Mapping+Quantification+DEGs test

###Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
###Tophat2: https://ccb.jhu.edu/software/tophat/index.shtml
###Cufflink: http://cole-trapnell-lab.github.io/cufflinks/install/
###SAMtools: https://sourceforge.net/projects/samtools/files/samtools/
###BlackOPs: https://sourceforge.net/projects/rnaseqvariantbl/
###Rsubread (R packages): http://bioconductor.org/packages/release/bioc/html/Rsubread.html
###DESeq2(R packages): https://bioconductor.org/packages/release/bioc/html/DESeq2.html
###STAR 2.5: https://github.com/alexdobin/STAR
###RSEM: https://github.com/bli25ucb/RSEM_tutorial
###EBSeq (R packages): http://bioconductor.org/packages/release/bioc/html/EBSeq.html
###Bowtie: http://bowtie-bio.sourceforge.net/index.shtml
###eXpress: http://bio.math.berkeley.edu/eXpress/index.html
###edgeR (R packages): https://bioconductor.org/packages/release/bioc/html/edgeR.html
###kallisto: https://pachterlab.github.io/kallisto/
###sleuth: https://github.com/pachterlab/sleuth
###HISAT: http://www.ccb.jhu.edu/software/hisat/index.shtml
###StringTie: http://www.ccb.jhu.edu/software/stringtie/index.shtml
###Ballgown (R packages): https://www.bioconductor.org/packages/release/bioc/html/ballgown.html



###Post-mapping

###ngs.plot.r: https://github.com/shenlab-sinai/ngsplot
###wigToBigWig: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
###preprocessCore (R packages): http://bioconductor.org/packages/release/bioc/html/preprocessCore.html
###ggfortify (R packages): https://cran.r-project.org/web/packages/ggfortify/index.html
###cluster (R packages): https://cran.r-project.org/web/packages/cluster/index.html
###gplots (R packages): https://cran.r-project.org/web/packages/gplots/index.html



###DNA variation

###GATK: https://www.broadinstitute.org/gatk/index.php
###FusionCatcher: https://github.com/ndaniel/fusioncatcher

#$opts{i}=$opts{i}

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

print Profile "mkdir -p $opts{i}\n\n";

###Perl:



#print Profile "wget http://www.cpan.org/src/5.0/perl-5.22.2.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/perl-5.22.2.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/perl-5.22.2.tar.gz\n\n";
print Profile "cd $opts{i}/perl-5.22.2\n\n";
print Profile "sh Configure -de -Dprefix=$opts{i}/perl-5.22.2/perl-5.22.2\n\n";
print Profile "make\n\n";
print Profile "make install\n\n";
print Profile "export PATH=$opts{i}/perl-5.22.2/perl-5.22.2/bin:\$PATH\n\n";



################################################################################################################


################################################################################################################
################################################################################################################



###Java SE Development Kit
#print Profile "wget --no-check-certificate --no-cookies --header \"Cookie: oraclelicense=accept-securebackup-cookie\" http://download.oracle.com/otn-pub/java/jdk/8u91-b14/jdk-8u91-linux-x64.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/jdk-8u91-linux-x64.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/jdk-8u91-linux-x64.tar.gz\n\n";
print Profile "export PATH=$opts{i}/jdk1.8.0_91/bin:\$PATH\n\n";



###Python:


#print Profile "wget https://www.python.org/ftp/python/2.7.11/Python-2.7.11.tgz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/Python-2.7.11.tgz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/Python-2.7.11.tgz\n\n";
print Profile "cd $opts{i}/Python-2.7.11\n\n";
print Profile "./configure --prefix=$opts{i}/Python-2.7.11/Python-2.7.11\n\n";
print Profile "make\n\n";
print Profile "make install\n\n";
print Profile "export PATH=$opts{i}/Python-2.7.11:\$PATH\n\n";
print Profile "export PATH=$opts{i}/Python-2.7.11/Python-2.7.11/bin:\$PATH\n\n";



###setuptools 21.0.0
print Profile "mkdir $opts{i}/setuptools\n\n";
print Profile "wget --no-check-certificate https://bootstrap.pypa.io/ez_setup.py -P $opts{i}/setuptools\n\n";
print Profile "cd $opts{i}/setuptools\n\n";
print Profile "python ez_setup.py --insecure\n\n";



###numpy 1.11.0


#print Profile "wget https://pypi.python.org/packages/1a/5c/57c6920bf4a1b1c11645b625e5483d778cedb3823ba21a017112730f0a12/numpy-1.11.0.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/numpy-1.11.0.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/numpy-1.11.0.tar.gz\n\n";
print Profile "cd $opts{i}/numpy-1.11.0\n\n";
print Profile "python setup.py build -j 4 install\n\n";



###R dependencies:


#print Profile "wget https://www.openssl.org/source/openssl-1.0.2h.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/openssl-1.0.2h.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/openssl-1.0.2h.tar.gz\n\n";
print Profile "cd $opts{i}/openssl-1.0.2h\n\n";
print Profile "./config --prefix=$opts{i}/openssl-1.0.2h/openssl-1.0.2h shared -fPIC\n\n";
print Profile "make depend\n\n";
print Profile "make\n\n";
print Profile "make install\n\n";
print Profile "export PATH=$opts{i}/openssl-1.0.2h/openssl-1.0.2h/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/openssl-1.0.2h:\$PATH\n\n";


#print Profile "wget http://zlib.net/zlib-1.2.8.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/zlib-1.2.8.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/zlib-1.2.8.tar.gz\n\n";
print Profile "cd $opts{i}/zlib-1.2.8\n\n";
print Profile "./configure --prefix=$opts{i}/zlib-1.2.8\n\n";
print Profile "make\n\n";
print Profile "make install\n\n";


#print Profile "wget http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/bzip2-1.0.6.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/bzip2-1.0.6.tar.gz\n\n";
print Profile "cd $opts{i}/bzip2-1.0.6\n\n";
print Profile "make -f Makefile-libbz2_so\n\n";
print Profile "make clean\n\n";
print Profile "make\n\n";
print Profile "make install PREFIX=$opts{i}/bzip2-1.0.6\n\n";


#print Profile "wget http://tukaani.org/xz/xz-5.2.2.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/xz-5.2.2.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/xz-5.2.2.tar.gz\n\n";
print Profile "cd $opts{i}/xz-5.2.2\n\n";
print Profile "./configure --prefix=$opts{i}/xz-5.2.2\n\n";
print Profile "make -j3\n\n";
print Profile "make install\n\n";


#print Profile "wget ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.38.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/pcre-8.38.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/pcre-8.38.tar.gz\n\n";
print Profile "cd $opts{i}/pcre-8.38\n\n";
print Profile "./configure --enable-utf8 --prefix=$opts{i}/pcre-8.38\n\n";
print Profile "make -j3\n\n";
print Profile "make install\n\n";


#print Profile "wget ftp://ftp.gnu.org/gnu/libtool/libtool-2.4.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/libtool-2.4.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/libtool-2.4.tar.gz\n\n";
print Profile "cd $opts{i}/libtool-2.4\n\n";
print Profile "./configure --prefix=$opts{i}/libtool-2.4\n\n";
print Profile "make\n\n";
print Profile "make install\n\n";
print Profile "export PATH=$opts{i}/libtool-2.4/bin:\$PATH\n\n";


print Profile "PKG_CONFIG_PATH=$opts{i}/openssl-1.0.2h/openssl-1.0.2h/bin:\$PKG_CONFIG_PATH\n\n";

#print Profile "wget https://curl.haxx.se/download/curl-7.48.0.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/curl-7.48.0.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/curl-7.48.0.tar.gz\n\n";
print Profile "cd $opts{i}/curl-7.48.0\n\n";
print Profile "./configure --prefix=$opts{i}/curl-7.48.0/curl-7.48.0 --with-ssl\n\n";
print Profile "make\n\n";
print Profile "make install\n\n";
print Profile "export PATH=$opts{i}/curl-7.48.0/curl-7.48.0/bin:\$PATH\n\n";



print Profile "export PATH=$opts{i}/bzip2-1.0.6/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/xz-5.2.2/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/zlib-1.2.8:\$PATH\n\n";
print Profile "export PATH=$opts{i}/pcre-8.38/bin:\$PATH\n\n";

print Profile "export LD_LIBRARY_PATH=$opts{i}/curl-7.48.0/curl-7.48.0/lib:\$LD_LIBRARY_PATH\n\n";
print Profile "export CFLAGS=\"-I$opts{i}/zlib-1.2.8/include -I$opts{i}/bzip2-1.0.6/include -I$opts{i}/xz-5.2.2/include -I$opts{i}/pcre-8.38/include -I$opts{i}/libtool-2.4/include -I$opts{i}/openssl-1.0.2h/openssl-1.0.2h/include\"\n\n";
print Profile "export LDFLAGS=\"-L$opts{i}/zlib-1.2.8/lib -L$opts{i}/bzip2-1.0.6/lib -L$opts{i}/xz-5.2.2/lib -L$opts{i}/pcre-8.38/lib -L$opts{i}/libtool-2.4/lib -L$opts{i}/openssl-1.0.2h/openssl-1.0.2h/lib\"\n\n";



###R:
#print Profile "wget https://cloud.r-project.org/src/base/R-3/R-3.3.0.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/R-3.3.0.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/R-3.3.0.tar.gz\n\n";
print Profile "cd $opts{i}/R-3.3.0\n\n";
print Profile "./configure --prefix=$opts{i}/R-3.3.0/R-3.3.0\n\n";
print Profile "make\n\n";
print Profile "make install\n\n";
print Profile "export PATH=$opts{i}/R-3.3.0/R-3.3.0/bin:\$PATH\n\n";



###############################################Till now, all the paramount programs have been installed###############################










#########################################################################################################################################







###Pre-mapping



###FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#print Profile "wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/fastqc_v0.11.5.zip -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "unzip $opts{i}/fastqc_v0.11.5.zip\n\n";
print Profile "cd $opts{i}/FastQC\n\n";
print Profile "chmod 0755 fastqc\n\n";
print Profile "export PATH=$opts{i}/FastQC:\$PATH\n\n";



###FASTX-Toolkit: http://hannonlab.cshl.edu/fastx_toolkit/


#print Profile "wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar jxf $opts{i}/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2\n\n";
print Profile "mv bin fastx_toolkit_0.0.13\n\n";
print Profile "export PATH=$opts{i}/fastx_toolkit_0.0.13:\$PATH\n\n";



###PRINSEQ-lite: https://sourceforge.net/projects/prinseq/files/standalone/
#print Profile "wget http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/prinseq-lite-0.20.4.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/prinseq-lite-0.20.4.tar.gz\n\n";


###Cutadapt: https://pypi.python.org/pypi/cutadapt


#print Profile "wget https://pypi.python.org/packages/b0/40/9c853c1dc32414263031f7f956f859ee7af145c84b8391640a1d79ce25df/cutadapt-1.9.1.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/cutadapt-1.9.1.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/cutadapt-1.9.1.tar.gz\n\n";
print Profile "cd $opts{i}/cutadapt-1.9.1\n\n";
print Profile "python setup.py install\n\n";
print Profile "export PATH=$opts{i}/Python-2.7.11/Python-2.7.11/bin:\$PATH\n\n";


###Mapping+Quantification+DEGs test



###Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml


#print Profile "wget http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/bowtie2-2.2.9-linux-x86_64.zip -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "unzip $opts{i}/bowtie2-2.2.9-linux-x86_64.zip\n\n";
print Profile "export PATH=$opts{i}/bowtie2-2.2.9:\$PATH\n\n";


###Tophat2: https://ccb.jhu.edu/software/tophat/index.shtml


#print Profile "wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/tophat-2.1.1.Linux_x86_64.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/tophat-2.1.1.Linux_x86_64.tar.gz\n\n";
print Profile "export PATH=$opts{i}/tophat-2.1.1.Linux_x86_64:\$PATH\n\n";



###Cufflink: http://cole-trapnell-lab.github.io/cufflinks/install/


#print Profile "wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/cufflinks-2.2.1.Linux_x86_64.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/cufflinks-2.2.1.Linux_x86_64.tar.gz\n\n";
print Profile "export PATH=$opts{i}/cufflinks-2.2.1.Linux_x86_64:\$PATH\n\n";



###SAMtools: https://sourceforge.net/projects/samtools/files/samtools/


#print Profile "wget http://downloads.sourceforge.net/project/samtools/samtools/1.3.1/samtools-1.3.1.tar.bz2 -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/samtools-1.3.1.tar.bz2 -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar jxf $opts{i}/samtools-1.3.1.tar.bz2\n\n";
print Profile "cd $opts{i}/samtools-1.3.1\n\n";
print Profile "./configure --prefix=$opts{i}/samtools-1.3.1/samtools-1.3.1\n\n";
print Profile "make\n\n";
print Profile "make install\n\n";
print Profile "export PATH=$opts{i}/samtools-1.3.1/samtools-1.3.1/bin:\$PATH\n\n";



###BlackOPs: https://sourceforge.net/projects/rnaseqvariantbl/


#print Profile "wget http://downloads.sourceforge.net/project/rnaseqvariantbl/BlackOPs_v1.06.zip -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/BlackOPs_v1.06.zip -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "unzip $opts{i}/BlackOPs_v1.06.zip\n\n";


###############################Installation of R packages start##################################


print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/R_package.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/R_package.tar.gz\n\n";





###Rsubread (R packages): http://bioconductor.org/packages/release/bioc/html/Rsubread.html
#print Profile "wget http://bioconductor.org/packages/3.3/bioc/src/contrib/Rsubread_1.22.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/Rsubread_1.22.0.tar.gz\n\n";


###BiocGenerics
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/BiocGenerics_0.18.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/BiocGenerics_0.18.0.tar.gz\n\n";


###S4Vectors
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/S4Vectors_0.10.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/S4Vectors_0.10.0.tar.gz\n\n";


###IRanges
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/IRanges_2.6.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/IRanges_2.6.0.tar.gz\n\n";


###GenomeInfoDb
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/GenomeInfoDb_1.8.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/GenomeInfoDb_1.8.0.tar.gz\n\n";


###zlibbioc
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/zlibbioc_1.18.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/zlibbioc_1.18.0.tar.gz\n\n";


###XVector
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/XVector_0.12.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/XVector_0.12.0.tar.gz\n\n";


###GenomicRanges
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/GenomicRanges_1.24.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/GenomicRanges_1.24.0.tar.gz\n\n";


###Biobase
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/Biobase_2.32.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/Biobase_2.32.0.tar.gz\n\n";



###SummarizedExperiment
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/SummarizedExperiment_1.2.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/SummarizedExperiment_1.2.0.tar.gz\n\n";



###lambda.r
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/lambda.r/lambda.r_1.1.6.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/lambda.r_1.1.6.tar.gz\n\n";



###futile.options
#print Profile "wget https://cran.rstudio.com/src/contrib/futile.options_1.0.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/futile.options_1.0.0.tar.gz\n\n";



###futile.logger
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/futile.logger/futile.logger_1.4.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/futile.logger_1.4.tar.gz\n\n";



###snow
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/snow/snow_0.3-13.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/snow_0.3-13.tar.gz\n\n";


###BiocParallel
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/BiocParallel_1.6.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/BiocParallel_1.6.0.tar.gz\n\n";


###DBI
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/DBI/DBI_0.3.1.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/DBI_0.3.1.tar.gz\n\n";



###RSQLite
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/RSQLite/RSQLite_0.11.4.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/RSQLite_0.11.4.tar.gz\n\n";

###AnnotationDbi
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/AnnotationDbi_1.34.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/AnnotationDbi_1.34.0.tar.gz\n\n";


###XML
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/XML/XML_3.98-1.3.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/XML_3.98-1.3.tar.gz\n\n";

###xtable
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/xtable/xtable_1.8-0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/xtable_1.8-0.tar.gz\n\n";


###annotate
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/annotate_1.50.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/annotate_1.50.0.tar.gz\n\n";


###genefilter
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/genefilter_1.54.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/genefilter_1.54.0.tar.gz\n\n";


###locfit
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/locfit_1.5-9.tar.gz\n\n";


###RColorBrewer
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/RColorBrewer/RColorBrewer_1.0-5.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/RColorBrewer_1.0-5.tar.gz\n\n";


###geneplotter
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/geneplotter_1.50.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/geneplotter_1.50.0.tar.gz\n\n";



###digest
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/digest/digest_0.6.8.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/digest_0.6.8.tar.gz\n\n";


###gtable
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/gtable/gtable_0.1.2.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/gtable_0.1.2.tar.gz\n\n";


###Rcpp
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_0.12.3.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/Rcpp_0.12.3.tar.gz\n\n";


###plyr
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/plyr/plyr_1.8.2.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/plyr_1.8.2.tar.gz\n\n";







###stringi
#print Profile "wget https://github.com/Rexamine/stringi/archive/master.zip -P $opts{i}/R_package\n\n";
#print Profile "cd $opts{i}/R_package\n\n";
#print Profile "unzip master.zip\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/stringi-master\n\n";
#print Profile "rm $opts{i}/R_package/master.zip\n\n";



###magrittr
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/magrittr/magrittr_1.0.1.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/magrittr_1.0.1.tar.gz\n\n";


###stringr
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/stringr/stringr_0.6.2.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/stringr_0.6.2.tar.gz\n\n";


###reshape2
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/reshape2/reshape2_1.4.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/reshape2_1.4.tar.gz\n\n";







###dichromat
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/dichromat/dichromat_1.2-4.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/dichromat_1.2-4.tar.gz\n\n";


###colorspace
#print Profile "wget https://cran.r-project.org/src/contrib/colorspace_1.2-6.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/colorspace_1.2-6.tar.gz\n\n";


###munsell
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/munsell/munsell_0.4.2.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/munsell_0.4.2.tar.gz\n\n";


###labeling
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/labeling/labeling_0.2.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/labeling_0.2.tar.gz\n\n";



###scales
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/scales/scales_0.3.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/scales_0.3.0.tar.gz\n\n";


###ggplot2
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_2.0.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/ggplot2_2.0.0.tar.gz\n\n";




###Formula
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/Formula/Formula_1.2-0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/Formula_1.2-0.tar.gz\n\n";



###latticeExtra
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/latticeExtra/latticeExtra_0.6-26.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/latticeExtra_0.6-26.tar.gz\n\n";



###acepack
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/acepack/acepack_1.3-3.2.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/acepack_1.3-3.2.tar.gz\n\n";



###gridExtra
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/gridExtra/gridExtra_2.2.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/gridExtra_2.2.0.tar.gz\n\n";



###chron
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/chron/chron_2.3-46.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/chron_2.3-46.tar.gz\n\n";



###data.table
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.9.4.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/data.table_1.9.4.tar.gz\n\n";


###Hmisc
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_3.17-3.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/Hmisc_3.17-3.tar.gz\n\n";



###RcppArmadillo
#print Profile "wget wget https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.4.500.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/RcppArmadillo_0.4.500.0.tar.gz\n\n";


###DESeq2(R packages): https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/DESeq2_1.12.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/DESeq2_1.12.0.tar.gz\n\n";



###STAR 2.5: https://github.com/alexdobin/STAR/releases


#print Profile "wget https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/2.5.1b.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/2.5.1b.tar.gz\n\n";
print Profile "export PATH=$opts{i}/STAR-2.5.1b/bin/Linux_x86_64_static:\$PATH\n\n";


###RSEM: https://github.com/bli25ucb/RSEM_tutorial



###clear the variables before RSEM installation
print Profile "export LD_LIBRARY_PATH=\n\n";
print Profile "export CFLAGS=\n\n";
print Profile "export LDFLAGS=\n\n";



#print Profile "wget https://github.com/deweylab/RSEM/archive/v1.2.30.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/v1.2.30.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/v1.2.30.tar.gz\n\n";
print Profile "cd $opts{i}/RSEM-1.2.30\n\n";
print Profile "make\n\n";
print Profile "make ebseq\n\n";
print Profile "export PATH=$opts{i}/RSEM-1.2.30:\$PATH\n\n";
print Profile "export PATH=$opts{i}/RSEM-1.2.30/EBSeq:\$PATH\n\n";



###blockmodeling
#print Profile "wget https://cran.rstudio.com/src/contrib/Archive/blockmodeling/blockmodeling_0.1.7.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/blockmodeling_0.1.7.tar.gz\n\n";



###gtools
#print Profile "wget https://cran.rstudio.com/src/contrib/Archive/gtools/gtools_3.4.2.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/gtools_3.4.2.tar.gz\n\n";



###gdata
#print Profile "wget https://cran.rstudio.com/src/contrib/Archive/gdata/gdata_2.16.1.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/gdata_2.16.1.tar.gz\n\n";



###bitops
#print Profile "wget https://cran.rstudio.com/src/contrib/Archive/bitops/bitops_1.0-5.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/bitops_1.0-5.tar.gz\n\n";



###caTools
#print Profile "wget https://cran.rstudio.com/src/contrib/Archive/caTools/caTools_1.17.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/caTools_1.17.tar.gz\n\n";



###gplots
#print Profile "wget https://cran.rstudio.com/src/contrib/Archive/gplots/gplots_3.0.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/gplots_3.0.0.tar.gz\n\n";


###memoise
#print Profile "wget https://cran.rstudio.com/src/contrib/Archive/memoise/memoise_0.2.1.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/memoise_0.2.1.tar.gz\n\n";


###crayon
#print Profile "wget https://cran.rstudio.com/src/contrib/Archive/crayon/crayon_1.3.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/crayon_1.3.0.tar.gz\n\n";


###praise
#print Profile "wget https://cran.rstudio.com/src/contrib/praise_1.0.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/praise_1.0.0.tar.gz\n\n";


###R6
#print Profile "wget https://cran.rstudio.com/src/contrib/Archive/R6/R6_2.1.1.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/R6_2.1.1.tar.gz\n\n";



###testthat
#print Profile "wget https://cran.rstudio.com/src/contrib/Archive/testthat/testthat_1.0.1.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/testthat_1.0.1.tar.gz\n\n";



###EBSeq (R packages): http://bioconductor.org/packages/release/bioc/html/EBSeq.html
#print Profile "wget http://bioconductor.org/packages/3.3/bioc/src/contrib/EBSeq_1.12.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/EBSeq_1.12.0.tar.gz\n\n";



###Bowtie: http://bowtie-bio.sourceforge.net/index.shtml


#print Profile "wget http://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/bowtie-1.1.2-linux-x86_64.zip -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "unzip $opts{i}/bowtie-1.1.2-linux-x86_64.zip\n\n";
print Profile "export PATH=$opts{i}/bowtie-1.1.2:\$PATH\n\n";







###eXpress: http://bio.math.berkeley.edu/eXpress/index.html


#print Profile "wget http://bio.math.berkeley.edu/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/express-1.5.1-linux_x86_64.tgz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/express-1.5.1-linux_x86_64.tgz\n\n";
print Profile "export PATH=$opts{i}/express-1.5.1-linux_x86_64:\$PATH\n\n";




###limma
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/limma_3.28.2.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/limma_3.28.2.tar.gz\n\n";



###edgeR (R packages): https://bioconductor.org/packages/release/bioc/html/edgeR.html
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/edgeR_3.14.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/edgeR_3.14.0.tar.gz\n\n";


###kallisto: https://pachterlab.github.io/kallisto/


#print Profile "wget https://github.com/pachterlab/kallisto/releases/download/v0.42.5/kallisto_linux-v0.42.5.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/kallisto_linux-v0.42.5.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/kallisto_linux-v0.42.5.tar.gz\n\n";
print Profile "export PATH=$opts{i}/kallisto_linux-v0.42.5:\$PATH\n\n";


###memoise
#print Profile "wget https://cran.r-project.org/src/contrib/memoise_1.0.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/memoise_1.0.0.tar.gz\n\n";



###git2r
print Profile "R CMD INSTALL $opts{i}/R_package/git2r_0.15.0.tar.gz\n\n";



###jsonlite
print Profile "R CMD INSTALL $opts{i}/R_package/jsonlite_1.0.tar.gz\n\n";



###rstudioapi
print Profile "R CMD INSTALL $opts{i}/R_package/rstudioapi_0.6.tar.gz\n\n";



###whisker
print Profile "R CMD INSTALL $opts{i}/R_package/whisker_0.3-2.tar.gz\n\n";



###withr
print Profile "R CMD INSTALL $opts{i}/R_package/withr_1.0.2.tar.gz\n\n";



###mime
print Profile "R CMD INSTALL $opts{i}/R_package/mime_0.5.tar.gz\n\n";



###curl
print Profile "R CMD INSTALL $opts{i}/R_package/curl_1.2.tar.gz\n\n";



###openssl
print Profile "R CMD INSTALL $opts{i}/R_package/openssl_0.9.4.tar.gz\n\n";



###httr
print Profile "R CMD INSTALL $opts{i}/R_package/httr_1.2.1.tar.gz\n\n";



###devtools
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/devtools/devtools_1.11.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/devtools_1.11.0.tar.gz\n\n";



###################################################################redo the paths to help R to communicate with internet#################################################



print Profile "export PATH=$opts{i}/perl-5.22.2/perl-5.22.2/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/jdk1.8.0_91/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/openssl-1.0.2h/openssl-1.0.2h/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/openssl-1.0.2h:\$PATH\n\n";
print Profile "export PATH=$opts{i}/libtool-2.4/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/curl-7.48.0/curl-7.48.0/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/bzip2-1.0.6/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/xz-5.2.2/bin:\$PATH\n\n";
print Profile "export PATH=$opts{i}/zlib-1.2.8:\$PATH\n\n";
print Profile "export PATH=$opts{i}/pcre-8.38/bin:\$PATH\n\n";
print Profile "export LD_LIBRARY_PATH=$opts{i}/curl-7.48.0/curl-7.48.0/lib:\$LD_LIBRARY_PATH\n\n";
print Profile "export CFLAGS=\"-I$opts{i}/zlib-1.2.8/include -I$opts{i}/bzip2-1.0.6/include -I$opts{i}/xz-5.2.2/include -I$opts{i}/pcre-8.38/include -I$opts{i}/libtool-2.4/include -I$opts{i}/openssl-1.0.2h/openssl-1.0.2h/include\"\n\n";
print Profile "export LDFLAGS=\"-L$opts{i}/zlib-1.2.8/lib -L$opts{i}/bzip2-1.0.6/lib -L$opts{i}/xz-5.2.2/lib -L$opts{i}/pcre-8.38/lib -L$opts{i}/libtool-2.4/lib -L$opts{i}/openssl-1.0.2h/openssl-1.0.2h/lib\"\n\n";
print Profile "export PATH=$opts{i}/Python-2.7.11:\$PATH\n\n";
print Profile "export PATH=$opts{i}/R-3.3.0/R-3.3.0/bin:\$PATH\n\n";



#########################################################################################################################################################################



###sleuth: https://github.com/pachterlab/sleuth
print Profile "mkdir $opts{i}/sleuth\n\n";
print Profile "cd $opts{i}/sleuth\n\n";



print Profile "echo \"source(\\\"http://bioconductor.org/biocLite.R\\\")\" >> sleuth_install.R\n\n";
print Profile "echo \"biocLite(\\\"rhdf5\\\")\" >> sleuth_install.R\n\n";
print Profile "echo \"devtools::install_github(\\\"pachterlab/sleuth\\\")\" >> sleuth_install.R\n\n";

print Profile "R --slave --vanilla < $opts{i}/sleuth/sleuth_install.R\n\n";



###HISAT: http://ccb.jhu.edu/software/hisat2/index.shtml


#print Profile "wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.3-beta-Linux_x86_64.zip -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/hisat2-2.0.3-beta-Linux_x86_64.zip -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "unzip $opts{i}/hisat2-2.0.3-beta-Linux_x86_64.zip\n\n";
print Profile "export PATH=$opts{i}/hisat2-2.0.3-beta:\$PATH\n\n";



###StringTie: http://www.ccb.jhu.edu/software/stringtie/index.shtml


#print Profile "wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.2.2.Linux_x86_64.tar.gz -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/stringtie-1.2.2.Linux_x86_64.tar.gz -P $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/stringtie-1.2.2.Linux_x86_64.tar.gz\n\n";
print Profile "export PATH=$opts{i}/stringtie-1.2.2.Linux_x86_64:\$PATH\n\n";



###sva
#print Profile "wget https://www.bioconductor.org/packages/3.3/bioc/src/contrib/sva_3.20.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/sva_3.20.0.tar.gz\n\n";



###Biostrings
#print Profile "wget https://www.bioconductor.org/packages/3.3/bioc/src/contrib/Biostrings_2.40.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/Biostrings_2.40.0.tar.gz\n\n";



###RCurl
#print Profile "wget https://cran.rstudio.com/src/contrib/Archive/RCurl/RCurl_1.95-4.7.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/RCurl_1.95-4.7.tar.gz\n\n";



###Rsamtools
#print Profile "wget https://www.bioconductor.org/packages/3.3/bioc/src/contrib/Rsamtools_1.24.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/Rsamtools_1.24.0.tar.gz\n\n";



###GenomicAlignments
#print Profile "wget https://www.bioconductor.org/packages/3.3/bioc/src/contrib/GenomicAlignments_1.8.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/GenomicAlignments_1.8.0.tar.gz\n\n";



###rtracklayer
#print Profile "wget https://www.bioconductor.org/packages/3.3/bioc/src/contrib/rtracklayer_1.32.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/rtracklayer_1.32.0.tar.gz\n\n";



###Ballgown (R packages): https://www.bioconductor.org/packages/release/bioc/html/ballgown.html
#print Profile "wget https://www.bioconductor.org/packages/3.3/bioc/src/contrib/ballgown_2.4.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/ballgown_2.4.0.tar.gz\n\n";



###Post-mapping



###iterators
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/iterators/iterators_1.0.7.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/iterators_1.0.7.tar.gz\n\n";


###foreach
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/foreach/foreach_1.4.2.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/foreach_1.4.2.tar.gz\n\n";


###doMC
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/doMC/doMC_1.3.3.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/doMC_1.3.3.tar.gz\n\n";


###caTools
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/caTools/caTools_1.17.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/caTools_1.17.tar.gz\n\n";


###BSgenome
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/BSgenome_1.40.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/BSgenome_1.40.0.tar.gz\n\n";


###hwriter
#print Profile "wget https://cran.rstudio.com/src/contrib/Archive/hwriter/hwriter_1.3.1.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/hwriter_1.3.1.tar.gz\n\n";


###ShortRead
#print Profile "wget https://bioconductor.org/packages/3.3/bioc/src/contrib/ShortRead_1.30.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/ShortRead_1.30.0.tar.gz\n\n";



###ngs.plot.r: https://github.com/shenlab-sinai/ngsplot


#print Profile "wget https://github.com/shenlab-sinai/ngsplot/archive/develop.zip -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/develop.zip -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "unzip $opts{i}/develop.zip\n\n";
print Profile "rm $opts{i}/develop.zip\n\n";
print Profile "export PATH=$opts{i}/ngsplot-develop/bin:\$PATH\n\n";
print Profile "export NGSPLOT=$opts{i}/ngsplot-develop;\n\n";
print Profile "cd $opts{i}/ngsplot-develop\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/ngsplotdb_hg19_75_3.00.tar.gz -P $opts{i}/ngsplot-develop\n\n";
print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/ngsplotdb_mm10_75_3.00.tar.gz -P $opts{i}/ngsplot-develop\n\n";

print Profile "ngsplotdb.py -y install $opts{i}/ngsplot-develop/ngsplotdb_hg19_75_3.00.tar.gz\n\n";
print Profile "ngsplotdb.py -y install $opts{i}/ngsplot-develop/ngsplotdb_mm10_75_3.00.tar.gz\n\n";



###wigToBigWig: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
print Profile "mkdir $opts{i}/wig2bigwig\n\n";
print Profile "wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig -P $opts{i}/wig2bigwig\n\n";
print Profile "cd $opts{i}/wig2bigwig\n\n";
print Profile "chmod 0755 wigToBigWig\n\n";
print Profile "export PATH=$opts{i}/wig2bigwig:\$PATH\n\n";



print Profile "wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes -P $opts{i}/wig2bigwig\n\n";
print Profile "chmod 0755 $opts{i}/wig2bigwig/fetchChromSizes\n\n";
print Profile "$opts{i}/wig2bigwig/fetchChromSizes hg19 > $opts{i}/wig2bigwig/hg19.chrom.sizes\n\n";
print Profile "$opts{i}/wig2bigwig/fetchChromSizes mm10 > $opts{i}/wig2bigwig/mm10.chrom.sizes\n\n";



###preprocessCore (R packages): http://bioconductor.org/packages/release/bioc/html/preprocessCore.html
#print Profile "wget http://bioconductor.org/packages/3.3/bioc/src/contrib/preprocessCore_1.34.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/preprocessCore_1.34.0.tar.gz\n\n";




###proto
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/proto/proto_0.3-9.2.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/proto_0.3-9.2.tar.gz\n\n";



###assertthat
#print Profile "wget https://cran.r-project.org/src/contrib/assertthat_0.1.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/assertthat_0.1.tar.gz\n\n";


###lazyeval
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/lazyeval/lazyeval_0.1.9.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/lazyeval_0.1.9.tar.gz\n\n";


###BH
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/BH/BH_1.60.0-1.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/BH_1.60.0-1.tar.gz\n\n";


###dplyr
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_0.4.2.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/dplyr_0.4.2.tar.gz\n\n";


###tidyr
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/tidyr/tidyr_0.4.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/tidyr_0.4.0.tar.gz\n\n";


###ggfortify (R packages): https://cran.r-project.org/web/packages/ggfortify/index.html
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/ggfortify/ggfortify_0.1.0.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/ggfortify_0.1.0.tar.gz\n\n";



###cluster (R packages): https://cran.r-project.org/web/packages/cluster/index.html
#print Profile "wget https://cran.r-project.org/src/contrib/Archive/cluster/cluster_2.0.3.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/cluster_2.0.3.tar.gz\n\n";





###lncRNA de novo



###iSeeRNA



print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/iSeeRNA-1.2.2.x86_64.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/iSeeRNA-1.2.2.x86_64.tar.gz\n\n";
print Profile "cd $opts{i}/iSeeRNA-1.2.2\n\n";
print Profile "sh auto_download_data.sh hg19\n\n";
print Profile "sh auto_download_data.sh mm10\n\n";



###CircRNA detection



###DCC package



print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/v0.4.4.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/v0.4.4.tar.gz\n\n";
print Profile "cd $opts{i}/DCC-0.4.4\n\n";
print Profile "python setup.py install\n\n";



###CircTest



###aod (dependence)
#print Profile "wget https://cran.r-project.org/src/contrib/aod_1.3.tar.gz -P $opts{i}/R_package\n\n";
print Profile "R CMD INSTALL $opts{i}/R_package/aod_1.3.tar.gz\n\n";



print Profile "mkdir $opts{i}/circTest\n\n";
print Profile "cd $opts{i}/circTest\n\n";



print Profile "echo \"require(devtools)\" >> circTest_install.R\n\n";
print Profile "echo \"install_github(\\\"dieterich-lab/CircTest\\\")\" >> circTest_install.R\n\n";
print Profile "echo \"library(CircTest)\" >> circTest_install.R\n\n";



print Profile "R --slave --vanilla < $opts{i}/circTest/circTest_install.R\n\n";



###DNA variation

###picard


#print Profile "wget https://github.com/broadinstitute/picard/releases/download/2.3.0/picard-tools-2.3.0.zip -P $opts{i}\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/picard-tools-2.3.0.zip -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "unzip $opts{i}/picard-tools-2.3.0.zip\n\n";



###GATK: https://www.broadinstitute.org/gatk/index.php
#print Profile "mkdir $opts{i}/GATK\n\n";
#print Profile "cd $opts{i}/GATK\n\n";
#print Profile "tar jxf $opts{i}/GATK/GenomeAnalysisTK-3.5.tar.bz2\n\n";



############################################################################install numpy required by fusioncatcher#############################

###################################################################################################################################################



###FusionCatcher: https://github.com/ndaniel/fusioncatcher
print Profile "mkdir $opts{i}/FusionCatcher\n\n";
print Profile "cd $opts{i}/FusionCatcher\n\n";
print Profile "wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py && python bootstrap.py -t --download -y\n\n";
print Profile "export PATH=$opts{i}/FusionCatcher/fusioncatcher/bin:\$PATH\n\n";

print Profile "fusioncatcher-build --ftp-ensembl-path=pub/release-80 -g mus_musculus -o $opts{i}/FusionCatcher/fusioncatcher/data/mus_musculus\n\n";
print Profile "fusioncatcher-build -g mus_musculus -o $opts{i}/FusionCatcher/fusioncatcher/data/mus_musculus\n\n";


###Snakemake Tutorial

print Profile "mkdir $opts{i}/snakemake\n\n";

print Profile "cd $opts{i}/snakemake\n\n";



#print Profile "wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/Miniconda3-latest-Linux-x86_64.sh\n\n";


print Profile "bash Miniconda3-latest-Linux-x86_64.sh -b -p $opts{i}/snakemake/Miniconda3\n\n";

print Profile "export PATH=$opts{i}/snakemake/Miniconda3/bin:\$PATH\n\n";

print Profile "mkdir $opts{i}/snakemake/snakemake-tutorial\n\n";

print Profile "cd $opts{i}/snakemake/snakemake-tutorial\n\n";



#print Profile "wget https://bitbucket.org/snakemake/snakemake/downloads/snakemake-tutorial-data.tar.gz\n\n";

print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/snakemake-tutorial-data.tar.gz\n\n";


print Profile "tar -xf snakemake-tutorial-data.tar.gz\n\n";

print Profile "$opts{i}/snakemake/Miniconda3/bin/conda create -y -n snakemake-tutorial -c bioconda --file requirements.txt\n\n";

print Profile "export PATH=$opts{i}/snakemake/Miniconda3/envs/snakemake-tutorial/bin:\$PATH\n\n";



###For my_perl


print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/my_perl.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/my_perl.tar.gz\n\n";



###For my_R


print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/my_R.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/my_R.tar.gz\n\n";



###For genome_data


print Profile "wget http://downloads.sourceforge.net/project/hpprna-dependencies/hppRNA_genome.tar.gz -P $opts{i}\n\n";
print Profile "cd $opts{i}\n\n";
print Profile "tar zxvf $opts{i}/hppRNA_genome.tar.gz\n\n";



print Profile "\n\n\n";



}



close Profile;



