mkdir -p /data/hppRNA_software

wget http://downloads.sourceforge.net/project/hpprna-dependencies/perl-5.22.2.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/perl-5.22.2.tar.gz

cd /data/hppRNA_software/perl-5.22.2

sh Configure -de -Dprefix=/data/hppRNA_software/perl-5.22.2/perl-5.22.2

make

make install

export PATH=/data/hppRNA_software/perl-5.22.2/perl-5.22.2/bin:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/jdk-8u91-linux-x64.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/jdk-8u91-linux-x64.tar.gz

export PATH=/data/hppRNA_software/jdk1.8.0_91/bin:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/Python-2.7.11.tgz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/Python-2.7.11.tgz

cd /data/hppRNA_software/Python-2.7.11

./configure --prefix=/data/hppRNA_software/Python-2.7.11/Python-2.7.11

make

make install

export PATH=/data/hppRNA_software/Python-2.7.11:$PATH

export PATH=/data/hppRNA_software/Python-2.7.11/Python-2.7.11/bin:$PATH

mkdir /data/hppRNA_software/setuptools

wget --no-check-certificate https://bootstrap.pypa.io/ez_setup.py -P /data/hppRNA_software/setuptools

cd /data/hppRNA_software/setuptools

python ez_setup.py --insecure

wget http://downloads.sourceforge.net/project/hpprna-dependencies/numpy-1.11.0.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/numpy-1.11.0.tar.gz

cd /data/hppRNA_software/numpy-1.11.0

python setup.py build -j 4 install

wget http://downloads.sourceforge.net/project/hpprna-dependencies/openssl-1.0.2h.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/openssl-1.0.2h.tar.gz

cd /data/hppRNA_software/openssl-1.0.2h

./config --prefix=/data/hppRNA_software/openssl-1.0.2h/openssl-1.0.2h shared -fPIC

make depend

make

make install

export PATH=/data/hppRNA_software/openssl-1.0.2h/openssl-1.0.2h/bin:$PATH

export PATH=/data/hppRNA_software/openssl-1.0.2h:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/zlib-1.2.8.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/zlib-1.2.8.tar.gz

cd /data/hppRNA_software/zlib-1.2.8

./configure --prefix=/data/hppRNA_software/zlib-1.2.8

make

make install

wget http://downloads.sourceforge.net/project/hpprna-dependencies/bzip2-1.0.6.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/bzip2-1.0.6.tar.gz

cd /data/hppRNA_software/bzip2-1.0.6

make -f Makefile-libbz2_so

make clean

make

make install PREFIX=/data/hppRNA_software/bzip2-1.0.6

wget http://downloads.sourceforge.net/project/hpprna-dependencies/xz-5.2.2.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/xz-5.2.2.tar.gz

cd /data/hppRNA_software/xz-5.2.2

./configure --prefix=/data/hppRNA_software/xz-5.2.2

make -j3

make install

wget http://downloads.sourceforge.net/project/hpprna-dependencies/pcre-8.38.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/pcre-8.38.tar.gz

cd /data/hppRNA_software/pcre-8.38

./configure --enable-utf8 --prefix=/data/hppRNA_software/pcre-8.38

make -j3

make install

wget http://downloads.sourceforge.net/project/hpprna-dependencies/libtool-2.4.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/libtool-2.4.tar.gz

cd /data/hppRNA_software/libtool-2.4

./configure --prefix=/data/hppRNA_software/libtool-2.4

make

make install

export PATH=/data/hppRNA_software/libtool-2.4/bin:$PATH

PKG_CONFIG_PATH=/data/hppRNA_software/openssl-1.0.2h/openssl-1.0.2h/bin:$PKG_CONFIG_PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/curl-7.48.0.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/curl-7.48.0.tar.gz

cd /data/hppRNA_software/curl-7.48.0

./configure --prefix=/data/hppRNA_software/curl-7.48.0/curl-7.48.0 --with-ssl

make

make install

export PATH=/data/hppRNA_software/curl-7.48.0/curl-7.48.0/bin:$PATH

export PATH=/data/hppRNA_software/bzip2-1.0.6/bin:$PATH

export PATH=/data/hppRNA_software/xz-5.2.2/bin:$PATH

export PATH=/data/hppRNA_software/zlib-1.2.8:$PATH

export PATH=/data/hppRNA_software/pcre-8.38/bin:$PATH

export LD_LIBRARY_PATH=/data/hppRNA_software/curl-7.48.0/curl-7.48.0/lib:$LD_LIBRARY_PATH

export CFLAGS="-I/data/hppRNA_software/zlib-1.2.8/include -I/data/hppRNA_software/bzip2-1.0.6/include -I/data/hppRNA_software/xz-5.2.2/include -I/data/hppRNA_software/pcre-8.38/include -I/data/hppRNA_software/libtool-2.4/include -I/data/hppRNA_software/openssl-1.0.2h/openssl-1.0.2h/include"

export LDFLAGS="-L/data/hppRNA_software/zlib-1.2.8/lib -L/data/hppRNA_software/bzip2-1.0.6/lib -L/data/hppRNA_software/xz-5.2.2/lib -L/data/hppRNA_software/pcre-8.38/lib -L/data/hppRNA_software/libtool-2.4/lib -L/data/hppRNA_software/openssl-1.0.2h/openssl-1.0.2h/lib"

wget http://downloads.sourceforge.net/project/hpprna-dependencies/R-3.3.0.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/R-3.3.0.tar.gz

cd /data/hppRNA_software/R-3.3.0

./configure --prefix=/data/hppRNA_software/R-3.3.0/R-3.3.0

make

make install

export PATH=/data/hppRNA_software/R-3.3.0/R-3.3.0/bin:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/fastqc_v0.11.5.zip -P /data/hppRNA_software

cd /data/hppRNA_software

unzip /data/hppRNA_software/fastqc_v0.11.5.zip

cd /data/hppRNA_software/FastQC

chmod 0755 fastqc

export PATH=/data/hppRNA_software/FastQC:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 -P /data/hppRNA_software

cd /data/hppRNA_software

tar jxf /data/hppRNA_software/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2

mv bin fastx_toolkit_0.0.13

export PATH=/data/hppRNA_software/fastx_toolkit_0.0.13:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/prinseq-lite-0.20.4.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/prinseq-lite-0.20.4.tar.gz

wget http://downloads.sourceforge.net/project/hpprna-dependencies/cutadapt-1.9.1.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/cutadapt-1.9.1.tar.gz

cd /data/hppRNA_software/cutadapt-1.9.1

python setup.py install

export PATH=/data/hppRNA_software/Python-2.7.11/Python-2.7.11/bin:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/bowtie2-2.2.9-linux-x86_64.zip -P /data/hppRNA_software

cd /data/hppRNA_software

unzip /data/hppRNA_software/bowtie2-2.2.9-linux-x86_64.zip

export PATH=/data/hppRNA_software/bowtie2-2.2.9:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/tophat-2.1.1.Linux_x86_64.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/tophat-2.1.1.Linux_x86_64.tar.gz

export PATH=/data/hppRNA_software/tophat-2.1.1.Linux_x86_64:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/cufflinks-2.2.1.Linux_x86_64.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/cufflinks-2.2.1.Linux_x86_64.tar.gz

export PATH=/data/hppRNA_software/cufflinks-2.2.1.Linux_x86_64:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/samtools-1.3.1.tar.bz2 -P /data/hppRNA_software

cd /data/hppRNA_software

tar jxf /data/hppRNA_software/samtools-1.3.1.tar.bz2

cd /data/hppRNA_software/samtools-1.3.1

./configure --prefix=/data/hppRNA_software/samtools-1.3.1/samtools-1.3.1

make

make install

export PATH=/data/hppRNA_software/samtools-1.3.1/samtools-1.3.1/bin:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/BlackOPs_v1.06.zip -P /data/hppRNA_software

cd /data/hppRNA_software

unzip /data/hppRNA_software/BlackOPs_v1.06.zip

wget http://downloads.sourceforge.net/project/hpprna-dependencies/R_package.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/R_package.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/Rsubread_1.22.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/BiocGenerics_0.18.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/S4Vectors_0.10.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/IRanges_2.6.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/GenomeInfoDb_1.8.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/zlibbioc_1.18.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/XVector_0.12.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/GenomicRanges_1.24.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/Biobase_2.32.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/SummarizedExperiment_1.2.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/lambda.r_1.1.6.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/futile.options_1.0.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/futile.logger_1.4.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/snow_0.3-13.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/BiocParallel_1.6.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/DBI_0.3.1.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/RSQLite_0.11.4.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/AnnotationDbi_1.34.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/XML_3.98-1.3.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/xtable_1.8-0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/annotate_1.50.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/genefilter_1.54.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/locfit_1.5-9.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/RColorBrewer_1.0-5.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/geneplotter_1.50.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/digest_0.6.8.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/gtable_0.1.2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/Rcpp_0.12.3.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/plyr_1.8.2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/stringi-master

R CMD INSTALL /data/hppRNA_software/R_package/magrittr_1.0.1.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/stringr_0.6.2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/reshape2_1.4.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/dichromat_1.2-4.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/colorspace_1.2-6.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/munsell_0.4.2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/labeling_0.2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/scales_0.3.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/ggplot2_2.0.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/Formula_1.2-0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/latticeExtra_0.6-26.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/acepack_1.3-3.2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/gridExtra_2.2.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/chron_2.3-46.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/data.table_1.9.4.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/Hmisc_3.17-3.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/RcppArmadillo_0.4.500.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/DESeq2_1.12.0.tar.gz

wget http://downloads.sourceforge.net/project/hpprna-dependencies/2.5.1b.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/2.5.1b.tar.gz

export PATH=/data/hppRNA_software/STAR-2.5.1b/bin/Linux_x86_64_static:$PATH

export LD_LIBRARY_PATH=

export CFLAGS=

export LDFLAGS=

wget http://downloads.sourceforge.net/project/hpprna-dependencies/v1.2.30.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/v1.2.30.tar.gz

cd /data/hppRNA_software/RSEM-1.2.30

make

make ebseq

export PATH=/data/hppRNA_software/RSEM-1.2.30:$PATH

export PATH=/data/hppRNA_software/RSEM-1.2.30/EBSeq:$PATH

R CMD INSTALL /data/hppRNA_software/R_package/blockmodeling_0.1.7.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/gtools_3.4.2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/gdata_2.16.1.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/bitops_1.0-5.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/caTools_1.17.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/gplots_3.0.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/memoise_0.2.1.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/crayon_1.3.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/praise_1.0.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/R6_2.1.1.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/testthat_1.0.1.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/EBSeq_1.12.0.tar.gz

wget http://downloads.sourceforge.net/project/hpprna-dependencies/bowtie-1.1.2-linux-x86_64.zip -P /data/hppRNA_software

cd /data/hppRNA_software

unzip /data/hppRNA_software/bowtie-1.1.2-linux-x86_64.zip

export PATH=/data/hppRNA_software/bowtie-1.1.2:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/express-1.5.1-linux_x86_64.tgz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/express-1.5.1-linux_x86_64.tgz

export PATH=/data/hppRNA_software/express-1.5.1-linux_x86_64:$PATH

R CMD INSTALL /data/hppRNA_software/R_package/limma_3.28.2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/edgeR_3.14.0.tar.gz

wget http://downloads.sourceforge.net/project/hpprna-dependencies/kallisto_linux-v0.42.5.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/kallisto_linux-v0.42.5.tar.gz

export PATH=/data/hppRNA_software/kallisto_linux-v0.42.5:$PATH

R CMD INSTALL /data/hppRNA_software/R_package/memoise_1.0.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/git2r_0.15.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/jsonlite_1.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/rstudioapi_0.6.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/whisker_0.3-2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/withr_1.0.2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/mime_0.5.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/curl_1.2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/openssl_0.9.4.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/httr_1.2.1.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/devtools_1.11.0.tar.gz

export PATH=/data/hppRNA_software/perl-5.22.2/perl-5.22.2/bin:$PATH

export PATH=/data/hppRNA_software/jdk1.8.0_91/bin:$PATH

export PATH=/data/hppRNA_software/openssl-1.0.2h/openssl-1.0.2h/bin:$PATH

export PATH=/data/hppRNA_software/openssl-1.0.2h:$PATH

export PATH=/data/hppRNA_software/libtool-2.4/bin:$PATH

export PATH=/data/hppRNA_software/curl-7.48.0/curl-7.48.0/bin:$PATH

export PATH=/data/hppRNA_software/bzip2-1.0.6/bin:$PATH

export PATH=/data/hppRNA_software/xz-5.2.2/bin:$PATH

export PATH=/data/hppRNA_software/zlib-1.2.8:$PATH

export PATH=/data/hppRNA_software/pcre-8.38/bin:$PATH

export LD_LIBRARY_PATH=/data/hppRNA_software/curl-7.48.0/curl-7.48.0/lib:$LD_LIBRARY_PATH

export CFLAGS="-I/data/hppRNA_software/zlib-1.2.8/include -I/data/hppRNA_software/bzip2-1.0.6/include -I/data/hppRNA_software/xz-5.2.2/include -I/data/hppRNA_software/pcre-8.38/include -I/data/hppRNA_software/libtool-2.4/include -I/data/hppRNA_software/openssl-1.0.2h/openssl-1.0.2h/include"

export LDFLAGS="-L/data/hppRNA_software/zlib-1.2.8/lib -L/data/hppRNA_software/bzip2-1.0.6/lib -L/data/hppRNA_software/xz-5.2.2/lib -L/data/hppRNA_software/pcre-8.38/lib -L/data/hppRNA_software/libtool-2.4/lib -L/data/hppRNA_software/openssl-1.0.2h/openssl-1.0.2h/lib"

export PATH=/data/hppRNA_software/Python-2.7.11:$PATH

export PATH=/data/hppRNA_software/R-3.3.0/R-3.3.0/bin:$PATH

mkdir /data/hppRNA_software/sleuth

cd /data/hppRNA_software/sleuth

echo "source(\"http://bioconductor.org/biocLite.R\")" >> sleuth_install.R

echo "biocLite(\"rhdf5\")" >> sleuth_install.R

echo "devtools::install_github(\"pachterlab/sleuth\")" >> sleuth_install.R

R --slave --vanilla < /data/hppRNA_software/sleuth/sleuth_install.R

wget http://downloads.sourceforge.net/project/hpprna-dependencies/hisat2-2.0.3-beta-Linux_x86_64.zip -P /data/hppRNA_software

cd /data/hppRNA_software

unzip /data/hppRNA_software/hisat2-2.0.3-beta-Linux_x86_64.zip

export PATH=/data/hppRNA_software/hisat2-2.0.3-beta:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/stringtie-1.2.2.Linux_x86_64.tar.gz -P /data/hppRNA_software

tar zxvf /data/hppRNA_software/stringtie-1.2.2.Linux_x86_64.tar.gz

export PATH=/data/hppRNA_software/stringtie-1.2.2.Linux_x86_64:$PATH

R CMD INSTALL /data/hppRNA_software/R_package/sva_3.20.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/Biostrings_2.40.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/RCurl_1.95-4.7.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/Rsamtools_1.24.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/GenomicAlignments_1.8.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/rtracklayer_1.32.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/ballgown_2.4.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/iterators_1.0.7.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/foreach_1.4.2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/doMC_1.3.3.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/caTools_1.17.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/BSgenome_1.40.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/hwriter_1.3.1.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/ShortRead_1.30.0.tar.gz

wget http://downloads.sourceforge.net/project/hpprna-dependencies/develop.zip -P /data/hppRNA_software

cd /data/hppRNA_software

unzip /data/hppRNA_software/develop.zip

rm /data/hppRNA_software/develop.zip

export PATH=/data/hppRNA_software/ngsplot-develop/bin:$PATH

export NGSPLOT=/data/hppRNA_software/ngsplot-develop;

cd /data/hppRNA_software/ngsplot-develop

wget http://downloads.sourceforge.net/project/hpprna-dependencies/ngsplotdb_hg19_75_3.00.tar.gz -P /data/hppRNA_software/ngsplot-develop

wget http://downloads.sourceforge.net/project/hpprna-dependencies/ngsplotdb_mm10_75_3.00.tar.gz -P /data/hppRNA_software/ngsplot-develop

ngsplotdb.py -y install /data/hppRNA_software/ngsplot-develop/ngsplotdb_hg19_75_3.00.tar.gz

ngsplotdb.py -y install /data/hppRNA_software/ngsplot-develop/ngsplotdb_mm10_75_3.00.tar.gz

mkdir /data/hppRNA_software/wig2bigwig

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig -P /data/hppRNA_software/wig2bigwig

cd /data/hppRNA_software/wig2bigwig

chmod 0755 wigToBigWig

export PATH=/data/hppRNA_software/wig2bigwig:$PATH

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes -P /data/hppRNA_software/wig2bigwig

chmod 0755 /data/hppRNA_software/wig2bigwig/fetchChromSizes

/data/hppRNA_software/wig2bigwig/fetchChromSizes hg19 > /data/hppRNA_software/wig2bigwig/hg19.chrom.sizes

/data/hppRNA_software/wig2bigwig/fetchChromSizes mm10 > /data/hppRNA_software/wig2bigwig/mm10.chrom.sizes

R CMD INSTALL /data/hppRNA_software/R_package/preprocessCore_1.34.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/proto_0.3-9.2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/assertthat_0.1.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/lazyeval_0.1.9.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/BH_1.60.0-1.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/dplyr_0.4.2.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/tidyr_0.4.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/ggfortify_0.1.0.tar.gz

R CMD INSTALL /data/hppRNA_software/R_package/cluster_2.0.3.tar.gz

wget http://downloads.sourceforge.net/project/hpprna-dependencies/iSeeRNA-1.2.2.x86_64.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/iSeeRNA-1.2.2.x86_64.tar.gz

cd /data/hppRNA_software/iSeeRNA-1.2.2

sh auto_download_data.sh hg19

sh auto_download_data.sh mm10

wget http://downloads.sourceforge.net/project/hpprna-dependencies/v0.4.4.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/v0.4.4.tar.gz

cd /data/hppRNA_software/DCC-0.4.4

python setup.py install

R CMD INSTALL /data/hppRNA_software/R_package/aod_1.3.tar.gz

mkdir /data/hppRNA_software/circTest

cd /data/hppRNA_software/circTest

echo "require(devtools)" >> circTest_install.R

echo "install_github(\"dieterich-lab/CircTest\")" >> circTest_install.R

echo "library(CircTest)" >> circTest_install.R

R --slave --vanilla < /data/hppRNA_software/circTest/circTest_install.R

wget http://downloads.sourceforge.net/project/hpprna-dependencies/picard-tools-2.3.0.zip -P /data/hppRNA_software

cd /data/hppRNA_software

unzip /data/hppRNA_software/picard-tools-2.3.0.zip

mkdir /data/hppRNA_software/FusionCatcher

cd /data/hppRNA_software/FusionCatcher

wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py && python bootstrap.py -t --download -y

export PATH=/data/hppRNA_software/FusionCatcher/fusioncatcher/bin:$PATH

fusioncatcher-build --ftp-ensembl-path=pub/release-80 -g mus_musculus -o /data/hppRNA_software/FusionCatcher/fusioncatcher/data/mus_musculus

fusioncatcher-build -g mus_musculus -o /data/hppRNA_software/FusionCatcher/fusioncatcher/data/mus_musculus

mkdir /data/hppRNA_software/snakemake

cd /data/hppRNA_software/snakemake

wget http://downloads.sourceforge.net/project/hpprna-dependencies/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh -b -p /data/hppRNA_software/snakemake/Miniconda3

export PATH=/data/hppRNA_software/snakemake/Miniconda3/bin:$PATH

mkdir /data/hppRNA_software/snakemake/snakemake-tutorial

cd /data/hppRNA_software/snakemake/snakemake-tutorial

wget http://downloads.sourceforge.net/project/hpprna-dependencies/snakemake-tutorial-data.tar.gz

tar -xf snakemake-tutorial-data.tar.gz

/data/hppRNA_software/snakemake/Miniconda3/bin/conda create -y -n snakemake-tutorial -c bioconda --file requirements.txt

export PATH=/data/hppRNA_software/snakemake/Miniconda3/envs/snakemake-tutorial/bin:$PATH

wget http://downloads.sourceforge.net/project/hpprna-dependencies/my_perl.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/my_perl.tar.gz

wget http://downloads.sourceforge.net/project/hpprna-dependencies/my_R.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/my_R.tar.gz

wget http://downloads.sourceforge.net/project/hpprna-dependencies/hppRNA_genome.tar.gz -P /data/hppRNA_software

cd /data/hppRNA_software

tar zxvf /data/hppRNA_software/hppRNA_genome.tar.gz




