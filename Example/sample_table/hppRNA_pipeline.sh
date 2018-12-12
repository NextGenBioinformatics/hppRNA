#!/bin/bash -l

# 1. Force bash as the executing shell.
#$ -S /bin/bash

# 2. Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=50:0:0

# 3. Request 1 gigabyte of RAM 
#$ -l mem=10G

# 4. Request 15 gigabyte of TMPDIR space
#$ -l tmpfs=10G

# 5. Set the name of the job.
#$ -N hppRNA_data

# 6. Set the working directory to somewhere in your space.
#$ -wd /location_name_example/output

# 7. Run the application.



export PATH=/data/hppRNA_software/snakemake/Miniconda3/envs/snakemake-tutorial/bin:$PATH

export PATH=/data/hppRNA_software/perl-5.22.2/perl-5.22.2/bin:$PATH

export PATH=/data/hppRNA_software/jdk1.8.0_91/bin:$PATH

export PATH=/data/hppRNA_software/Python-2.7.11:$PATH

export PATH=/data/hppRNA_software/R-3.3.0/R-3.3.0/bin:$PATH

export PATH=/data/hppRNA_software/FastQC:$PATH

export PATH=/data/hppRNA_software/fastx_toolkit_0.0.13:$PATH

export PATH=/data/hppRNA_software/Python-2.7.11/Python-2.7.11/bin:$PATH

export PATH=/data/hppRNA_software/bowtie2-2.2.9:$PATH

export PATH=/data/hppRNA_software/tophat-2.1.1.Linux_x86_64:$PATH

export PATH=/data/hppRNA_software/cufflinks-2.2.1.Linux_x86_64:$PATH

export PATH=/data/hppRNA_software/samtools-1.3.1/samtools-1.3.1/bin:$PATH

export PATH=/data/hppRNA_software/ngsplot-develop/bin:$PATH

export NGSPLOT=/data/hppRNA_software/ngsplot-develop;

export PATH=/data/hppRNA_software/hisat2-2.0.3-beta:$PATH

export PATH=/data/hppRNA_software/stringtie-1.2.2.Linux_x86_64:$PATH

export PATH=/data/hppRNA_software/kallisto_linux-v0.42.5:$PATH

export PATH=/data/hppRNA_software/STAR-2.5.1b/bin/Linux_x86_64_static:$PATH

export PATH=/data/hppRNA_software/RSEM-1.2.30:$PATH

export PATH=/data/hppRNA_software/RSEM-1.2.30/EBSeq:$PATH

export PATH=/data/hppRNA_software/express-1.5.1-linux_x86_64:$PATH

export PATH=/data/hppRNA_software/FusionCatcher/fusioncatcher/bin:$PATH

export PATH=/data/hppRNA_software/bowtie-1.1.2:$PATH

export PATH=/data/hppRNA_software/iSeeRNA-1.2.2:$PATH



snakemake -p --cores 32 --snakefile /data/hppRNA_data/workflow_1_protein_coding_paired.snakemake



