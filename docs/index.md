hppRNA - a Snakemake-based handy parameter-free pipeline for RNA-Seq analysis of numerous samples


hppRNA package is dedicated to the RNA-Seq analysis for a large number of samples simultaneously from the very beginning to the very end, which is formulated in Snakemake pipeline management system. It starts from fastq files and will produce gene/isoform expression matrix, differentially-expressed-genes, sample clusters as well as detection of SNP and fusion genes by combination of the state-of-the-art software.

The first version handles protein-coding genes, lncRNAs and circRNAs and includes six core-workflows.

* (1) Tophat - Cufflink - Cuffdiff;
* (2) Subread - featureCounts - DESeq2;
* (3) STAR - RSEM - EBSeq;
* (4) Bowtie - eXpress - edgeR;
* (5) kallisto - sleuth;
* (6) HISAT - StringTie - Ballgown.

Please cite the following paper when using this package:

Dapeng Wang. hppRNA—a Snakemake-based handy parameter-free pipeline for RNA-Seq analysis of numerous samples. Briefings in Bioinformatics, Volume 19, Issue 4, 20 July 2018, Pages 622–626.

