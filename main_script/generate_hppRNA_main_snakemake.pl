#!/usr/bin/perl -w
# Program name:	generate_hppRNA_main_snakemake.pl
# Programmer:	Dapeng Wang
# Email:	wangdp123@gmail.com
# Date:	2016-12-01
# Last update:	2018-08-04
# Version:	v1.3.6



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
	-i	:	sample_table
	-o	:	NGS analysis snakemake file


USAGE

die $usage unless $opts{i};
die $usage unless $opts{o};


#####################################################################################Essential parameter definition start##########################################################################################



my @Sample=();
my @Sample_name=();
my @real_sample_name=();
my @real_group_name=();
my @read_type=();


my %read_length=();

my $real_project_name="";
my $Project="";
my $Seq="";
my $Type="";
my $Assembly_Version="";
my $Working_FOLDER="";
my $Analysis_type="";
my $Core_workflow="";
my $Core_workflow_name="";

my $Adaptor="";
my $Fusion="";
my $SNP="";



my $PERL_FOLDER="";
my $R_FOLDER="";
my $Contaminants_FILE="";
my $Adapters_FILE="";
my $Prinseq_lite_FILE="";
my $Bowtie2Index_FILE="";
my $hisat2Index_FILE="";
my $Extract_unique_FILE="";
my $BowtieIndex_FILE="";

my $kallisto_FILE="";
my $kallisto_transcript2gene="";
my $FusionCatcher_data="";

my $hppRNA_genome_fasta_file="";
my $hppRNA_genome_gtf_file="";
my $hppRNA_transcript_fasta_file="";

my $Software_FOLDER="";
my $GATK_FOLDER="";



my $Thread="";
my $Tophat_mate_inner_dist="";
my $Tophat_mate_std_dev="";
my $Tophat_Cufflinks_Cuffdiff_library_type="";
my $Subread_minFragLength="";
my $Subread_maxFragLength="";
my $Subread_PE_orientation="";
my $featureCounts_strandSpecific="";
my $Bowtie_HISAT_minins="";
my $Bowtie_HISAT_maxins="";
my $Bowtie_HISAT_PE_orientation="";
my $Kallisto_l="";
my $Kallisto_s="";
my $RSEM_fragment_length_min="";
my $RSEM_fragment_length_max="";
my $RSEM_fragment_length_mean="";
my $RSEM_fragment_length_sd="";
my $RSEM_forward_prob="";
my $Bowtie_HISAT_nofw_norc="";
my $eXpress_orientation="";
my $DCC_strand="";



my $Tophat_fusion_post_FILE="";

my $Species="";
my $Directory="";
my $result_folder="";
my $prefix="";



my %species_version=();
my %DESeq2_species_specific_gene_name_list=();
my %DESeq2_featurecount_gtf_file=();

my %bowtie1_index_species=();
my %bowtie2_index_species=();
my %gtf_file_species=();

my %bed_file=();

#my $build_folders="";
#my $delete_intermediate_files="";


my @name_fastq_file=();
my @name_adapt_fastq_file=();
my @name_adapt_qc_fastq_file=();

my @name_adapt_qc_collapser_fasta_file=();

my @name_qc_bam_file=();
my @name_qc_sam_file=();
my @name_qc_10_bam_file=();
my @name_qc_10_sorted_file=();
my @name_qc_10_sorted_bam_file=();

my @treatment=();
my @control=();

my %is_control=();

my $treatment_bam="";
my $control_bam="";

############################################
############################################

my @platform=();

my @name_thout_file=();

my @name_fusion_thout_file=();

my @name_clout_file=();
my @name_align_summary_file=();

my @name_sorted_file=();
my @name_sorted_bam_file=();
my @name_sorted_sam_file=();
my @name_sorted_rmdup_bam_file=();
my @name_sorted_unique_sam_file=();
my @name_sorted_unique_bam_file=();

my @name_sorted_unique_rmdup_bam_file=();

my @name_pre_snp_bam_file=();

my @name_snp_bam_file=();

my @snp_1st_file=();
my @snp_2nd_file=();



my @name_wig_file=();
my @name_wig_2_file=();
my @name_bw_file=();



my %replicate=();






my @file_name=();
my @sample_name=();
my @group_name=();
my %group_name_hash=();
my %group=();





open I,"$opts{i}";

while (<I>)

{

chomp;

my @t=split(/\t/,$_);

if ($t[0] eq "Project") {$Project=$t[1];}

if ($t[0] eq "Seq") {$Seq=$t[1];}

if ($Seq eq "RNA-Seq") {$prefix="R";}

if ($t[0] eq "Type") {$Type=$t[1];}

if ($t[0] eq "Assembly_Version") {$Assembly_Version=$t[1];}

if ($t[0] eq "Adaptor") {$Adaptor=$t[1];}

if ($t[0] eq "Analysis_type") {$Analysis_type=$t[1];}

if ($t[0] eq "Core_workflow") {$Core_workflow=$t[1];}



if ($Core_workflow eq 1) {$Core_workflow_name="Tophat_Cufflink_Cuffdiff";}
if ($Core_workflow eq 2) {$Core_workflow_name="Subread_featureCounts_DESeq2";}
if ($Core_workflow eq 3) {$Core_workflow_name="STAR_RSEM_EBSeq";}
if ($Core_workflow eq 4) {$Core_workflow_name="Bowtie_eXpress_edgeR";}
if ($Core_workflow eq 5) {$Core_workflow_name="Kallisto_sleuth";}
if ($Core_workflow eq 6) {$Core_workflow_name="HISAT_StringTie_Ballgown";}



if ($t[0] eq "Fusion") {$Fusion=$t[1];}

if ($t[0] eq "SNP") {$SNP=$t[1];}

if ($t[0] eq "Species") {$Species=$t[1];}

if ($t[0] eq "Working_FOLDER") {$Working_FOLDER=$t[1];}

if ($t[0] eq "Output") {$result_folder=$t[1];}

if ($t[0] eq "Software_FOLDER") {$Software_FOLDER=$t[1];}

if ($t[0] eq "GATK_FOLDER") {$GATK_FOLDER=$t[1];}



############################start the arguments for each mapping software



if ($t[0] eq "Thread") {$Thread=$t[1];}
if ($t[0] eq "Tophat: mate-inner-dist") {$Tophat_mate_inner_dist=$t[1];}
if ($t[0] eq "Tophat: mate-std-dev") {$Tophat_mate_std_dev=$t[1];}



if ($t[0] eq "Tophat/Cufflinks/Cuffdiff: library-type") {$Tophat_Cufflinks_Cuffdiff_library_type=$t[1];}



if ($t[0] eq "Subread: minFragLength") {$Subread_minFragLength=$t[1];}
if ($t[0] eq "Subread: maxFragLength") {$Subread_maxFragLength=$t[1];}
if ($t[0] eq "Subread: PE_orientation") {$Subread_PE_orientation=$t[1];}



if ($t[0] eq "featureCounts: strandSpecific") {$featureCounts_strandSpecific=$t[1];}



if ($t[0] eq "Bowtie/HISAT: minins") {$Bowtie_HISAT_minins=$t[1];}
if ($t[0] eq "Bowtie/HISAT: maxins") {$Bowtie_HISAT_maxins=$t[1];}
if ($t[0] eq "Bowtie/HISAT: PE_orientation") {$Bowtie_HISAT_PE_orientation=$t[1];}



if ($t[0] eq "Bowtie/HISAT: nofw/norc") {

if (!(defined $t[1])) {

$Bowtie_HISAT_nofw_norc=" ";

}

else

{

$Bowtie_HISAT_nofw_norc=$t[1];

}

}



if ($t[0] eq "eXpress: orientation") {

if (!(defined $t[1])) {

$eXpress_orientation=" ";

}

else

{

$eXpress_orientation=$t[1];

}

}



if ($t[0] eq "Kallisto: l") {$Kallisto_l=$t[1];}
if ($t[0] eq "Kallisto: s") {$Kallisto_s=$t[1];}

if ($t[0] eq "DCC: strand") {$DCC_strand=$t[1];}

if ($t[0] eq "RSEM: fragment-length-min") {$RSEM_fragment_length_min=$t[1];}
if ($t[0] eq "RSEM: fragment-length-max") {$RSEM_fragment_length_max=$t[1];}
if ($t[0] eq "RSEM: fragment-length-mean") {$RSEM_fragment_length_mean=$t[1];}
if ($t[0] eq "RSEM: fragment-length-sd") {$RSEM_fragment_length_sd=$t[1];}
if ($t[0] eq "RSEM: forward-prob") {$RSEM_forward_prob=$t[1];}



############################################################################



#if ($t[0] eq "Build_folders") {$build_folders=$t[1];}
#if ($t[0] eq "Delete_intermediate_files") {$delete_intermediate_files=$t[1];}



if ($t[0] eq "Sample") {

my $sample_full_name=$prefix."_".$t[6];

push(@file_name,$t[5]);
push(@sample_name,$sample_full_name);
push(@group_name,$t[8]);
push(@{$replicate{$t[8]}},$sample_full_name);
$group_name_hash{$t[8]}=1;
$group{$sample_full_name}=$t[8];


} #if ($t[0] eq "Sample") {


if ($t[0] eq "DEG") {

push(@control,$t[1]);

push(@treatment,$t[2]);

}  # if ($t[0] eq "DEG") {



} #end of while

close I;



#####################################################################################Essential parameter definition end##########################################################################################


#####################################################################################Snakemake start##########################################################################################

###############################################Variables for snakemake setup start############################################################

###########################################Sample information start#################################################

open Profile,">$opts{o}";



{



###for samples



################################################################################
my $samples_1="\"".$file_name[0]."\"";


if (@file_name>1) {

for (my $i=1;$i<@file_name;$i++) {$samples_1.=", \"".$file_name[$i]."\"";}

}
################################################################################
################################################################################
my $samples_2="\"".$sample_name[0]."\"";

if (@sample_name>1) {

for (my $i=1;$i<@sample_name;$i++) {$samples_2.=", \"".$sample_name[$i]."\"";}

}
################################################################################
################################################################################
my $samples_join=$sample_name[0].";";


if (@sample_name>1) {

for (my $i=1;$i<@sample_name;$i++) {$samples_join.=$sample_name[$i].";";}

}
################################################################################


################################################################################
my $samples_group_join=$group{$sample_name[0]}.";";


if (@sample_name>1) {

for (my $i=1;$i<@sample_name;$i++) {$samples_group_join.=$group{$sample_name[$i]}.";";}

}
################################################################################



print Profile "\n\n\n";

print Profile "SAMPLES_1 = [".$samples_1."]\n\n";

print Profile "SAMPLES_2 = [".$samples_2."]\n\n";

print Profile "SAMPLE_join = \"".$samples_join."\"\n\n";

print Profile "SAMPLE_group_join = \"".$samples_group_join."\"\n\n";



#SAMPLES_1 = ["ERR315477", "ERR315455", "ERR315432", "ERR315352", "ERR315456", "ERR315391"]

#SAMPLES_2 = ["R_brain_3b", "R_brain_3c", "R_brain_a", "R_testis_7a", "R_testis_7b", "R_testis_7c"]

#SAMPLE_join = "R_brain_3b;R_brain_3c;R_brain_a;R_testis_7a;R_testis_7b;R_testis_7c;"

#SAMPLE_group_join = "brain;brain;brain;testis;testis;testis;"



print Profile "\n\n\n";



foreach my $key (keys %group_name_hash) {



################################################################################
my $samples_array_1="\"".$replicate{$key}[0]."\"";

if (@{$replicate{$key}}>1) {

for (my $j=1;$j<@{$replicate{$key}};$j++) {$samples_array_1.=", \"".$replicate{$key}[$j]."\"";}

}
################################################################################



print Profile "CLASS_".$key." = [".$samples_array_1."]\n\n";



#CLASS_testis = ["R_testis_7a", "R_testis_7b", "R_testis_7c"]

#CLASS_brain = ["R_brain_3b", "R_brain_3c", "R_brain_a"]



################################################################################
my $samples_array_2=$replicate{$key}[0];


if (@{$replicate{$key}}>1) {

for (my $j=1;$j<@{$replicate{$key}};$j++) {$samples_array_2.=";".$replicate{$key}[$j];}

}
################################################################################


################################################################################



print Profile "CLASS_".$key."_2 = \"".$samples_array_2."\"\n\n";



#CLASS_testis_2 = "R_testis_7a;R_testis_7b;R_testis_7c"

#CLASS_brain_2 = "R_brain_3b;R_brain_3c;R_brain_a"



}###foreach my $key (keys %group_name_hash) {



###for comparisons



for (my $i=0;$i<@treatment;$i++) {



print Profile "CLASS_".$control[$i]."_vs_".$treatment[$i]." = CLASS_".$control[$i]." + CLASS_".$treatment[$i]."\n\n";



#PAIRS = ["testis_vs_brain"]

#CLASS_testis_vs_brain = CLASS_testis + CLASS_brain



################################################################################
################################################################################
my $samples_join_treatment=$replicate{$treatment[$i]}[0].";";


if (@{$replicate{$treatment[$i]}}>1) {

for (my $j=1;$j<@{$replicate{$treatment[$i]}};$j++) {$samples_join_treatment.=$replicate{$treatment[$i]}[$j].";";}

}
################################################################################
################################################################################
################################################################################
my $samples_join_control=$replicate{$control[$i]}[0].";";


if (@{$replicate{$control[$i]}}>1) {

for (my $j=1;$j<@{$replicate{$control[$i]}};$j++) {$samples_join_control.=$replicate{$control[$i]}[$j].";";}

}
################################################################################

my $samples_join_control_vs_treatment=$samples_join_control.$samples_join_treatment;

################################################################################




################################################################################
################################################################################
my $samples_group_join_treatment=$group{$replicate{$treatment[$i]}[0]}.";";


if (@{$replicate{$treatment[$i]}}>1) {

for (my $j=1;$j<@{$replicate{$treatment[$i]}};$j++) {$samples_group_join_treatment.=$group{$replicate{$treatment[$i]}[$j]}.";";}

}
################################################################################
################################################################################
################################################################################
my $samples_group_join_control=$group{$replicate{$control[$i]}[0]}.";";


if (@{$replicate{$control[$i]}}>1) {

for (my $j=1;$j<@{$replicate{$control[$i]}};$j++) {$samples_group_join_control.=$group{$replicate{$control[$i]}[$j]}.";";}

}
################################################################################

my $samples_group_join_control_vs_treatment=$samples_group_join_control.$samples_group_join_treatment;

################################################################################



print Profile "SAMPLE_join_".$control[$i]."_vs_".$treatment[$i]." = \"".$samples_join_control_vs_treatment."\"\n\n";

print Profile "SAMPLE_group_join_".$control[$i]."_vs_".$treatment[$i]." = \"".$samples_group_join_control_vs_treatment."\"\n\n";



#SAMPLE_join_testis_vs_brain = "R_brain_3b;R_brain_3c;R_brain_a;R_testis_7a;R_testis_7b;R_testis_7c;"

#SAMPLE_group_join_testis_vs_brain = "brain;brain;brain;testis;testis;testis;"



} ###for (my $i=0;$i<@treatment;$i++) {



my $pair_array="\"".$control[0]."_vs_".$treatment[0]."\"";



if (@treatment>1) {

for (my $k=1;$k<@treatment;$k++) {$pair_array.=", \"".$control[$k]."_vs_".$treatment[$k]."\"";}

}



print Profile "PAIRS = [".$pair_array."]\n\n";



###########################################Sample information end#################################################

###########################################Variables definition start#############################################

print Profile "\n\n\n";

print Profile "PROJECT_NAME = \"".$prefix."_".$Project."\"\n\n";

print Profile "WORKFLOW_NAME = \"".$Core_workflow_name."\"\n\n";

print Profile "WORKING_FOLDER = \"".$Working_FOLDER."\"\n\n";

print Profile "HPPRNA_SOFTWARE_FOLDER = \"".$Software_FOLDER."\"\n\n";

print Profile "ADAPTER_SEQUENCE = \"".$Adaptor."\"\n\n";

print Profile "ASSEMBLY_VERSION = \"".$Assembly_Version."\"\n\n";

print Profile "GATK_FOLDER = \"".$GATK_FOLDER."\"\n\n";



if (($Analysis_type eq "protein-coding")||($Analysis_type eq "circRNA")) {

print Profile '

GTF = HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + "/genes.selected.pc.gtf"

MAP_FILE = HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + "/transcript_to_gene.txt"

';

}



else {

print Profile '

GTF = HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + "/mRNA_lncRNA.gtf"

MAP_FILE = HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + "/mRNA_lncRNA.transcript_to_gene.txt"

';

}



print Profile "\n\n\n";

###########################################Variables definition end #############################################

###########################################Constants definition start############################################
print Profile '

FASTQC_OUTPUT_FOLDER = WORKING_FOLDER + "/standard_results/FastQC"

CONTAMINANT_LIST_FILE = HPPRNA_SOFTWARE_FOLDER + "/FastQC/Configuration/contaminant_list.txt"

ADAPTER_LIST_FILE = HPPRNA_SOFTWARE_FOLDER + "/FastQC/Configuration/adapter_list.txt"

PRINSEQ_LITE_FILE = HPPRNA_SOFTWARE_FOLDER + "/prinseq-lite-0.20.4/prinseq-lite.pl"

WIGTOBIGWIG = HPPRNA_SOFTWARE_FOLDER + "/wig2bigwig/wigToBigWig"

CHROM_SIZES = HPPRNA_SOFTWARE_FOLDER + "/wig2bigwig/" + ASSEMBLY_VERSION + ".chrom.sizes"

GENOME_FA_FILE = HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + "/genome.fa"

REPEAT_GTF = HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + "/repeats.gtf"

TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE = WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION

TOPHAT_LNCRNA_INDEX_FILE = WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Index/" + ASSEMBLY_VERSION

EXTRACT_UNIQUE = HPPRNA_SOFTWARE_FOLDER + "/BlackOPs_v1.06/perl/extract_unique.pl"

';


if ($Assembly_Version eq "hg19") {

print Profile "FUSIONCATCHER_DATA_FOLDER = HPPRNA_SOFTWARE_FOLDER + \"/FusionCatcher/fusioncatcher/data/current/\"\n\n";

}


if ($Assembly_Version eq "mm10") {

print Profile "FUSIONCATCHER_DATA_FOLDER = HPPRNA_SOFTWARE_FOLDER + \"/FusionCatcher/fusioncatcher/data/mus_musculus/\"\n\n";

}



###########################################Constents definition end############################################



###########################################Software parameter definition start############################################



print Profile "THREADS = ".$Thread."\n\n";
print Profile "TOPHAT_MATE_INNER_DIST = ".$Tophat_mate_inner_dist."\n\n";
print Profile "TOPHAT_MATE_STD_DEV = ".$Tophat_mate_std_dev."\n\n";
print Profile "TOPHAT_CUFFLINKS_CUFFDIFF_LIBRARY_TYPE = \"".$Tophat_Cufflinks_Cuffdiff_library_type."\"\n\n";
print Profile "SUBREAD_MINFRAGLENGTH = ".$Subread_minFragLength."\n\n";
print Profile "SUBREAD_MAXFRAGLENGTH = ".$Subread_maxFragLength."\n\n";
print Profile "SUBREAD_PE_ORIENTATION = \"".$Subread_PE_orientation."\"\n\n";
print Profile "FEATURECOUNTS_STRANDSPECIFIC = \"".$featureCounts_strandSpecific."\"\n\n";
print Profile "BOWTIE_HISAT_MININS = ".$Bowtie_HISAT_minins."\n\n";
print Profile "BOWTIE_HISAT_MAXINS = ".$Bowtie_HISAT_maxins."\n\n";
print Profile "BOWTIE_HISAT_PE_ORIENTATION = \"".$Bowtie_HISAT_PE_orientation."\"\n\n";
print Profile "KALLISTO_L = ".$Kallisto_l."\n\n";
print Profile "KALLISTO_S = ".$Kallisto_s."\n\n";
print Profile "RSEM_FRAGMENT_LENGTH_MIN = ".$RSEM_fragment_length_min."\n\n";
print Profile "RSEM_FRAGMENT_LENGTH_MAX = ".$RSEM_fragment_length_max."\n\n";
print Profile "RSEM_FRAGMENT_LENGTH_MEAN = ".$RSEM_fragment_length_mean."\n\n";
print Profile "RSEM_FRAGMENT_LENGTH_SD = ".$RSEM_fragment_length_sd."\n\n";
print Profile "RSEM_FORWARD_PROB = ".$RSEM_forward_prob."\n\n";



print Profile "BOWTIE_HISAT_NOFW_NORC = \"".$Bowtie_HISAT_nofw_norc."\"\n\n";
print Profile "EXPRESS_ORIENTATION = \"".$eXpress_orientation."\"\n\n";
print Profile "DCC_STRAND = \"".$DCC_strand."\"\n\n";



###########################################Software parameter definition end############################################



###############################################Variables for snakemake setup end############################################################



###############################################Write up all the desired results start############################################################



###fastq (paired vs. single)



if ($Type eq "Paired-End") {



print Profile '



rule all:
    input:
        expand(WORKING_FOLDER + "/standard_results/Raw_FASTQ/{sample}.R1.fastq", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/Raw_FASTQ/{sample}.R2.fastq", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/FastQC/{sample}.R1_fastqc.zip", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/FastQC/{sample}.R1_fastqc.html", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/FastQC/{sample}.R2_fastqc.zip", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/FastQC/{sample}.R2_fastqc.html", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/Processed_FASTQ/{sample}_1.fastq", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/Processed_FASTQ/{sample}_2.fastq", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/FastQC/{sample}_1_fastqc.zip", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/FastQC/{sample}_1_fastqc.html", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/FastQC/{sample}_2_fastqc.zip", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/FastQC/{sample}_2_fastqc.html", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/Fasta/{sample}_1.collapser.fasta", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/Fasta/{sample}_2.collapser.fasta", sample = SAMPLES_2),
';



}



if ($Type eq "Single-End") {



print Profile '



rule all:
    input:
        expand(WORKING_FOLDER + "/standard_results/Raw_FASTQ/{sample}.fastq", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/FastQC/{sample}_fastqc.zip", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/FastQC/{sample}_fastqc.html", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/Processed_FASTQ/{sample}.adapt.qc.fastq", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/FastQC/{sample}.adapt.qc_fastqc.zip", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/FastQC/{sample}.adapt.qc_fastqc.html", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/Fasta/{sample}.adapt.qc.collapser.fasta", sample = SAMPLES_2),
';



}



###bam (workflow 1-6)
###matrix (workflow 1-6)
###NGSplot (workflow 1-6)
###bw (workflow 1-6)
###cluster (workflow 1-6)
###DEGs (workflow 1-6)



if (($Analysis_type eq "protein-coding")||($Analysis_type eq "known lncRNA")) {



if ($Core_workflow eq 1) {



print Profile '        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}_thout/accepted_hits.bam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}_thout/align_summary.txt", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}_align_summary.txt", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.sorted.rmdup.bam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.sorted.unique.bam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.sorted.unique.bam.bai", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}_clout/genes.fpkm_tracking", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}_clout/isoforms.fpkm_tracking", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/cufflink/{sample}.genes.fpkm_tracking", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/cufflink_isoform/{sample}.isoforms.fpkm_tracking", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{pair}_diff_out/gene_exp.diff", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{pair}_diff_out/isoform_exp.diff", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.gene.DEG.txt", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.transcript.DEG.txt", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/cuffdiff/{pair}.gene_exp.diff", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/cuffdiff_isoform/{pair}.isoform_exp.diff", pair = PAIRS),
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".transcript.FPKM.txt",
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{sample}.genebody.avgprof.pdf", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{sample}.genebody.heatmap.pdf", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{sample}.genebody.zip", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{sample}.bw", sample = SAMPLES_2),
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.heatmap.R",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.heatmap.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.PCA.R",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.matrix.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.PCA.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.rotation.txt"';



}



if ($Core_workflow eq 2) {



print Profile '        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.sorted.unique.bam", work_folder = WORKING_FOLDER, workflow = WORKFLOW_NAME, sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.sorted.unique.bam.bai", sample = SAMPLES_2),
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.counts.csv",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.counts.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.csv",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.txt",
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.DESeq2.csv", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.DESeq2.txt", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{sample}.genebody.avgprof.pdf", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{sample}.genebody.heatmap.pdf", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{sample}.genebody.zip", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{sample}.bw", sample = SAMPLES_2),
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.heatmap.R",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.txt.1.heatmap.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.PCA.R",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.txt.1.matrix.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.txt.1.PCA.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.txt.1.rotation.txt"';



}



if ($Core_workflow eq 3) {



print Profile '        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".transcript.FPKM.txt",
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.gene.DEG.txt", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.transcript.DEG.txt", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.genome.sorted.bam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.genome.sorted.bam.bai", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{sample}.genebody.avgprof.pdf", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{sample}.genebody.heatmap.pdf", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{sample}.genebody.zip", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{sample}.bw", sample = SAMPLES_2),
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.heatmap.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.matrix.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.PCA.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.rotation.txt"';



}



if ($Core_workflow eq 4) {



print Profile '        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.bam", sample = SAMPLES_2),
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".transcript.FPKM.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.txt",
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.gene.DEG.txt", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.transcript.DEG.txt", pair = PAIRS),
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.heatmap.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.matrix.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.PCA.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.rotation.txt"';



}



if ($Core_workflow eq 5) {



print Profile '        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.bam", sample = SAMPLES_2),
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".transcript.TPM.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.txt",
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.transcript.DEG.txt", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.gene.DEG.txt", pair = PAIRS),
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.clean.txt.1.heatmap.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.clean.txt.1.PCA.pdf"';



}



if ($Core_workflow eq 6) {



print Profile '        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.sorted.unique.bam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.sorted.unique.bam.bai", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.sorted.rmdup.bam", sample = SAMPLES_2),
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".transcript.FPKM.txt",
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.gene.DEG.txt", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.transcript.DEG.txt", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{sample}.genebody.avgprof.pdf", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{sample}.genebody.heatmap.pdf", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{sample}.genebody.zip", sample = SAMPLES_2),
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.heatmap.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.PCA.pdf"';



}



}#if (($Analysis_type eq "protein-coding")||($Analysis_type eq "known lncRNA")) {



if ($Analysis_type eq "novel lncRNA") {



print Profile '        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{sample}_thout/accepted_hits.bam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{sample}_thout/align_summary.txt", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{sample}_align_summary.txt", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{sample}.sorted.bam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{sample}.sorted.rmdup.bam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{sample}.sorted.sam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{sample}.sorted.unique.sam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{sample}.sorted.unique.bam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{sample}.sorted.unique.bam.bai", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/{sample}_clout/transcripts.gtf", sample = SAMPLES_2),
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/assembly_list.txt",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/merged.gtf",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp.combined.gtf",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp.combined.selected.gtf",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA/cuffcmp.combined.selected.result",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA/cuffcmp.combined.selected.noncoding.gtf",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA/denovo_lncRNA.transcript_to_gene.txt",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/all_lncRNA_mRNA.transcript_to_gene.txt",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/all_lncRNA_mRNA.gtf",
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Processed_BAM/{sample}.bam", sample = SAMPLES_2),
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/" + PROJECT_NAME + ".lncRNA_denovo.transcript.TPM.txt",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.txt",
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/" + PROJECT_NAME + ".lncRNA_denovo.{pair}.transcript.DEG.txt", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/" + PROJECT_NAME + ".lncRNA_denovo.{pair}.gene.DEG.txt", pair = PAIRS),
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt.1.heatmap.pdf",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt.1.PCA.pdf"';



}



if ($Analysis_type eq "circRNA") {



if ($Type eq "Paired-End") {



print Profile '        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/joint_mapping/{sample}.Aligned.sortedByCoord.out.bam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/joint_mapping/{sample}.Aligned.sortedByCoord.out.bam.bai", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/joint_mapping/{sample}.Chimeric.out.junction", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/mate1_mapping/{sample}.Chimeric.out.junction", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/mate2_mapping/{sample}.Chimeric.out.junction", sample = SAMPLES_2),
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/CircCoordinates",
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/CircRNACount",
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/CircSkipJunctions",
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/LinearCount",
        expand(WORKING_FOLDER + "/standard_results/CircRNA/CircTest/{pair}/CircCoordinates", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/CircTest/{pair}/CircRNACount", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/CircTest/{pair}/CircSkipJunctions", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/CircTest/{pair}/LinearCount", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/CircTest/{pair}/" + PROJECT_NAME + ".STAR_DCC_circTest.{pair}.circTest.csv", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/CircTest/{pair}/" + PROJECT_NAME + ".STAR_DCC_circTest.{pair}.circTest.txt", pair = PAIRS)';


}



if ($Type eq "Single-End") {



print Profile '        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/single_mapping/{sample}.Aligned.sortedByCoord.out.bam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/single_mapping/{sample}.Aligned.sortedByCoord.out.bam.bai", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/single_mapping/{sample}.Chimeric.out.junction", sample = SAMPLES_2),
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/CircCoordinates",
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/CircRNACount",
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/CircSkipJunctions",
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/LinearCount",
        expand(WORKING_FOLDER + "/standard_results/CircRNA/CircTest/{pair}/CircCoordinates", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/CircTest/{pair}/CircRNACount", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/CircTest/{pair}/CircSkipJunctions", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/CircTest/{pair}/LinearCount", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/CircTest/{pair}/" + PROJECT_NAME + ".STAR_DCC_circTest.{pair}.circTest.csv", pair = PAIRS),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/CircTest/{pair}/" + PROJECT_NAME + ".STAR_DCC_circTest.{pair}.circTest.txt", pair = PAIRS)';



}



}



if (($SNP eq "Yes")||($Fusion eq "Yes")) {print Profile ",\n";}



else {print Profile "\n";}



if (($SNP eq "Yes")&&($Fusion eq "No")) {



print Profile '        expand(WORKING_FOLDER + "/standard_results/SNP/{sample}/vcf/{sample}.filter.vcf", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/SNP/{sample}/1pass/{sample}.Aligned.out.sam", sample = SAMPLES_2)
';



}



if (($SNP eq "No")&&($Fusion eq "Yes")) {



print Profile '        expand(WORKING_FOLDER + "/standard_results/Fusion_gene/Fusion_result/{sample}/summary_candidate_fusions.txt", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/Fusion_gene/Fusion_result/{sample}/final-list_candidate-fusion-genes.txt", sample = SAMPLES_2)
';



}



if (($SNP eq "Yes")&&($Fusion eq "Yes")) {



print Profile '        expand(WORKING_FOLDER + "/standard_results/SNP/{sample}/vcf/{sample}.filter.vcf", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/SNP/{sample}/1pass/{sample}.Aligned.out.sam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/Fusion_gene/Fusion_result/{sample}/summary_candidate_fusions.txt", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/Fusion_gene/Fusion_result/{sample}/final-list_candidate-fusion-genes.txt", sample = SAMPLES_2)
';



}



###############################################Write up all the desired results end############################################################


###############################################Before mapping start############################################################

############################Single with adaptor start#############################



if (($Type eq "Single-End")&&($Adaptor ne "No")) {



print Profile '
rule rename_fastq:
        input:
                expand(WORKING_FOLDER + "/{sample}.fastq", sample = SAMPLES_1)
        output:
                expand(WORKING_FOLDER + "/standard_results/Raw_FASTQ/{sample}.fastq", sample = SAMPLES_2)
        run:
                shell("mkdir -p {WORKING_FOLDER}/standard_results")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Raw_FASTQ")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Processed_FASTQ")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Fasta")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/FastQC")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/SNP")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Snakemake_logs")
                for i in range(len(SAMPLES_1)):
                    single_in = "".join(["{WORKING_FOLDER}/", SAMPLES_1[i], ".fastq"])
                    single_out = "".join(["{WORKING_FOLDER}/standard_results/Raw_FASTQ/", SAMPLES_2[i], ".fastq"])
                    cmd_single = " ".join(["cp", single_in, single_out])
                    shell(cmd_single)



rule fastqc_fisrt:
    input:
        WORKING_FOLDER + "/standard_results/Raw_FASTQ/{smp}.fastq"
    output:
        WORKING_FOLDER + "/standard_results/FastQC/{smp}_fastqc.zip",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}_fastqc.html"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_fisrt.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_fisrt.stderr"
    threads: THREADS
    shell: """
           fastqc -o {FASTQC_OUTPUT_FOLDER} --contaminants {CONTAMINANT_LIST_FILE} --adapters {ADAPTER_LIST_FILE} --threads {threads} {input} 2> {log.err} 1> {log.out}
           """



rule cutadapt:
    input:
        WORKING_FOLDER + "/standard_results/Raw_FASTQ/{smp}.fastq"
    output:
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.fastq"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.cutadapt.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.cutadapt.stderr"
    shell: """
           cutadapt -a {ADAPTER_SEQUENCE} -m 20 -o {output} {input} 2> {log.err} 1> {log.out}
           """



rule fastq_quality_filter:
    input:
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.fastq"
    output:
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastq_quality_filter.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastq_quality_filter.stderr"
    shell: """
           fastq_quality_filter -q 20 -p 90 -Q 33 -i {input} -o {output} 2> {log.err} 1> {log.out}
           """



rule fastqc_second:
    input:
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq"
    output:
        WORKING_FOLDER + "/standard_results/FastQC/{smp}.adapt.qc_fastqc.zip",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}.adapt.qc_fastqc.html"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_second.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_second.stderr"
    threads: THREADS
    shell: """
           fastqc -o {FASTQC_OUTPUT_FOLDER} --contaminants {CONTAMINANT_LIST_FILE} --adapters {ADAPTER_LIST_FILE} --threads {threads} {input} 2> {log.err} 1> {log.out}
           """



rule fastx_collapser:
    input:
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq"
    output:
        WORKING_FOLDER + "/standard_results/Fasta/{smp}.adapt.qc.collapser.fasta"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastx_collapser.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastx_collapser.stderr"
    shell: """
           fastx_collapser -Q 33 -v -i {input} -o {output} 2> {log.err} 1> {log.out}
           """



';



} #if (($Type eq "Single-End")&&($Adaptor ne "No")) {



############################Single with adaptor end##################################

############################Single without adaptor start#############################



if (($Type eq "Single-End")&&($Adaptor eq "No")) {



print Profile '
rule rename_fastq:
        input:
                expand(WORKING_FOLDER + "/{sample}.fastq", sample = SAMPLES_1)
        output:
                expand(WORKING_FOLDER + "/standard_results/Raw_FASTQ/{sample}.fastq", sample = SAMPLES_2)
        run:
                shell("mkdir -p {WORKING_FOLDER}/standard_results")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Raw_FASTQ")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Processed_FASTQ")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Fasta")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/FastQC")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/SNP")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Snakemake_logs")
                for i in range(len(SAMPLES_1)):
                    single_in = "".join(["{WORKING_FOLDER}/", SAMPLES_1[i], ".fastq"])
                    single_out = "".join(["{WORKING_FOLDER}/standard_results/Raw_FASTQ/", SAMPLES_2[i], ".fastq"])
                    cmd_single = " ".join(["cp", single_in, single_out])
                    shell(cmd_single)



rule fastqc_fisrt:
    input:
        WORKING_FOLDER + "/standard_results/Raw_FASTQ/{smp}.fastq"
    output:
        WORKING_FOLDER + "/standard_results/FastQC/{smp}_fastqc.zip",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}_fastqc.html"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_fisrt.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_fisrt.stderr"
    threads: THREADS
    shell: """
           fastqc -o {FASTQC_OUTPUT_FOLDER} --contaminants {CONTAMINANT_LIST_FILE} --adapters {ADAPTER_LIST_FILE} --threads {threads} {input} 2> {log.err} 1> {log.out}
           """



rule fastq_quality_filter:
    input:
        WORKING_FOLDER + "/standard_results/Raw_FASTQ/{smp}.fastq"
    output:
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastq_quality_filter.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastq_quality_filter.stderr"
    shell: """
           fastq_quality_filter -q 20 -p 90 -Q 33 -i {input} -o {output} 2> {log.err} 1> {log.out}
           """



rule fastqc_second:
    input:
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq"
    output:
        WORKING_FOLDER + "/standard_results/FastQC/{smp}.adapt.qc_fastqc.zip",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}.adapt.qc_fastqc.html"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_second.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_second.stderr"
    threads: THREADS
    shell: """
           fastqc -o {FASTQC_OUTPUT_FOLDER} --contaminants {CONTAMINANT_LIST_FILE} --adapters {ADAPTER_LIST_FILE} --threads {threads} {input} 2> {log.err} 1> {log.out}
           """



rule fastx_collapser:
    input:
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq"
    output:
        WORKING_FOLDER + "/standard_results/Fasta/{smp}.adapt.qc.collapser.fasta"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastx_collapser.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastx_collapser.stderr"
    shell: """
           fastx_collapser -Q 33 -v -i {input} -o {output} 2> {log.err} 1> {log.out}
           """



';



}#if (($Type eq "Single-End")&&($Adaptor eq "No")) {



############################Single without adaptor end#############################

############################Paired with adaptor start#############################



if (($Type eq "Paired-End")&&($Adaptor ne "No")) {



print Profile '



rule rename_fastq:
        input:
                R1=expand(WORKING_FOLDER + "/{sample}.R1.fastq", sample = SAMPLES_1), 
                R2=expand(WORKING_FOLDER + "/{sample}.R2.fastq", sample = SAMPLES_1)
        output:
                R1=expand(WORKING_FOLDER + "/standard_results/Raw_FASTQ/{sample}.R1.fastq", sample = SAMPLES_2),
                R2=expand(WORKING_FOLDER + "/standard_results/Raw_FASTQ/{sample}.R2.fastq", sample = SAMPLES_2)
        run:
                shell("mkdir -p {WORKING_FOLDER}/standard_results")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Raw_FASTQ")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Processed_FASTQ")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Fasta")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/FastQC")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/SNP")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Snakemake_logs")
                for i in range(len(SAMPLES_1)):
                    left_in = "".join(["{WORKING_FOLDER}/", SAMPLES_1[i], ".R1.fastq"])
                    right_in = "".join(["{WORKING_FOLDER}/", SAMPLES_1[i], ".R2.fastq"])
                    left_out = "".join(["{WORKING_FOLDER}/standard_results/Raw_FASTQ/", SAMPLES_2[i], ".R1.fastq"])
                    right_out = "".join(["{WORKING_FOLDER}/standard_results/Raw_FASTQ/", SAMPLES_2[i], ".R2.fastq"])
                    cmd_left = " ".join(["cp", left_in, left_out])
                    cmd_right = " ".join(["cp", right_in, right_out])
                    shell(cmd_left)
                    shell(cmd_right)



rule fastqc_fisrt:
    input:
        R1=WORKING_FOLDER + "/standard_results/Raw_FASTQ/{smp}.R1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Raw_FASTQ/{smp}.R2.fastq"
    output:
        WORKING_FOLDER + "/standard_results/FastQC/{smp}.R1_fastqc.zip",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}.R1_fastqc.html",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}.R2_fastqc.zip",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}.R2_fastqc.html"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_fisrt.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_fisrt.stderr"
    threads: THREADS
    shell: """
           fastqc -o {FASTQC_OUTPUT_FOLDER} --contaminants {CONTAMINANT_LIST_FILE} --adapters {ADAPTER_LIST_FILE} --threads {threads} {input.R1} 2> {log.err} 1> {log.out}
           fastqc -o {FASTQC_OUTPUT_FOLDER} --contaminants {CONTAMINANT_LIST_FILE} --adapters {ADAPTER_LIST_FILE} --threads {threads} {input.R2} 2>> {log.err} 1>> {log.out}
    """



rule cutadapt:
    input:
        R1=WORKING_FOLDER + "/standard_results/Raw_FASTQ/{smp}.R1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Raw_FASTQ/{smp}.R2.fastq"
    output:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.R1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.R2.fastq"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.cutadapt.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.cutadapt.stderr"
    shell: """
           cutadapt -a {ADAPTER_SEQUENCE} -A {ADAPTER_SEQUENCE} -m 20 -o {output.R1} -p {output.R2} {input.R1} {input.R2} 2> {log.err} 1> {log.out}
           """



rule prinseq:
    input:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.R1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.R2.fastq"
    output:
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq",
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1_singletons.fastq",
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2_singletons.fastq"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.prinseq.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.prinseq.stderr"
    shell: """
           perl {PRINSEQ_LITE_FILE} -trim_qual_left 20 -trim_qual_right 20 -min_len 20 -min_qual_mean 20 -ns_max_n 2 -out_format 3 -no_qual_header -fastq {input.R1} -fastq2 {input.R2} -out_good {WORKING_FOLDER}/standard_results/Processed_FASTQ/{wildcards.smp} -out_bad null 2> {log.err} 1> {log.out}
           """



rule fastqc_second:
    input:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq"
    output:
        WORKING_FOLDER + "/standard_results/FastQC/{smp}_1_fastqc.zip",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}_1_fastqc.html",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}_2_fastqc.zip",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}_2_fastqc.html"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_second.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_second.stderr"
    threads: THREADS
    shell: """
           fastqc -o {FASTQC_OUTPUT_FOLDER} --contaminants {CONTAMINANT_LIST_FILE} --adapters {ADAPTER_LIST_FILE} --threads {threads} {input.R1} 2> {log.err} 1> {log.out}
           fastqc -o {FASTQC_OUTPUT_FOLDER} --contaminants {CONTAMINANT_LIST_FILE} --adapters {ADAPTER_LIST_FILE} --threads {threads} {input.R2} 2>> {log.err} 1>> {log.out}
           """



rule fastx_collapser:
    input:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq"
    output:
        R1=WORKING_FOLDER + "/standard_results/Fasta/{smp}_1.collapser.fasta",
        R2=WORKING_FOLDER + "/standard_results/Fasta/{smp}_2.collapser.fasta"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastx_collapser.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastx_collapser.stderr"
    shell: """
           fastx_collapser -Q 33 -v -i {input.R1} -o {output.R1} 2> {log.err} 1> {log.out}
           fastx_collapser -Q 33 -v -i {input.R2} -o {output.R2} 2>> {log.err} 1>> {log.out}
           """



';



}#if (($Type eq "Paired-End")&&($Adaptor ne "No")) {



############################Paired with adaptor end#############################

############################Paired without adaptor start#############################



if (($Type eq "Paired-End")&&($Adaptor eq "No")) {



print Profile '
rule rename_fastq:
        input:
                R1=expand(WORKING_FOLDER + "/{sample}.R1.fastq", sample = SAMPLES_1), 
                R2=expand(WORKING_FOLDER + "/{sample}.R2.fastq", sample = SAMPLES_1)
        output:
                R1=expand(WORKING_FOLDER + "/standard_results/Raw_FASTQ/{sample}.R1.fastq", sample = SAMPLES_2),
                R2=expand(WORKING_FOLDER + "/standard_results/Raw_FASTQ/{sample}.R2.fastq", sample = SAMPLES_2)
        run:
                shell("mkdir -p {WORKING_FOLDER}/standard_results")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Raw_FASTQ")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Processed_FASTQ")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Fasta")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/FastQC")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/SNP")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene")
                shell("mkdir -p {WORKING_FOLDER}/standard_results/Snakemake_logs")
                for i in range(len(SAMPLES_1)):
                    left_in = "".join(["{WORKING_FOLDER}/", SAMPLES_1[i], ".R1.fastq"])
                    right_in = "".join(["{WORKING_FOLDER}/", SAMPLES_1[i], ".R2.fastq"])
                    left_out = "".join(["{WORKING_FOLDER}/standard_results/Raw_FASTQ/", SAMPLES_2[i], ".R1.fastq"])
                    right_out = "".join(["{WORKING_FOLDER}/standard_results/Raw_FASTQ/", SAMPLES_2[i], ".R2.fastq"])
                    cmd_left = " ".join(["cp", left_in, left_out])
                    cmd_right = " ".join(["cp", right_in, right_out])
                    shell(cmd_left)
                    shell(cmd_right)



rule fastqc_fisrt:
    input:
        R1=WORKING_FOLDER + "/standard_results/Raw_FASTQ/{smp}.R1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Raw_FASTQ/{smp}.R2.fastq"
    output:
        WORKING_FOLDER + "/standard_results/FastQC/{smp}.R1_fastqc.zip",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}.R1_fastqc.html",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}.R2_fastqc.zip",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}.R2_fastqc.html"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_fisrt.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_fisrt.stderr"
    threads: THREADS
    shell: """
           fastqc -o {FASTQC_OUTPUT_FOLDER} --contaminants {CONTAMINANT_LIST_FILE} --adapters {ADAPTER_LIST_FILE} --threads {threads} {input.R1} 2> {log.err} 1> {log.out}
           fastqc -o {FASTQC_OUTPUT_FOLDER} --contaminants {CONTAMINANT_LIST_FILE} --adapters {ADAPTER_LIST_FILE} --threads {threads} {input.R2} 2>> {log.err} 1>> {log.out}
    """



rule prinseq:
    input:
        R1=WORKING_FOLDER + "/standard_results/Raw_FASTQ/{smp}.R1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Raw_FASTQ/{smp}.R2.fastq"
    output:
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq",
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1_singletons.fastq",
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2_singletons.fastq"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.prinseq.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.prinseq.stderr"
    shell: """
           perl {PRINSEQ_LITE_FILE} -trim_qual_left 20 -trim_qual_right 20 -min_len 20 -min_qual_mean 20 -ns_max_n 2 -out_format 3 -no_qual_header -fastq {input.R1} -fastq2 {input.R2} -out_good {WORKING_FOLDER}/standard_results/Processed_FASTQ/{wildcards.smp} -out_bad null 2> {log.err} 1> {log.out}
           """



rule fastqc_second:
    input:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq"
    output:
        WORKING_FOLDER + "/standard_results/FastQC/{smp}_1_fastqc.zip",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}_1_fastqc.html",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}_2_fastqc.zip",
        WORKING_FOLDER + "/standard_results/FastQC/{smp}_2_fastqc.html"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_second.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastqc_second.stderr"
    threads: THREADS
    shell: """
           fastqc -o {FASTQC_OUTPUT_FOLDER} --contaminants {CONTAMINANT_LIST_FILE} --adapters {ADAPTER_LIST_FILE} --threads {threads} {input.R1} 2> {log.err} 1> {log.out}
           fastqc -o {FASTQC_OUTPUT_FOLDER} --contaminants {CONTAMINANT_LIST_FILE} --adapters {ADAPTER_LIST_FILE} --threads {threads} {input.R2} 2>> {log.err} 1>> {log.out}
           """



rule fastx_collapser:
    input:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq"
    output:
        R1=WORKING_FOLDER + "/standard_results/Fasta/{smp}_1.collapser.fasta",
        R2=WORKING_FOLDER + "/standard_results/Fasta/{smp}_2.collapser.fasta"
    log:
        out=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastx_collapser.stdout",
        err=WORKING_FOLDER + "/standard_results/Snakemake_logs/{smp}.fastx_collapser.stderr"
    shell: """
           fastx_collapser -Q 33 -v -i {input.R1} -o {output.R1} 2> {log.err} 1> {log.out}
           fastx_collapser -Q 33 -v -i {input.R2} -o {output.R2} 2>> {log.err} 1>> {log.out}
           """



';



}#if (($Type eq "Paired-End")&&($Adaptor eq "No")) {

############################Paired without adaptor end#############################

###############################################Before mapping end############################################################



if (($Analysis_type eq "protein-coding")||($Analysis_type eq "known lncRNA")) {



###############################################Core workflow start############################################################

############################Workflow 1: Tophat_Cufflink_Cuffdiff start#################################



if ($Core_workflow eq 1) {



print Profile '
rule bowtie2_build:
    input:
        GENOME_FA_FILE
    output:
        fa=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".fa",
        index_1=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".1.bt2",
        index_2=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".2.bt2",
        index_3=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".3.bt2",
        index_4=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".4.bt2",
        index_5=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".rev.1.bt2",
        index_6=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".rev.2.bt2"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/bowtie2_build.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/bowtie2_build.stderr"
    threads: THREADS
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/NGS_plot
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Bigwig
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Snakemake_logs
           cp {input} {output.fa} 2> {log.err} 1> {log.out}
           bowtie2-build --threads {threads} {output.fa} {TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE} 2>> {log.err} 1>> {log.out}
           """



';



if ($Type eq "Paired-End") {



print Profile '
rule tophat2:
    input:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq",
        index_1=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".1.bt2",
        index_2=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".2.bt2",
        index_3=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".3.bt2",
        index_4=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".4.bt2",
        index_5=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".rev.1.bt2",
        index_6=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".rev.2.bt2"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}_thout/accepted_hits.bam",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}_thout/align_summary.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.tophat2.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.tophat2.stderr"
    threads: THREADS
    shell: """
           tophat2 -p {threads} --mate-inner-dist {TOPHAT_MATE_INNER_DIST} --mate-std-dev {TOPHAT_MATE_STD_DEV} --library-type {TOPHAT_CUFFLINKS_CUFFDIFF_LIBRARY_TYPE} -G {GTF} -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/{wildcards.smp}_thout --no-novel-juncs {TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE} {input.R1} {input.R2} 2> {log.err} 1> {log.out}
           """



';



}



if ($Type eq "Single-End") {



print Profile '
rule tophat2:
    input:
        fastq=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq",
        index_1=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".1.bt2",
        index_2=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".2.bt2",
        index_3=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".3.bt2",
        index_4=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".4.bt2",
        index_5=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".rev.1.bt2",
        index_6=TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE + ".rev.2.bt2"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}_thout/accepted_hits.bam",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}_thout/align_summary.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.tophat2.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.tophat2.stderr"
    threads: THREADS
    shell: """
           tophat2 -p {threads} --library-type {TOPHAT_CUFFLINKS_CUFFDIFF_LIBRARY_TYPE} -G {GTF} -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/{wildcards.smp}_thout --no-novel-juncs {TOPHAT_CUFFLINK_CUFFDIFF_INDEX_FILE} {input.fastq} 2> {log.err} 1> {log.out}
           """



';



}



print Profile '
rule copy_summary:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}_thout/align_summary.txt"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}_align_summary.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.copy_summary.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.copy_summary.stderr"
    shell: """
           cp {input} {output} 2> {log.err} 1> {log.out}
           """



rule sort_bam:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}_thout/accepted_hits.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sort_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sort_bam.stderr"
    shell: """
           samtools sort -o {output} {input} 2> {log.err} 1> {log.out}
           """



';



if ($Type eq "Paired-End") {



print Profile '
rule rmdup_bam:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.rmdup.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.rmdup_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.rmdup_bam.stderr"
    shell: """
           samtools rmdup {input} {output} 2> {log.err} 1> {log.out}
           """



';



}



if ($Type eq "Single-End") {



print Profile '
rule rmdup_bam:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.rmdup.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.rmdup_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.rmdup_bam.stderr"
    shell: """
           samtools rmdup -s {input} {output} 2> {log.err} 1> {log.out}
           """



';



}



print Profile '
rule bam2sam:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bam2sam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bam2sam.stderr"
    shell: """
           samtools view -h {input} > {output} 2> {log.err}
           """



rule select_unique_mapping:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.sam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.select_unique_mapping.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.select_unique_mapping.stderr"
    shell: """
           perl {EXTRACT_UNIQUE} --input_sam {input} --output {output} --aligner T 2> {log.err} 1> {log.out}
           """



rule sam2bam:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.sam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sam2bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sam2bam.stderr"
    shell: """
           samtools view -S -b -o {output} {input} 2> {log.err} 1> {log.out}
           """



rule index_bam:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam.bai"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.index_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.index_bam.stderr"
    shell: """
           samtools index {input} 2> {log.err} 1> {log.out}
           """



rule cufflinks:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}_clout/genes.fpkm_tracking",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}_clout/isoforms.fpkm_tracking"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.cufflinks.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.cufflinks.stderr"
    threads: THREADS
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/cufflink
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/cufflink_isoform
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/cuffdiff
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/cuffdiff_isoform
           cufflinks -p {threads} --library-type {TOPHAT_CUFFLINKS_CUFFDIFF_LIBRARY_TYPE} -G {GTF} -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/{wildcards.smp}_clout {input} 2> {log.err} 1> {log.out}
           """



rule copy_fpkm:
    input:
        genes=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}_clout/genes.fpkm_tracking",
        isoforms=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}_clout/isoforms.fpkm_tracking"
    output:
        genes=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/cufflink/{smp}.genes.fpkm_tracking",
        isoforms=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/cufflink_isoform/{smp}.isoforms.fpkm_tracking"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.copy_fpkm.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.copy_fpkm.stderr"
    shell: """
           cp {input.genes} {output.genes} 2> {log.err} 1> {log.out}
           cp {input.isoforms} {output.isoforms} 2>> {log.err} 1>> {log.out}
           """



';



###DEGs



print Profile "\n\n\n";



###DEGs



for (my $i=0;$i<@treatment;$i++) {



#$control[$i]
#$treatment[$i]



print Profile "rule cuffdiff_".$control[$i]."_vs_".$treatment[$i].":\n";
print Profile "    input:\n";
print Profile "        expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Processed_BAM/{sample}.sorted.unique.bam\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i].")\n";
print Profile "    output:\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Processed_BAM/".$control[$i]."_vs_".$treatment[$i]."_diff_out/gene_exp.diff\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Processed_BAM/".$control[$i]."_vs_".$treatment[$i]."_diff_out/isoform_exp.diff\"\n";
print Profile "    log:\n";
print Profile "        out=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/cuffdiff_".$control[$i]."_vs_".$treatment[$i].".stdout\",\n";
print Profile "        err=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/cuffdiff_".$control[$i]."_vs_".$treatment[$i].".stderr\"\n";
print Profile "    threads: THREADS\n";
print Profile "    run:\n";
print Profile "        classes = [\",\".join(expand(\"{work_folder}/standard_results/{workflow}/Processed_BAM/{sample}.sorted.unique.bam\", work_folder = WORKING_FOLDER, workflow = WORKFLOW_NAME, sample=cls)) for cls in (CLASS_".$control[$i].", CLASS_".$treatment[$i].")]\n";
print Profile "        shell(\"cuffdiff -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/".$control[$i]."_vs_".$treatment[$i]."_diff_out -p {THREADS} --library-type {TOPHAT_CUFFLINKS_CUFFDIFF_LIBRARY_TYPE} -L ".$control[$i].",".$treatment[$i]." {GTF} {classes[0]} {classes[1]} 2> {log.err} 1> {log.out}\")\n\n\n";



}###for (my $i=0;$i<@treatment;$i++) {



print Profile '
rule copy_cuffdiff:
    input:
        gene=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{pair}_diff_out/gene_exp.diff",
        isoform=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{pair}_diff_out/isoform_exp.diff"
    output:
        gene_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.gene.DEG.txt",
        transcript_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/DEG/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".{pair}.transcript.DEG.txt",
        gene_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/cuffdiff/{pair}.gene_exp.diff",
        transcript_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/cuffdiff_isoform/{pair}.isoform_exp.diff"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{pair}.copy_cuffdiff.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{pair}.copy_cuffdiff.stderr"
    shell: """
           cp {input.gene} {output.gene_1} 2> {log.err} 1> {log.out}
           cp {input.isoform} {output.transcript_1} 2>> {log.err} 1>> {log.out}
           cp {input.gene} {output.gene_2} 2>> {log.err} 1>> {log.out}
           cp {input.isoform} {output.transcript_2} 2>> {log.err} 1>> {log.out}
           """



rule perl_gene_level:
    input:
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/cufflink/{sample}.genes.fpkm_tracking", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/cuffdiff/{pair}.gene_exp.diff", pair = PAIRS)
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/perl_gene_level.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/perl_gene_level.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/combine_cufflinks_results_new.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/cufflink -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cufflink.combine 2> {log.err} 1> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/combine_cuffdiff_results_new.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/cuffdiff -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff.combine 2>> {log.err} 1>> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/combine_cuffdiff_cufflink_new.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff.combine -j {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cufflink.combine -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink.combine 2>> {log.err} 1>> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/classify_result_new_new.pl -i {GTF} -j {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink.combine -o1 {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink.combine.class -o2 {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink.combine.signature 2>> {log.err} 1>> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/only_pc.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink.combine.class -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink.combine.class.pc 2>> {log.err} 1>> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/only_FPKM.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink.combine.class.pc -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.{WORKFLOW_NAME}.gene.FPKM.txt 2>> {log.err} 1>> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/clean_FPKM.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.{WORKFLOW_NAME}.gene.FPKM.txt -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster/{PROJECT_NAME}.{WORKFLOW_NAME}.gene.FPKM.clean.txt 2>> {log.err} 1>> {log.out}
           """



rule perl_isoform_level:
    input:
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/cufflink_isoform/{sample}.isoforms.fpkm_tracking", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/cuffdiff_isoform/{pair}.isoform_exp.diff", pair = PAIRS)
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".transcript.FPKM.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/perl_isoform_level.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/perl_isoform_level.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/combine_cufflinks_isoform_results_new.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/cufflink_isoform -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cufflink_isoform.combine 2> {log.err} 1> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/combine_cuffdiff_results_new.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/cuffdiff_isoform -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_isoform.combine 2>> {log.err} 1>> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/combine_cuffdiff_cufflink_isoform_new.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_isoform.combine -j {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cufflink_isoform.combine -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink_isoform.combine 2>> {log.err} 1>> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/classify_result_new_new_isoform.pl -i {GTF} -j {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink_isoform.combine -o1 {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink_isoform.combine.class -o2 {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink_isoform.combine.signature 2>> {log.err} 1>> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/only_pc.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink_isoform.combine.class -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink_isoform.combine.class.pc 2>> {log.err} 1>> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/only_FPKM_isoform.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.cuffdiff_cufflink_isoform.combine.class.pc -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/{PROJECT_NAME}.{WORKFLOW_NAME}.transcript.FPKM.txt 2>> {log.err} 1>> {log.out}
           """



rule ngs_plot:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{smp}.genebody.avgprof.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{smp}.genebody.heatmap.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{smp}.genebody.zip"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.ngs_plot.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.ngs_plot.stderr"
    threads: THREADS
    shell: """
           ngs.plot.r -G {ASSEMBLY_VERSION} -R genebody -C {input} -O {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/NGS_plot/NGSplot.{wildcards.smp}.genebody -T {wildcards.smp} -L 3000 -RB 0.05 -F rnaseq 2> {log.err} 1> {log.out}
           """



rule bam2wig:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.gz"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bam2wig.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bam2wig.stderr"
    shell: """
           samtools mpileup -BQ0 {input} | perl -ne \'BEGIN{{print "track type=wiggle_0 name={wildcards.smp} maxHeightPixels=64:64:11 color=31,120,180 visibility=full\\\\n"}};($c, $start, undef, $depth) = split; if ($c ne $lastC) {{ print "variableStep chrom=$c\\\\n"; }};$lastC=$c;next unless $. % 5 ==0;print "$start\\\\t$depth\\\\n" unless $depth<3;\' | gzip -c > {output} 2> {log.err}
           """



rule wig2wig2:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.gz"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.2.gz"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig2wig2.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig2wig2.stderr"
    shell: """
           zcat {input} | sed \'1d\' | gzip -c > {output} 2> {log.err}
           """



rule wig22bw:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.2.gz"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.bw"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig22bw.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig22bw.stderr"
    shell: """
           {WIGTOBIGWIG} {input} {CHROM_SIZES} {output} 2> {log.err} 1> {log.out}
           """



rule heatmap:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt"
    output:
        R=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.heatmap.R",
        pdf=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.heatmap.pdf"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/heatmap.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/heatmap.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_heatmap_Rscript.pl -i {PROJECT_NAME}.{WORKFLOW_NAME}.gene.FPKM.clean.txt -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



rule pca:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt"
    output:
        R=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.PCA.R",
        matrix=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.matrix.txt",
        pdf=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.PCA.pdf",
        rotation=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.rotation.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/pca.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/pca.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_PCA_Rscript.pl -i {PROJECT_NAME}.{WORKFLOW_NAME}.gene.FPKM.clean.txt -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



';



}



############################Workflow 1: Tophat_Cufflink_Cuffdiff end#################################

############################Workflow 2: Subread_featureCounts_DESeq2 strat#############################



if ($Core_workflow eq 2) {



print Profile '
rule build_index:
    input:
        HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + "/genome.fa",
        HPPRNA_SOFTWARE_FOLDER + "/my_R/" + WORKFLOW_NAME + ".buildindex.R"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".00.b.array",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".00.b.tab",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".reads"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/build_index.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/build_index.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/NGS_plot
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Bigwig
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Snakemake_logs
           R --slave --vanilla --args "{WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION}" "{HPPRNA_SOFTWARE_FOLDER}/hppRNA_genome/{ASSEMBLY_VERSION}/genome.fa" < {HPPRNA_SOFTWARE_FOLDER}/my_R/{WORKFLOW_NAME}.buildindex.R 2> {log.err} 1> {log.out}
           """



';



if ($Type eq "Paired-End") {



print Profile '
rule align:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".00.b.array",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".00.b.tab",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".reads",
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq",
        HPPRNA_SOFTWARE_FOLDER + "/my_R/" + WORKFLOW_NAME + ".align.R"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.unique.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.align.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.align.stderr"
    shell: """
           R --slave --vanilla --args Paired-End {wildcards.smp} {WORKING_FOLDER}/standard_results {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION} {SUBREAD_MINFRAGLENGTH} {SUBREAD_MAXFRAGLENGTH} {SUBREAD_PE_ORIENTATION} < {HPPRNA_SOFTWARE_FOLDER}/my_R/{WORKFLOW_NAME}.align.R 2> {log.err} 1> {log.out}
           """



';



}



if ($Type eq "Single-End") {



print Profile '
rule align:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".00.b.array",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".00.b.tab",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".reads",
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.unique.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.align.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.align.stderr"
    shell: """
           R --slave --vanilla --args Single-End {wildcards.smp} {WORKING_FOLDER}/standard_results {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION} < {HPPRNA_SOFTWARE_FOLDER}/my_R/{WORKFLOW_NAME}.align.R 2> {log.err} 1> {log.out}
           """



';



}



print Profile '
rule sort_bam:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.unique.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sort_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sort_bam.stderr"
    shell: """
           samtools sort -o {output} {input} 2> {log.err} 1> {log.out}
           """



rule index_bam:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam.bai"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.index_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.index_bam.stderr"
    shell: """
           samtools index {input} 2> {log.err} 1> {log.out}
           """



';



if ($Type eq "Paired-End") {



print Profile '
rule feature_counts:
    input:
        GTF,
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.sorted.unique.bam", sample = SAMPLES_2)
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.counts.csv",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.counts.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.csv",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/feature_counts.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/feature_counts.stderr"
    threads: THREADS
    shell: """
           R --slave --vanilla --args Paired-End {WORKING_FOLDER}/standard_results {GTF} "{SAMPLE_join}" {PROJECT_NAME} {threads} {FEATURECOUNTS_STRANDSPECIFIC} {SUBREAD_MINFRAGLENGTH} {SUBREAD_MAXFRAGLENGTH} < {HPPRNA_SOFTWARE_FOLDER}/my_R/{WORKFLOW_NAME}.featureCounts.R 2> {log.err} 1> {log.out}
           """



';



}



if ($Type eq "Single-End") {



print Profile '
rule feature_counts:
    input:
        GTF,
        expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{sample}.sorted.unique.bam", sample = SAMPLES_2)
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.counts.csv",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.counts.txt",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.csv",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/feature_counts.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/feature_counts.stderr"
    threads: THREADS
    shell: """
           R --slave --vanilla --args Single-End {WORKING_FOLDER}/standard_results {GTF} "{SAMPLE_join}" {PROJECT_NAME} {threads} {FEATURECOUNTS_STRANDSPECIFIC} < {HPPRNA_SOFTWARE_FOLDER}/my_R/{WORKFLOW_NAME}.featureCounts.R 2> {log.err} 1> {log.out}
           """



';



}



###DEGs



#CLASS_testis_2 = "R_testis_7a;R_testis_7b;R_testis_7c"

#CLASS_brain_2 = "R_brain_3b;R_brain_3c;R_brain_a"



#$control[$i]
#$treatment[$i]



if ($Type eq "Paired-End") {



for (my $i=0;$i<@treatment;$i++) {


print Profile "\n\n\n";

print Profile "rule deseq2_".$control[$i]."_vs_".$treatment[$i].":\n";
print Profile "    input:\n";
print Profile "        expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Processed_BAM/{sample}.sorted.unique.bam\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "        GTF\n";
print Profile "    output:\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".col.csv\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".DESeq2.csv\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".DESeq2.txt\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".gene.counts.csv\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".gene.counts.txt\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".gene.RPKM.csv\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".gene.RPKM.txt\"\n";
print Profile "    log:\n";
print Profile "        out=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/deseq2_".$control[$i]."_vs_".$treatment[$i].".stdout\",\n";
print Profile "        err=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/deseq2_".$control[$i]."_vs_".$treatment[$i].".stderr\"\n";
print Profile "    threads: THREADS\n";
print Profile "    shell: \"\"\"\n";
print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_DESeq2_configure.pl -c1 ".$control[$i]." -c2 ".$treatment[$i]." -b1 \"{CLASS_".$control[$i]."_2}\" -b2 \"{CLASS_".$treatment[$i]."_2}\" -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/{PROJECT_NAME}.{WORKFLOW_NAME}.".$control[$i]."_vs_".$treatment[$i].".col.csv 2> {log.err} 1> {log.out}\n";
print Profile "           R --slave --vanilla --args Paired-End ".$control[$i]." ".$treatment[$i]." \"{CLASS_".$control[$i]."_2}\" \"{CLASS_".$treatment[$i]."_2}\" {GTF} {WORKING_FOLDER}/standard_results {PROJECT_NAME} {threads} {FEATURECOUNTS_STRANDSPECIFIC} {SUBREAD_MINFRAGLENGTH} {SUBREAD_MAXFRAGLENGTH} < {HPPRNA_SOFTWARE_FOLDER}/my_R/{WORKFLOW_NAME}.DESeq2.R 2>> {log.err} 1>> {log.out}\n";
print Profile "           \"\"\"\n\n\n";



}#for (my $i=0;$i<@treatment;$i++) {



}#if ($Type eq "Paired-End") {



if ($Type eq "Single-End") {



#$control[$i]
#$treatment[$i]



for (my $i=0;$i<@treatment;$i++) {



print Profile "\n\n\n";



print Profile "rule deseq2_".$control[$i]."_vs_".$treatment[$i].":\n";
print Profile "    input:\n";
print Profile "        expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Processed_BAM/{sample}.sorted.unique.bam\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "        GTF\n";
print Profile "    output:\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".col.csv\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".DESeq2.csv\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".DESeq2.txt\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".gene.counts.csv\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".gene.counts.txt\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".gene.RPKM.csv\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".gene.RPKM.txt\"\n";
print Profile "    log:\n";
print Profile "        out=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/deseq2_".$control[$i]."_vs_".$treatment[$i].".stdout\",\n";
print Profile "        err=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/deseq2_".$control[$i]."_vs_".$treatment[$i].".stderr\"\n";
print Profile "    threads: THREADS\n";
print Profile "    shell: \"\"\"\n";
print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_DESeq2_configure.pl -c1 ".$control[$i]." -c2 ".$treatment[$i]." -b1 \"{CLASS_".$control[$i]."_2}\" -b2 \"{CLASS_".$treatment[$i]."_2}\" -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/{PROJECT_NAME}.{WORKFLOW_NAME}.".$control[$i]."_vs_".$treatment[$i].".col.csv 2> {log.err} 1> {log.out}\n";
print Profile "           R --slave --vanilla --args Single-End ".$control[$i]." ".$treatment[$i]." \"{CLASS_".$control[$i]."_2}\" \"{CLASS_".$treatment[$i]."_2}\" {GTF} {WORKING_FOLDER}/standard_results {PROJECT_NAME} {threads} {FEATURECOUNTS_STRANDSPECIFIC} < {HPPRNA_SOFTWARE_FOLDER}/my_R/{WORKFLOW_NAME}.DESeq2.R 2>> {log.err} 1>> {log.out}\n";
print Profile "           \"\"\"\n\n\n";



}#for (my $i=0;$i<@treatment;$i++) {



}#if ($Type eq "Single-End") {



print Profile '
rule ngs_plot:
    input:
        bam=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam",
        bai=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam.bai"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{smp}.genebody.avgprof.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{smp}.genebody.heatmap.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{smp}.genebody.zip"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.ngs_plot.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.ngs_plot.stderr"
    shell: """
           ngs.plot.r -G {ASSEMBLY_VERSION} -R genebody -C {input.bam} -O {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/NGS_plot/NGSplot.{wildcards.smp}.genebody -T {wildcards.smp} -L 3000 -RB 0.05 -F rnaseq 2> {log.err} 1> {log.out}
           """



rule bam2wig:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.gz"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bam2wig.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bam2wig.stderr"
    shell: """
           samtools mpileup -BQ0 {input} | perl -ne \'BEGIN{{print "track type=wiggle_0 name={wildcards.smp} maxHeightPixels=64:64:11 color=31,120,180 visibility=full\\\\n"}};($c, $start, undef, $depth) = split; if ($c ne $lastC) {{ print "variableStep chrom=$c\\\\n"; }};$lastC=$c;next unless $. % 5 ==0;print "$start\\\\t$depth\\\\n" unless $depth<3;\' | gzip -c > {output} 2> {log.err}
           """



rule wig2wig2:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.gz"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.2.gz"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig2wig2.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig2wig2.stderr"
    shell: """
           zcat {input} | sed \'1d\' | gzip -c > {output} 2> {log.err}
           """



rule wig22bw:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.2.gz"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.bw"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig22bw.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig22bw.stderr"
    shell: """
           {WIGTOBIGWIG} {input} {CHROM_SIZES} {output} 2> {log.err} 1> {log.out}
           """



rule prepare_rpkm:
    input:
        perl=HPPRNA_SOFTWARE_FOLDER + "/my_perl/clean_RPKM_" + WORKFLOW_NAME + ".pl",
        txt=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.txt"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/prepare_rpkm.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/prepare_rpkm.stderr"
    shell: """
           perl {input.perl} -i {input.txt} -o {output} 2> {log.err} 1> {log.out}
           """



rule heatmap:
    input:
        HPPRNA_SOFTWARE_FOLDER + "/my_perl/generate_heatmap_Rscript.pl",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.txt"
    output:
        R=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.heatmap.R",
        pdf=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.txt.1.heatmap.pdf"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/heatmap.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/heatmap.stderr"
    shell: """
         perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_heatmap_Rscript.pl -i {PROJECT_NAME}.{WORKFLOW_NAME}.gene.RPKM.clean.txt -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster -o {output.R} 2> {log.err} 1> {log.out}
         R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



rule pca:
    input:
        txt=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.txt",
        perl=HPPRNA_SOFTWARE_FOLDER + "/my_perl/generate_PCA_Rscript.pl"
    output:
        R=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.PCA.R",
        matrix=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.txt.1.matrix.txt",
        pdf=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.txt.1.PCA.pdf",
        rotation=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.RPKM.clean.txt.1.rotation.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/pca.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/pca.stderr"
    shell: """
           perl {input.perl} -i {PROJECT_NAME}.{WORKFLOW_NAME}.gene.RPKM.clean.txt -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



';



}



############################Workflow 2: Subread_featureCounts_DESeq2 end#############################


############################Workflow 3: STAR_RSEM_EBSeq start##########################################



if ($Core_workflow eq 3) {



print Profile '
rule rsem_prepare_reference:
    input:
          GTF,
          GENOME_FA_FILE
    output:
          WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/SA",
          WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/SAindex"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/rsem_prepare_reference.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/rsem_prepare_reference.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/NGS_plot
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Bigwig
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Snakemake_logs
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Summary_gene
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Summary_transcript
           rsem-prepare-reference --gtf {GTF} --star -p 1 {GENOME_FA_FILE} {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION} 2> {log.err} 1> {log.out}
           """



';



if ($Type eq "Paired-End") {



print Profile '
rule rsem_calculate_expression:
    input:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq",
        index_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/SA",
        index_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/SAindex"
    output:
        bam=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.genome.bam",
        gene_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.genes.results",
        isoform_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.isoforms.results",
        gene_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Summary_gene/{smp}.genes.results",
        isoform_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Summary_transcript/{smp}.isoforms.results"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.rsem_calculate_expression.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.rsem_calculate_expression.stderr"
    threads: THREADS
    shell: """
           rsem-calculate-expression -p {threads} --forward-prob {RSEM_FORWARD_PROB} --paired-end --star --estimate-rspd --fragment-length-min {RSEM_FRAGMENT_LENGTH_MIN} --fragment-length-max {RSEM_FRAGMENT_LENGTH_MAX} --output-genome-bam {input.R1} {input.R2} {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION} {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/{wildcards.smp} 2> {log.err} 1> {log.out}
           cp {output.gene_1} {output.gene_2} 2>> {log.err} 1>> {log.out}
           cp {output.isoform_1} {output.isoform_2} 2>> {log.err} 1>> {log.out}
           """



';



}



if ($Type eq "Single-End") {



print Profile '
rule rsem_calculate_expression:
    input:
        fastq=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq",
        index_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/SA",
        index_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/SAindex"
    output:
        bam=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.genome.bam",
        gene_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.genes.results",
        isoform_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.isoforms.results",
        gene_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Summary_gene/{smp}.genes.results",
        isoform_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Summary_transcript/{smp}.isoforms.results"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.rsem_calculate_expression.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.rsem_calculate_expression.stderr"
    threads: THREADS
    shell: """
           rsem-calculate-expression -p {threads} --forward-prob {RSEM_FORWARD_PROB} --star --estimate-rspd --append-names --output-genome-bam --fragment-length-mean {RSEM_FRAGMENT_LENGTH_MEAN} --fragment-length-sd {RSEM_FRAGMENT_LENGTH_SD} {input.fastq} {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION} {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/{wildcards.smp} 2> {log.err} 1> {log.out}
           cp {output.gene_1} {output.gene_2} 2>> {log.err} 1>> {log.out}
           cp {output.isoform_1} {output.isoform_2} 2>> {log.err} 1>> {log.out}
           """



';



}



print Profile '
rule fpkm_matrix:
    input:
         perl_gene=HPPRNA_SOFTWARE_FOLDER + "/my_perl/combine_RSEM_gene_results.pl",
         perl_transcript=HPPRNA_SOFTWARE_FOLDER + "/my_perl/combine_RSEM_transcript_results.pl",
         gene=expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Summary_gene/{sample}.genes.results", sample = SAMPLES_2),
         isoform=expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Summary_transcript/{sample}.isoforms.results", sample = SAMPLES_2)
    output:
         gene=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.txt",
         transcript=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".transcript.FPKM.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/fpkm_matrix.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/fpkm_matrix.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/combine_RSEM_gene_results.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Summary_gene -o {output.gene} 2> {log.err} 1> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/combine_RSEM_transcript_results.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Summary_transcript -o {output.transcript} 2>> {log.err} 1>> {log.out}
           """



';



print Profile "\n\n\n";



###DEGs



for (my $i=0;$i<@treatment;$i++) {



my $EBSeq_gene_level="";
my $EBSeq_transcript_level="";
my $treatment_count=0;
my $control_count=0;



for (my $j=0;$j<@{$replicate{$treatment[$i]}};$j++) {

$EBSeq_gene_level.="{WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/".$replicate{$treatment[$i]}[$j].".genes.results ";
$EBSeq_transcript_level.="{WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/".$replicate{$treatment[$i]}[$j].".isoforms.results ";
$treatment_count++;

}


for (my $j=0;$j<@{$replicate{$control[$i]}};$j++) {


$EBSeq_gene_level.="{WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/".$replicate{$control[$i]}[$j].".genes.results ";
$EBSeq_transcript_level.="{WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM/".$replicate{$control[$i]}[$j].".isoforms.results ";
$control_count++;

}



print Profile "rule deg_".$control[$i]."_vs_".$treatment[$i].":\n";
print Profile "    input:\n";
print Profile "         expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Processed_BAM/{sample}.genes.results\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i].")\n";
print Profile "    output:\n";
print Profile "         WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".GeneMat.txt\",\n";
print Profile "         WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".gene.DEG.txt\",\n";
print Profile "         WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".transcript.DEG.txt\"\n";



print Profile "    log:\n";
print Profile "         out=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/deg_".$control[$i]."_vs_".$treatment[$i].".stdout\",\n";
print Profile "         err=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/deg_".$control[$i]."_vs_".$treatment[$i].".stderr\"\n";



print Profile "    shell: \"\"\"\n";
print Profile "           rsem-generate-data-matrix ".$EBSeq_gene_level."> {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/{PROJECT_NAME}.{WORKFLOW_NAME}.".$control[$i]."_vs_".$treatment[$i].".GeneMat.txt\n";
print Profile "           rsem-run-ebseq {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/{PROJECT_NAME}.{WORKFLOW_NAME}.".$control[$i]."_vs_".$treatment[$i].".GeneMat.txt ".$treatment_count.",".$control_count." {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/{PROJECT_NAME}.{WORKFLOW_NAME}.".$control[$i]."_vs_".$treatment[$i].".gene.DEG.txt\n";



#rsem-control-fdr GeneMat.results 0.05 GeneMat.de.txt



print Profile "           rsem-generate-ngvector {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION}.transcripts.fa {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION}\n";
print Profile "           rsem-generate-data-matrix ".$EBSeq_transcript_level."> {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/{PROJECT_NAME}.{WORKFLOW_NAME}.".$control[$i]."_vs_".$treatment[$i].".IsoMat.txt\n";
print Profile "           rsem-run-ebseq --ngvector {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION}.ngvec {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/{PROJECT_NAME}.{WORKFLOW_NAME}.".$control[$i]."_vs_".$treatment[$i].".IsoMat.txt 3,3 {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/{PROJECT_NAME}.{WORKFLOW_NAME}.".$control[$i]."_vs_".$treatment[$i].".transcript.DEG.txt\n";
print Profile "           \"\"\"\n";



#rsem-control-fdr IsoMat.results 0.05 IsoMat.de.txt



} # end of for (my $i=0;$i<@treatment;$i++) {



print Profile "\n\n\n";



print Profile '
rule sort_bam:
    input:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.genome.bam"
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.genome.sorted.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sort_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sort_bam.stderr"
    shell: """
           samtools sort -o {output} {input} 2> {log.err} 1> {log.out}
           """



rule index_bam:
    input:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.genome.sorted.bam"
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.genome.sorted.bam.bai"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.index_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.index_bam.stderr"
    shell: """
           samtools index {input} 2> {log.err} 1> {log.out}
           """



rule ngs_plot:
    input:
         bam=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.genome.sorted.bam",
         bai=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.genome.sorted.bam.bai"
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{smp}.genebody.avgprof.pdf",
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{smp}.genebody.heatmap.pdf",
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{smp}.genebody.zip"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.ngs_plot.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.ngs_plot.stderr"
    shell: """
           ngs.plot.r -G {ASSEMBLY_VERSION} -R genebody -C {input.bam} -O {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/NGS_plot/NGSplot.{wildcards.smp}.genebody -T {wildcards.smp} -L 3000 -RB 0.05 -F rnaseq 2> {log.err} 1> {log.out}
           """



rule bam2wig:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.genome.sorted.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.gz"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bam2wig.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bam2wig.stderr"
    shell: """
           samtools mpileup -BQ0 {input} | perl -ne \'BEGIN{{print "track type=wiggle_0 name={wildcards.smp} maxHeightPixels=64:64:11 color=31,120,180 visibility=full\\\\n"}};($c, $start, undef, $depth) = split; if ($c ne $lastC) {{ print "variableStep chrom=$c\\\\n"; }};$lastC=$c;next unless $. % 5 ==0;print "$start\\\\t$depth\\\\n" unless $depth<3;\' | gzip -c > {output}
           """



rule wig2wig2:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.gz"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.2.gz"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig2wig2.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig2wig2.stderr"
    shell: """
           zcat {input} | sed \'1d\' | gzip -c > {output}
           """



rule wig22bw:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.2.gz"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.bw"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig22bw.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig22bw.stderr"
    shell: """
           {WIGTOBIGWIG} {input} {CHROM_SIZES} {output} 2> {log.err} 1> {log.out}
           """



rule prepare_fpkm:
    input:
         perl=HPPRNA_SOFTWARE_FOLDER + "/my_perl/clean_FPKM_" + WORKFLOW_NAME + ".pl",
         txt=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.txt"
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/prepare_fpkm.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/prepare_fpkm.stderr"
    shell: """
           perl {input.perl} -i {input.txt} -o {output} 2> {log.err} 1> {log.out}
           """



rule heatmap:
    input:
         HPPRNA_SOFTWARE_FOLDER + "/my_perl/generate_heatmap_Rscript.pl",
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt"
    output:
         R=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.heatmap.R",
         pdf=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.heatmap.pdf"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/heatmap.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/heatmap.stderr"
    shell: """
         perl /data/dops-tree-genome/dops0659/hppRNA_software/my_perl/generate_heatmap_Rscript.pl -i {PROJECT_NAME}.{WORKFLOW_NAME}.gene.FPKM.clean.txt -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster -o {output.R} 2> {log.err} 1> {log.out}
         R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



rule pca:
    input:
         txt=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt",
         perl=HPPRNA_SOFTWARE_FOLDER + "/my_perl/generate_PCA_Rscript.pl"
    output:
         R=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.PCA.R",
         matrix=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.matrix.txt",
         pdf=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.PCA.pdf",
         rotation=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.rotation.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/pca.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/pca.stderr"
    shell: """
           perl {input.perl} -i {PROJECT_NAME}.{WORKFLOW_NAME}.gene.FPKM.clean.txt -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



';



}#if ($Core_workflow eq 3) {



############################Workflow 3: STAR_RSEM_EBSeq end##########################################

############################Workflow 4: Bowtie_eXpress_edgeR start#####################################



if ($Core_workflow eq 4) {



print Profile '
rule bowtie2_build:
    input:
         HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + /genome.fa",
         GTF
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".1.bt2",
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".2.bt2",
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".3.bt2",
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".4.bt2",
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".rev.1.bt2",
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".rev.2.bt2",
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".transcript.fa"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/bowtie2_build.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/bowtie2_build.stderr"
    threads: THREADS
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/NGS_plot
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Bigwig
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Snakemake_logs
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Quantification
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Summary
           gffread -w {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION}.transcript.fa -g {HPPRNA_SOFTWARE_FOLDER}/hppRNA_genome/{ASSEMBLY_VERSION}/genome.fa {GTF} 2> {log.err} 1> {log.out}
           bowtie2-build --threads {threads} --offrate 1 {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION}.transcript.fa {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION} 2>> {log.err} 1>> {log.out}
           """



';



if ($Type eq "Paired-End") {



print Profile '
rule bowtie2:
    input:
         R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
         R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq",
         index_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".1.bt2",
         index_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".2.bt2",
         index_3=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".3.bt2",
         index_4=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".4.bt2",
         index_5=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".rev.1.bt2",
         index_6=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".rev.2.bt2",
         index_7=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".transcript.fa"
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bowtie2.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bowtie2.stderr"
    threads: THREADS
    shell: """
           bowtie2 -p {threads} -a -I {BOWTIE_HISAT_MININS} -X {BOWTIE_HISAT_MAXINS} {BOWTIE_HISAT_PE_ORIENTATION} {BOWTIE_HISAT_NOFW_NORC} -x {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION} -1 {input.R1} -2 {input.R2} -S {output} 2> {log.err} 1> {log.out}
           """



';



}



if ($Type eq "Single-End") {



print Profile '
rule bowtie2:
    input:
         fastq=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq",
         index_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".1.bt2",
         index_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".2.bt2",
         index_3=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".3.bt2",
         index_4=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".4.bt2",
         index_5=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".rev.1.bt2",
         index_6=WORKING_FOLDER + "/standard_results/ " + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".rev.2.bt2",
         index_7=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".transcript.fa"
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bowtie2.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bowtie2.stderr"
    threads: THREADS
    shell: """
           bowtie2 -p {threads} {BOWTIE_HISAT_NOFW_NORC} -a -x {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION} -U {input.fastq} -S {output} 2> {log.err} 1> {log.out}
           """



';



}



print Profile '
rule sam2bam:
    input:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sam"
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sam2bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sam2bam.stderr"
    shell: """
           samtools view -S -b -o {output} {input} 2> {log.err} 1> {log.out}
           """



rule express:
    input:
         fa=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".transcript.fa",
         bam=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.bam"
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Quantification/{smp}/params.xprs",
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Quantification/{smp}/results.xprs"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.express.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.express.stderr"
    shell: """
           express {EXPRESS_ORIENTATION} -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Quantification/{wildcards.smp} {input.fa} {input.bam} 2> {log.err} 1> {log.out}
           """



rule express_copy:
    input:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Quantification/{smp}/results.xprs"
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Summary/{smp}.results.xprs"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.express_copy.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.express_copy.stderr"
    shell: """
           cp {input} {output} 2> {log.err} 1> {log.out}
           """



rule matrix:
    input:
          xprs=expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Summary/{sample}.results.xprs", sample = SAMPLES_2),
          transcript2gene=MAP_FILE
    output:
          transcript=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".transcript.FPKM.txt",
          gene=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/matrix.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/matrix.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/combine_eXpress_results.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Summary -j {input.transcript2gene} -o {output.transcript} 2> {log.err} 1> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/sleuth_transcript2gene.pl -i {input.transcript2gene} -j {output.transcript} -o {output.gene} 2>> {log.err} 1>> {log.out}
           """



';



###DEGs



#CLASS_testis_2 = "R_testis_7a;R_testis_7b;R_testis_7c"

#CLASS_brain_2 = "R_brain_3b;R_brain_3c;R_brain_a"



#$control[$i]
#$treatment[$i]



for (my $i=0;$i<@treatment;$i++) {



#$control[$i]
#$treatment[$i]



print Profile "\n\n\n";



print Profile "rule deg_".$control[$i]."_vs_".$treatment[$i].":\n";
print Profile "    input:\n";
print Profile "          expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Gene_expression_matrix/Summary/{sample}.results.xprs\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i].")\n";
print Profile "    output:\n";
print Profile "          expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/".$control[$i]."_vs_".$treatment[$i]."/{sample}.results.xprs\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "          WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".eff_counts.txt\",\n";
print Profile "          WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".gene.DEG.txt\",\n";
print Profile "          WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".transcript.DEG.csv\",\n";
print Profile "          WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".transcript.DEG.txt\"\n";
print Profile "    log:\n";
print Profile "          out=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/deg_".$control[$i]."_vs_".$treatment[$i].".stdout\",\n";
print Profile "          err=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/deg_".$control[$i]."_vs_".$treatment[$i].".stderr\"\n";
print Profile "    shell: \"\"\"\n";
print Profile "           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/".$control[$i]."_vs_".$treatment[$i]."\n";



for (my $j=0;$j<@{$replicate{$control[$i]}};$j++) {



print Profile "           cp {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Summary/".$replicate{$control[$i]}[$j].".results.xprs {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/".$control[$i]."_vs_".$treatment[$i]."/\n";



}



for (my $j=0;$j<@{$replicate{$treatment[$i]}};$j++) {



print Profile "           cp {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Summary/".$replicate{$treatment[$i]}[$j].".results.xprs {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/".$control[$i]."_vs_".$treatment[$i]."/\n";



}



print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/combine_eXpress_counts.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/".$control[$i]."_vs_".$treatment[$i]." -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/{PROJECT_NAME}.{WORKFLOW_NAME}.".$control[$i]."_vs_".$treatment[$i].".eff_counts.txt\n";
print Profile "           R --slave --vanilla --args ".$control[$i]." ".$treatment[$i]." \"{CLASS_".$control[$i]."_2}\" \"{CLASS_".$treatment[$i]."_2}\" {WORKING_FOLDER}/standard_results {PROJECT_NAME} < {HPPRNA_SOFTWARE_FOLDER}/my_R/{WORKFLOW_NAME}.edgeR.R\n";
print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/edgeR_DEG_gene.pl -i {MAP_FILE} -j {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/{PROJECT_NAME}.{WORKFLOW_NAME}.".$control[$i]."_vs_".$treatment[$i].".transcript.DEG.txt -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/{PROJECT_NAME}.{WORKFLOW_NAME}.".$control[$i]."_vs_".$treatment[$i].".gene.DEG.txt\n";
print Profile "           \"\"\"\n\n\n";



}#for (my $i=0;$i<@treatment;$i++) {



print Profile '
rule clean_fpkm:
    input:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.txt"
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/clean_fpkm.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/clean_fpkm.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/clean_FPKM_STAR_RSEM_EBSeq.pl -i {input} -o {output} 2> {log.err} 1> {log.out}
           """



rule heatmap:
    input:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt"
    output:
         R=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.heatmap.R",
         pdf=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.heatmap.pdf"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/heatmap.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/heatmap.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_heatmap_Rscript.pl -i {PROJECT_NAME}.{WORKFLOW_NAME}.gene.FPKM.clean.txt -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster/{PROJECT_NAME}.{WORKFLOW_NAME}.gene.FPKM.clean.heatmap.R 2> {log.err} 1> {log.out}
           R --slave --vanilla < {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster/{PROJECT_NAME}.{WORKFLOW_NAME}.gene.FPKM.clean.heatmap.R 2>> {log.err} 1>> {log.out}
           """



rule pca:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt"
    output:
        R=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.PCA.R",
        matrix=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.matrix.txt",
        pdf=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.PCA.pdf",
        rotation=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.rotation.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/pca.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/pca.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_PCA_Rscript.pl -i {PROJECT_NAME}.{WORKFLOW_NAME}.gene.FPKM.clean.txt -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



';



}



############################Workflow 4: Bowtie_eXpress_edgeR end#####################################


############################Workflow 5: Kallisto_sleuth start##########################################



if ($Core_workflow eq 5) {



print Profile '
rule kallisto_index:
    input:
         fa=HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + "/genome.fa",
         gtf=GTF
    output:
         fa=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".transcript.fa",
         idx=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".transcript.idx"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/kallisto_index.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/kallisto_index.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/NGS_plot
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Bigwig
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Snakemake_logs
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Quantification
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Summary
           gffread -w {output.fa} -g {input.fa} {input.gtf} 2> {log.err} 1> {log.out}
           kallisto index -i {output.idx} {output.fa} 2>> {log.err} 1>> {log.out}
           """



';



if ($Type eq "Paired-End") {



print Profile '
rule quant:
    input:
          R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
          R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq",
          idx=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".transcript.idx"
    output:
          sam=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sam",
          h5=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Quantification/{smp}/kallisto/abundance.h5",
          tsv_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Quantification/{smp}/kallisto/abundance.tsv",
          json=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Quantification/{smp}/kallisto/run_info.json",
          tsv_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Summary/{smp}.abundance.tsv"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.quant.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.quant.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Quantification/{wildcards.smp}
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Quantification/{wildcards.smp}/kallisto
           kallisto quant --pseudobam -i {input.idx} -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Quantification/{wildcards.smp}/kallisto -b 100 {input.R1} {input.R2} > {output.sam}
           cp {output.tsv_1} {output.tsv_2}
           """



';



}



if ($Type eq "Single-End") {



print Profile '



rule quant:
    input:
          fastq=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq",
          idx=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".transcript.idx"
    output:
          sam=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sam",
          h5=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Quantification/{smp}/kallisto/abundance.h5",
          tsv_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Quantification/{smp}/kallisto/abundance.tsv",
          json=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Quantification/{smp}/kallisto/run_info.json",
          tsv_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Summary/{smp}.abundance.tsv"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.quant.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.quant.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Quantification/{wildcards.smp}
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Quantification/{wildcards.smp}/kallisto
           kallisto quant --pseudobam -i {input.idx} -o {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Quantification/{wildcards.smp}/kallisto -b 100 --single -l {KALLISTO_L} -s {KALLISTO_S} {input.fastq} > {output.sam}
           cp {output.tsv_1} {output.tsv_2}
           """



';



}



print Profile '
rule sam2bam:
    input:
          WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sam"
    output:
          WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sam2bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sam2bam.stderr"
    shell: """
           samtools view -S -b -o {output} {input} 2> {log.err} 1> {log.out}
           """



rule matrix:
    input:
          tsv=expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Summary/{sample}.abundance.tsv", sample = SAMPLES_2),
          txt=MAP_FILE
    output:
          transcript=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".transcript.TPM.txt",
          gene=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/matrix.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/matrix.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/combine_kallisto_results.pl -i {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Summary -j {input.txt} -o {output.transcript} 2> {log.err} 1> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/sleuth_transcript2gene.pl -i {input.txt} -j {output.transcript} -o {output.gene} 2>> {log.err} 1>> {log.out}
           """



';



###DEGs



#$control[$i]
#$treatment[$i]



for (my $i=0;$i<@treatment;$i++) {



print Profile "rule deg_".$control[$i]."_vs_".$treatment[$i].":\n";
print Profile "    input:\n";
print Profile "         expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Gene_expression_matrix/Quantification/{sample}/kallisto/abundance.h5\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Gene_expression_matrix/Quantification/{sample}/kallisto/abundance.tsv\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Gene_expression_matrix/Quantification/{sample}/kallisto/run_info.json\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i].")\n";
print Profile "    output:\n";
print Profile "         R=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".sleuth.R\",\n";
print Profile "         transcript=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".transcript.DEG.txt\",\n";
print Profile "         gene=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".gene.DEG.txt\"\n";
print Profile "    log:\n";
print Profile "         out=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/deg_".$control[$i]."_vs_".$treatment[$i].".stdout\",\n";
print Profile "         err=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/deg_".$control[$i]."_vs_".$treatment[$i].".stderr\"\n";
print Profile "    shell: \"\"\"\n";
print Profile "           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/".$control[$i]."_vs_".$treatment[$i]."\n";



for (my $j=0;$j<@{$replicate{$control[$i]}};$j++) {



print Profile "           cp -r {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Quantification/".$replicate{$control[$i]}[$j]." {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/".$control[$i]."_vs_".$treatment[$i]."/\n";



}



for (my $j=0;$j<@{$replicate{$treatment[$i]}};$j++) {



print Profile "           cp -r {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Quantification/".$replicate{$treatment[$i]}[$j]." {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/".$control[$i]."_vs_".$treatment[$i]."/\n";



}



print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_comparison_info.pl -p {PROJECT_NAME} -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG -s \"{SAMPLE_join_".$control[$i]."_vs_".$treatment[$i]."}\" -g \"{SAMPLE_group_join_".$control[$i]."_vs_".$treatment[$i]."}\" -c \"".$control[$i].":".$treatment[$i].";\"\n";
print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_sleuth_Rscript.pl -p {PROJECT_NAME} -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME} -c \"".$control[$i].":".$treatment[$i].";\" -o {output.R}\n";
print Profile "           R --slave --vanilla < {output.R}\n";
print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/sleuth_DEG_gene.pl -i {MAP_FILE} -j {output.transcript} -o {output.gene}\n";
print Profile "           \"\"\"\n\n\n";



#SAMPLE_join = "R_brain_3b;R_brain_3c;R_brain_a;R_testis_7a;R_testis_7b;R_testis_7c;"

#SAMPLE_group_join = "brain;brain;brain;testis;testis;testis;"



}#for (my $i=0;$i<@treatment;$i++) {



print Profile '



rule clean_tpm:
    input:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.txt"
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.clean.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/clean_tpm.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/clean_tpm.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/clean_FPKM_STAR_RSEM_EBSeq.pl -i {input} -o {output} 2> {log.err} 1> {log.out}
           """



rule heatmap:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.clean.txt"
    output:
        R=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.clean.heatmap.R",
        pdf=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.clean.txt.1.heatmap.pdf"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/heatmap.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/heatmap.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_heatmap_Rscript.pl -i {PROJECT_NAME}.{WORKFLOW_NAME}.gene.TPM.clean.txt -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



rule pca:
    input:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.clean.txt"
    output:
         R=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.clean.PCA.R",
         matrix=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.clean.txt.1.matrix.txt",
         pdf=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.clean.txt.1.PCA.pdf",
         rotation=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.TPM.clean.txt.1.rotation.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/pca.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/pca.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_PCA_Rscript.pl -i {PROJECT_NAME}.{WORKFLOW_NAME}.gene.TPM.clean.txt -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



';



}



############################Workflow 5: Kallisto_sleuth end##########################################

############################Workflow 6: HISAT_StringTie_Ballgown start#################################



if ($Core_workflow eq 6) {



print Profile '
rule hisat2_build:
    input:
        HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + "/genome.fa"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".1.ht2",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".2.ht2",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".3.ht2",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".4.ht2",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".5.ht2",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".6.ht2",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".7.ht2",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".8.ht2"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/hisat2_build.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/hisat2_build.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Processed_BAM
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/NGS_plot
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Bigwig
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index
           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Snakemake_logs
           hisat2-build {input} {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION} 2> {log.err} 1> {log.out}
           """



';



if ($Type eq "Paired-End") {



print Profile '
rule hisat2:
    input:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq",
        index_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".1.ht2",
        index_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".2.ht2",
        index_3=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".3.ht2",
        index_4=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".4.ht2",
        index_5=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".5.ht2",
        index_6=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".6.ht2",
        index_7=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".7.ht2",
        index_8=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".8.ht2"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.hisat2.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.hisat2.stderr"
    threads: THREADS
    shell: """
           hisat2 -p {threads} -I {BOWTIE_HISAT_MININS} -X {BOWTIE_HISAT_MAXINS} {BOWTIE_HISAT_PE_ORIENTATION} -x {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION} -1 {input.R1} -2 {input.R2} -S {output} 2> {log.err} 1> {log.out}
           """



';



}#if ($Type eq "Paired-End") {



if ($Type eq "Single-End") {



print Profile '
rule hisat2:
    input:
        fastq=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq",
        index_1=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".1.ht2",
        index_2=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".2.ht2",
        index_3=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".3.ht2",
        index_4=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".4.ht2",
        index_5=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".5.ht2",
        index_6=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".6.ht2",
        index_7=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".7.ht2",
        index_8=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Index/" + ASSEMBLY_VERSION + ".8.ht2"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.hisat2.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.hisat2.stderr"
    threads: THREADS
    shell: """
           hisat2 -p {threads} -x {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Index/{ASSEMBLY_VERSION} -U {input.fastq} -S {output} 2> {log.err} 1> {log.out}
           """



';



}#if ($Type eq "Single-End") {



print Profile '
rule sam2bam:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.1.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sam2bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sam2bam.stderr"
    threads: THREADS
    shell: """
           samtools view -S -b -o {output} {input} 2> {log.err} 1> {log.out}
           """



rule sort_bam:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.1.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sort_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sort_bam.stderr"
    shell: """
           samtools sort -o {output} {input} 2> {log.err} 1> {log.out}
           """



';



if ($Type eq "Paired-End") {



print Profile '
rule rmdup:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.rmdup.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.rmdup.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.rmdup.stderr"
    shell: """
           samtools rmdup {input} {output} 2> {log.err} 1> {log.out}
           """



';



}#if ($Type eq "Paired-End") {



if ($Type eq "Single-End") {



print Profile '
rule rmdup:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.rmdup.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.rmdup.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.rmdup.stderr"
    shell: """
           samtools rmdup -s {input} {output} 2> {log.err} 1> {log.out}
           """



';



}#if ($Type eq "Single-End") {



print Profile '
rule bam2sam:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bam2sam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bam2sam.stderr"
    shell: """
           samtools view -h {input} > {output}
           """



rule unique_mapping:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.sam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.unique_mapping.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.unique_mapping.stderr"
    shell: """
           grep -E \'NH:i:1|@\' {input} > {output}
           """



rule sam2bam_2:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.sam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sam2bam_2.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.sam2bam_2.stderr"
    shell: """
           samtools view -S -b -o {output} {input} 2> {log.err} 1> {log.out}
           """



rule index_bam:
    input:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam"
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam.bai"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.index_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.index_bam.stderr"
    threads: THREADS
    shell: """
           samtools index {input} 2> {log.err} 1> {log.out}
           """



rule stringtie:
    input:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam"
    output:
         gff=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/gff/{smp}.stringtie.gff",
         e2t=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Ballgown_input/{smp}/e2t.ctab",
         e=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Ballgown_input/{smp}/e_data.ctab",
         i2t=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Ballgown_input/{smp}/i2t.ctab",
         i=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Ballgown_input/{smp}/i_data.ctab",
         t=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Ballgown_input/{smp}/t_data.ctab"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.stringtie.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.stringtie.stderr"
    threads: THREADS
    shell: """
           stringtie -G {GTF} -b {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Ballgown_input/{wildcards.smp} -p {threads} {input} -o {output.gff} 2> {log.err} 1> {log.out}
           """



rule matrix:
    input:
         e2t=expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Ballgown_input/{sample}/e2t.ctab", sample = SAMPLES_2),
         e=expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Ballgown_input/{sample}/e_data.ctab", sample = SAMPLES_2),
         i2t=expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Ballgown_input/{sample}/i2t.ctab", sample = SAMPLES_2),
         i=expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Ballgown_input/{sample}/i_data.ctab", sample = SAMPLES_2),
         t=expand(WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/Ballgown_input/{sample}/t_data.ctab", sample = SAMPLES_2)
    output:
         R=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + ".Ballgown.matrix.R",
         gene=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.txt",
         transcript=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".transcript.FPKM.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/matrix.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/matrix.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_Ballgown_matrix_Rscript.pl -p {PROJECT_NAME} -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME} -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



';



###DEGs



#$control[$i]
#$treatment[$i]



for (my $i=0;$i<@treatment;$i++) {



print Profile "rule deg_".$control[$i]."_vs_".$treatment[$i].":\n";
print Profile "    input:\n";
print Profile "         e2t=expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Gene_expression_matrix/Ballgown_input/{sample}/e2t.ctab\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         e=expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Gene_expression_matrix/Ballgown_input/{sample}/e_data.ctab\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         i2t=expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Gene_expression_matrix/Ballgown_input/{sample}/i2t.ctab\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         i=expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Gene_expression_matrix/Ballgown_input/{sample}/i_data.ctab\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         t=expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Gene_expression_matrix/Ballgown_input/{sample}/t_data.ctab\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i].")\n";
print Profile "    output:\n";
print Profile "         e2t=expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/".$control[$i]."_vs_".$treatment[$i]."/{sample}/e2t.ctab\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         e=expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/".$control[$i]."_vs_".$treatment[$i]."/{sample}/e_data.ctab\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         i2t=expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/".$control[$i]."_vs_".$treatment[$i]."/{sample}/i2t.ctab\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         i=expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/".$control[$i]."_vs_".$treatment[$i]."/{sample}/i_data.ctab\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         t=expand(WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/".$control[$i]."_vs_".$treatment[$i]."/{sample}/t_data.ctab\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         R=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".Ballgown.DEG.R\",\n";
print Profile "         gene=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".$control[$i]_vs_$treatment[$i].gene.DEG.txt\",\n";
print Profile "         transcript=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/DEG/\" + PROJECT_NAME + \".\" + WORKFLOW_NAME + \".".$control[$i]."_vs_".$treatment[$i].".transcript.DEG.txt\"\n";
print Profile "    log:\n";
print Profile "         out=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/deg_".$control[$i]."_vs_".$treatment[$i].".stdout\",\n";
print Profile "         err=WORKING_FOLDER + \"/standard_results/\" + WORKFLOW_NAME + \"/Snakemake_logs/deg_".$control[$i]."_vs_".$treatment[$i].".stderr\"\n";
print Profile "    shell: \"\"\"\n";
print Profile "           mkdir -p {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/".$control[$i]."_vs_".$treatment[$i]."\n";



for (my $j=0;$j<@{$replicate{$control[$i]}};$j++) {



print Profile "           cp -r {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Ballgown_input/".$replicate{$control[$i]}[$j]." {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/".$control[$i]."_vs_".$treatment[$i]."/\n";



}



for (my $j=0;$j<@{$replicate{$treatment[$i]}};$j++) {



print Profile "           cp -r {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Gene_expression_matrix/Ballgown_input/".$replicate{$treatment[$i]}[$j]." {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/DEG/".$control[$i]."_vs_".$treatment[$i]."/\n";



}



print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_Ballgown_DEG_Rscript.pl -p {PROJECT_NAME} -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME} -s \"{SAMPLE_join_".$control[$i]."_vs_".$treatment[$i]."}\" -g \"{SAMPLE_group_join_".$control[$i]."_vs_".$treatment[$i]."}\" -c \"".$control[$i].":".$treatment[$i].";\" -o {output.R}\n";
print Profile "           R --slave --vanilla < {output.R}\n";
print Profile "           \"\"\"\n";



#SAMPLE_join = "R_brain_3b;R_brain_3c;R_brain_a;R_testis_7a;R_testis_7b;R_testis_7c;"

#SAMPLE_group_join = "brain;brain;brain;testis;testis;testis;"



}#for (my $i=0;$i<@treatment;$i++) {



print Profile '
rule clean_fpkm:
    input:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Gene_expression_matrix/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.txt"
    output:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/clean_fpkm.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/clean_fpkm.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/Ballgown_fix_FPKM_title.pl -i {input} -o {output} 2> {log.err} 1> {log.out}
           """



rule ngs_plot:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{smp}.genebody.avgprof.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{smp}.genebody.heatmap.pdf",
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/NGS_plot/NGSplot.{smp}.genebody.zip"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.ngs_plot.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.ngs_plot.stderr"
    threads: THREADS
    shell: """
           ngs.plot.r -G {ASSEMBLY_VERSION} -R genebody -C {input} -O {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/NGS_plot/NGSplot.{wildcards.smp}.genebody -T {wildcards.smp} -L 3000 -RB 0.05 -F rnaseq 2> {log.err} 1> {log.out}
           """



rule bam2wig:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Processed_BAM/{smp}.sorted.unique.bam"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.gz"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bam2wig.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.bam2wig.stderr"
    shell: """
           samtools mpileup -BQ0 {input} | perl -ne \'BEGIN{{print "track type=wiggle_0 name={wildcards.smp} maxHeightPixels=64:64:11 color=31,120,180 visibility=full\\\\n"}};($c, $start, undef, $depth) = split; if ($c ne $lastC) {{ print "variableStep chrom=$c\\\\n"; }};$lastC=$c;next unless $. % 5 ==0;print "$start\\\\t$depth\\\\n" unless $depth<3;\' | gzip -c > {output}
           """



rule wig2wig2:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.gz"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.2.gz"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig2wig2.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig2wig2.stderr"
    shell: """
           zcat {input} | sed \'1d\' | gzip -c > {output}
           """



rule wig22bw:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.wig.2.gz"
    output:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Bigwig/{smp}.bw"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig22bw.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/{smp}.wig22bw.stderr"
    shell: """
           {WIGTOBIGWIG} {input} {CHROM_SIZES} {output} 2> {log.err} 1> {log.out}
           """



rule heatmap:
    input:
         WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt"
    output:
         R=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.heatmap.R",
         pdf=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.heatmap.pdf"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/heatmap.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/heatmap.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_heatmap_Rscript.pl -i {PROJECT_NAME}.{WORKFLOW_NAME}.gene.FPKM.clean.txt -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



rule pca:
    input:
        WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt"
    output:
        R=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.PCA.R",
        matrix=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.matrix.txt",
        pdf=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.PCA.pdf",
        rotation=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Cluster/" + PROJECT_NAME + "." + WORKFLOW_NAME + ".gene.FPKM.clean.txt.1.rotation.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/pca.stdout",
        err=WORKING_FOLDER + "/standard_results/" + WORKFLOW_NAME + "/Snakemake_logs/pca.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_PCA_Rscript.pl -i {PROJECT_NAME}.{WORKFLOW_NAME}.gene.FPKM.clean.txt -d {WORKING_FOLDER}/standard_results/{WORKFLOW_NAME}/Cluster -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



';



}#if ($Core_workflow eq 6) {



############################Workflow 6: HISAT_StringTie_Ballgown end#################################


###############################################Core workflow end############################################################



}#if (($Analysis_type eq "protein-coding")||($Analysis_type eq "known lncRNA")) {




###############################################lncRNA de novo workflow start############################################################



if ($Analysis_type eq "novel lncRNA") {



##################################Paired-End start###################################



if ($Type eq "Paired-End") {



print Profile '
rule bowtie2_build:
    input:
        GENOME_FA_FILE
    output:
        fa=TOPHAT_LNCRNA_INDEX_FILE + ".fa",
        index_1=TOPHAT_LNCRNA_INDEX_FILE + ".1.bt2",
        index_2=TOPHAT_LNCRNA_INDEX_FILE + ".2.bt2",
        index_3=TOPHAT_LNCRNA_INDEX_FILE + ".3.bt2",
        index_4=TOPHAT_LNCRNA_INDEX_FILE + ".4.bt2",
        index_5=TOPHAT_LNCRNA_INDEX_FILE + ".rev.1.bt2",
        index_6=TOPHAT_LNCRNA_INDEX_FILE + ".rev.2.bt2"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/bowtie2_build.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/bowtie2_build.stderr"
    threads: THREADS
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Index
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/Snakemake_logs
           cp {input} {output.fa} 2> {log.err} 1> {log.out}
           bowtie2-build --threads {threads} {output.fa} {TOPHAT_LNCRNA_INDEX_FILE} 2>> {log.err} 1>> {log.out}
           """



rule tophat2:
    input:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq",
        index_1=TOPHAT_LNCRNA_INDEX_FILE + ".1.bt2",
        index_2=TOPHAT_LNCRNA_INDEX_FILE + ".2.bt2",
        index_3=TOPHAT_LNCRNA_INDEX_FILE + ".3.bt2",
        index_4=TOPHAT_LNCRNA_INDEX_FILE + ".4.bt2",
        index_5=TOPHAT_LNCRNA_INDEX_FILE + ".rev.1.bt2",
        index_6=TOPHAT_LNCRNA_INDEX_FILE + ".rev.2.bt2"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}_thout/accepted_hits.bam",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}_thout/align_summary.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.tophat2.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.tophat2.stderr"
    threads: THREADS
    shell: """
           tophat2 -p {threads} --mate-inner-dist {TOPHAT_MATE_INNER_DIST} --mate-std-dev {TOPHAT_MATE_STD_DEV} --library-type {TOPHAT_LIBRARY_TYPE} -G {GTF} -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{wildcards.smp}_thout {TOPHAT_LNCRNA_INDEX_FILE} {input.R1} {input.R2} 2> {log.err} 1> {log.out}
           """



rule copy_summary:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}_thout/align_summary.txt"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}_align_summary.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.copy_summary.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.copy_summary.stderr"
    shell: """
           cp {input} {output} 2> {log.err} 1> {log.out}
           """



rule sort_bam:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}_thout/accepted_hits.bam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.sort_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.sort_bam.stderr"
    shell: """
           samtools sort -o {output} {input} 2> {log.err} 1> {log.out}
           """



rule rmdup_bam:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.bam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.rmdup.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.rmdup_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.rmdup_bam.stderr"
    shell: """
           samtools rmdup {input} {output} 2> {log.err} 1> {log.out}
           """



rule bam2sam:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.bam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.bam2sam.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.bam2sam.stderr"
    shell: """
           samtools view -h {input} > {output}
           """



rule select_unique_mapping:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.sam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.unique.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.select_unique_mapping.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.select_unique_mapping.stderr"
    shell: """
           perl {EXTRACT_UNIQUE} --input_sam {input} --output {output} --aligner T 2> {log.err} 1> {log.out}
           """



rule sam2bam:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.unique.sam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.unique.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.sam2bam.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.sam2bam.stderr"
    shell: """
           samtools view -S -b -o {output} {input} 2> {log.err} 1> {log.out}
           """



rule index_bam:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.unique.bam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.unique.bam.bai"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.index_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.index_bam.stderr"
    shell: """
           samtools index {input} 2> {log.err} 1> {log.out}
           """



rule cufflinks:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.unique.bam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/{smp}_clout/transcripts.gtf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.cufflinks.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.cufflinks.stderr"
    threads: THREADS
    shell: """
           cufflinks -p {threads} --library-type {CUFFLINKS_LIBRARY_TYPE} -g {GTF} -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/{wildcards.smp}_clout {input} 2> {log.err} 1> {log.out}
           """



rule cuffmerge:
    input:
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/{sample}_clout/transcripts.gtf", sample = SAMPLES_2),
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/assembly_list.txt",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/merged.gtf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/cuffmerge.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/cuffmerge.stderr"
    threads: THREADS
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_cuffmerge_assembly_list.pl -d {WORKING_FOLDER} -s "{SAMPLE_join}" -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/assembly_list.txt 2> {log.err} 1> {log.out}
           cuffmerge -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly -g {GTF} -s {GENOME_FA_FILE} -p {threads} {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/assembly_list.txt 2>> {log.err} 1>> {log.out}
           """



rule cuffcompare:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/merged.gtf"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp.combined.gtf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/cuffcompare.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/cuffcompare.stderr"
    threads: THREADS
    shell: """
           cuffcompare -r {GTF} -s {GENOME_FA_FILE} -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/merged.gtf 2> {log.err} 1> {log.out}
           """



rule select_cuffcompare_class_code:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp.combined.gtf"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp.combined.selected.gtf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/select_cuffcompare_class_code.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/select_cuffcompare_class_code.stderr"
    threads: THREADS
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/select_cuffcompare_class_code.pl -i {input} -o {output} 2> {log.err} 1> {log.out}
           """



rule iSeeRNA:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp.combined.selected.gtf"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA/cuffcmp.combined.selected.result"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/iSeeRNA.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/iSeeRNA.stderr"
    threads: THREADS
    shell: """
           iSeeRNA -c {HPPRNA_SOFTWARE_FOLDER}/iSeeRNA-1.2.2/conf/{ASSEMBLY_VERSION}.conf -i {input} -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA 2> {log.err} 1> {log.out}
           cd {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA 2>> {log.err} 1>> {log.out}
           make 2>> {log.err} 1>> {log.out}
           """



rule identify_noncoding:
    input:
        result=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA/cuffcmp.combined.selected.result",
        gtf=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp.combined.selected.gtf"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA/cuffcmp.combined.selected.noncoding.gtf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/identify_noncoding.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/identify_noncoding.stderr"
    threads: THREADS
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/identify_noncoding.pl -i {input.result} -j {input.gtf} -o {output} 2> {log.err} 1> {log.out}
           """



rule compile_gtf:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA/cuffcmp.combined.selected.noncoding.gtf"
    output:
        denovo_map=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA/denovo_lncRNA.transcript_to_gene.txt",
        all_map=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/all_lncRNA_mRNA.transcript_to_gene.txt",
        all_gtf=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/all_lncRNA_mRNA.gtf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/compile_gtf.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/compile_gtf.stderr"
    threads: THREADS
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/get_transcript_to_gene_table.pl -i {input} -o {output.denovo_map} 2> {log.err} 1> {log.out}
           cat {output.denovo_map} {MAP_FILE} > {output.all_map}
           cat {input} {GTF} > {output.all_gtf}
           """



rule kallisto_index:
    input:
         fa=HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + "/genome.fa",
         gtf=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/all_lncRNA_mRNA.gtf"
    output:
         fa=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/" + ASSEMBLY_VERSION + ".transcript.fa",
         idx=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/" + ASSEMBLY_VERSION + ".transcript.idx"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/kallisto_index.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/kallisto_index.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Processed_BAM
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/DEG
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/NGS_plot
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Bigwig
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Index
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Summary
           gffread -w {output.fa} -g {input.fa} {input.gtf} 2> {log.err} 1> {log.out}
           kallisto index -i {output.idx} {output.fa} 2>> {log.err} 1>> {log.out}
           """



rule quant:
    input:
          R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
          R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq",
          idx=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/" + ASSEMBLY_VERSION + ".transcript.idx"
    output:
          sam=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Processed_BAM/{smp}.sam",
          h5=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{smp}/kallisto/abundance.h5",
          tsv_1=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{smp}/kallisto/abundance.tsv",
          json=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{smp}/kallisto/run_info.json",
          tsv_2=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Summary/{smp}.abundance.tsv"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.quant.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.quant.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{wildcards.smp}
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{wildcards.smp}/kallisto
           kallisto quant --pseudobam -i {input.idx} -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{wildcards.smp}/kallisto -b 100 {input.R1} {input.R2} > {output.sam}
           cp {output.tsv_1} {output.tsv_2}
           """



rule sam2bam_2:
    input:
          WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Processed_BAM/{smp}.sam"
    output:
          WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Processed_BAM/{smp}.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.sam2bam_2.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.sam2bam_2.stderr"
    shell: """
           samtools view -S -b -o {output} {input} 2> {log.err} 1> {log.out}
           """



rule matrix:
    input:
          tsv=expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Summary/{sample}.abundance.tsv", sample = SAMPLES_2),
          txt=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/all_lncRNA_mRNA.transcript_to_gene.txt"
    output:
          transcript=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/" + PROJECT_NAME + ".lncRNA_denovo.transcript.TPM.txt",
          gene=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/matrix.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/matrix.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/combine_kallisto_results.pl -i {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Summary -j {input.txt} -o {output.transcript} 2> {log.err} 1> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/sleuth_transcript2gene.pl -i {input.txt} -j {output.transcript} -o {output.gene} 2>> {log.err} 1>> {log.out}
           """



';



for (my $i=0;$i<@treatment;$i++) {



print Profile "rule deg_".$control[$i]."_vs_".$treatment[$i].":\n";
print Profile "    input:\n";
print Profile "         map=WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/Index/all_lncRNA_mRNA.transcript_to_gene.txt\",\n";
print Profile "         h5=expand(WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{sample}/kallisto/abundance.h5\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         tsv=expand(WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{sample}/kallisto/abundance.tsv\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         json=expand(WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{sample}/kallisto/run_info.json\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i].")\n";
print Profile "    output:\n";
print Profile "         R=WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/\" + PROJECT_NAME + \".lncRNA_denovo.".$control[$i]."_vs_".$treatment[$i].".sleuth.R\",\n";
print Profile "         transcript=WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/\" + PROJECT_NAME + \".lncRNA_denovo.".$control[$i]."_vs_".$treatment[$i].".transcript.DEG.txt\",\n";
print Profile "         gene=WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/\" + PROJECT_NAME + \".lncRNA_denovo.".$control[$i]."_vs_".$treatment[$i].".gene.DEG.txt\"\n";
print Profile "    log:\n";
print Profile "         out=WORKING_FOLDER + \"/standard_results/lncRNA_denovo/Snakemake_logs/deg_".$control[$i]."_vs_".$treatment[$i].".stdout\",\n";
print Profile "         err=WORKING_FOLDER + \"/standard_results/lncRNA_denovo/Snakemake_logs/deg_".$control[$i]."_vs_".$treatment[$i].".stderr\"\n";
print Profile "    shell: \"\"\"\n";
print Profile "           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/".$control[$i]."_vs_".$treatment[$i]."\n";



for (my $j=0;$j<@{$replicate{$control[$i]}};$j++) {



print Profile "           cp -r {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/".$replicate{$control[$i]}[$j]." {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/".$control[$i]."_vs_".$treatment[$i]."/\n";



}



for (my $j=0;$j<@{$replicate{$treatment[$i]}};$j++) {



print Profile "           cp -r {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/".$replicate{$treatment[$i]}[$j]." {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/".$control[$i]."_vs_".$treatment[$i]."/\n";



}



print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_comparison_info_lncRNA.pl -p {PROJECT_NAME} -d {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/DEG -s \"{SAMPLE_join_".$control[$i]."_vs_".$treatment[$i]."}\" -g \"{SAMPLE_group_join_".$control[$i]."_vs_".$treatment[$i]."}\" -c \"".$control[$i].":".$treatment[$i].";\"\n";
print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_sleuth_Rscript_lncRNA.pl -p {PROJECT_NAME} -d {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification -c \"".$control[$i].":".$treatment[$i].";\" -o {output.R}\n";
print Profile "           R --slave --vanilla < {output.R}\n";
print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/sleuth_DEG_gene.pl -i {input.map} -j {output.transcript} -o {output.gene}\n";
print Profile "           \"\"\"\n\n\n";



}#for (my $i=0;$i<@treatment;$i++) {



print Profile '
rule clean_tpm:
    input:
         WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.txt"
    output:
         WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/clean_tpm.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/clean_tpm.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/clean_FPKM_STAR_RSEM_EBSeq.pl -i {input} -o {output} 2> {log.err} 1> {log.out}
           """



rule heatmap:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt"
    output:
        R=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.heatmap.R",
        pdf=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt.1.heatmap.pdf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/heatmap.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/heatmap.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_heatmap_Rscript.pl -i {PROJECT_NAME}.lncRNA_denovo.gene.TPM.clean.txt -d {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



rule pca:
    input:
         WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt"
    output:
         R=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.PCA.R",
         matrix=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt.1.matrix.txt",
         pdf=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt.1.PCA.pdf",
         rotation=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt.1.rotation.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/pca.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/pca.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_PCA_Rscript.pl -i {PROJECT_NAME}.lncRNA_denovo.gene.TPM.clean.txt -d {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



';



}#if ($Type eq "Paired-End") {



##################################Paired-End end###################################




##################################Single-End start###################################



if ($Type eq "Single-End") {



print Profile '
rule bowtie2_build:
    input:
        GENOME_FA_FILE
    output:
        fa=TOPHAT_LNCRNA_INDEX_FILE + ".fa",
        index_1=TOPHAT_LNCRNA_INDEX_FILE + ".1.bt2",
        index_2=TOPHAT_LNCRNA_INDEX_FILE + ".2.bt2",
        index_3=TOPHAT_LNCRNA_INDEX_FILE + ".3.bt2",
        index_4=TOPHAT_LNCRNA_INDEX_FILE + ".4.bt2",
        index_5=TOPHAT_LNCRNA_INDEX_FILE + ".rev.1.bt2",
        index_6=TOPHAT_LNCRNA_INDEX_FILE + ".rev.2.bt2"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/bowtie2_build.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/bowtie2_build.stderr"
    threads: THREADS
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Index
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/Snakemake_logs
           cp {input} {output.fa} 2> {log.err} 1> {log.out}
           bowtie2-build --threads {threads} {output.fa} {TOPHAT_LNCRNA_INDEX_FILE} 2>> {log.err} 1>> {log.out}
           """



rule tophat2:
    input:
        fastq=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq",
        index_1=TOPHAT_LNCRNA_INDEX_FILE + ".1.bt2",
        index_2=TOPHAT_LNCRNA_INDEX_FILE + ".2.bt2",
        index_3=TOPHAT_LNCRNA_INDEX_FILE + ".3.bt2",
        index_4=TOPHAT_LNCRNA_INDEX_FILE + ".4.bt2",
        index_5=TOPHAT_LNCRNA_INDEX_FILE + ".rev.1.bt2",
        index_6=TOPHAT_LNCRNA_INDEX_FILE + ".rev.2.bt2"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}_thout/accepted_hits.bam",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}_thout/align_summary.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.tophat2.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.tophat2.stderr"
    threads: THREADS
    shell: """
           tophat2 -p {threads} --library-type {TOPHAT_LIBRARY_TYPE} -G {GTF} -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{wildcards.smp}_thout {TOPHAT_LNCRNA_INDEX_FILE} {input.fastq} 2> {log.err} 1> {log.out}
           """



rule copy_summary:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}_thout/align_summary.txt"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}_align_summary.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.copy_summary.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.copy_summary.stderr"
    shell: """
           cp {input} {output} 2> {log.err} 1> {log.out}
           """



rule sort_bam:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}_thout/accepted_hits.bam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.sort_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.sort_bam.stderr"
    shell: """
           samtools sort -o {output} {input} 2> {log.err} 1> {log.out}
           """




rule rmdup_bam:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.bam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.rmdup.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.rmdup_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.rmdup_bam.stderr"
    shell: """
           samtools rmdup -s {input} {output} 2> {log.err} 1> {log.out}
           """



rule bam2sam:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.bam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.bam2sam.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.bam2sam.stderr"
    shell: """
           samtools view -h {input} > {output}
           """



rule select_unique_mapping:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.sam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.unique.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.select_unique_mapping.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.select_unique_mapping.stderr"
    shell: """
           perl {EXTRACT_UNIQUE} --input_sam {input} --output {output} --aligner T 2> {log.err} 1> {log.out}
           """



rule sam2bam:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.unique.sam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.unique.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.sam2bam.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.sam2bam.stderr"
    shell: """
           samtools view -S -b -o {output} {input} 2> {log.err} 1> {log.out}
           """



rule index_bam:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.unique.bam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.unique.bam.bai"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.index_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.index_bam.stderr"
    shell: """
           samtools index {input} 2> {log.err} 1> {log.out}
           """



rule cufflinks:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Processed_BAM/{smp}.sorted.unique.bam"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/{smp}_clout/transcripts.gtf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.cufflinks.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.cufflinks.stderr"
    threads: THREADS
    shell: """
           cufflinks -p {threads} --library-type {CUFFLINKS_LIBRARY_TYPE} -g {GTF} -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/{wildcards.smp}_clout {input} 2> {log.err} 1> {log.out}
           """



rule cuffmerge:
    input:
        expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/{sample}_clout/transcripts.gtf", sample = SAMPLES_2)
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/assembly_list.txt",
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/merged.gtf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/cuffmerge.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/cuffmerge.stderr"
    threads: THREADS
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_cuffmerge_assembly_list.pl -d {WORKING_FOLDER} -s "{SAMPLE_join}" -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/assembly_list.txt 2> {log.err} 1> {log.out}
           cuffmerge -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly -g {GTF} -s {GENOME_FA_FILE} -p {threads} {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/assembly_list.txt 2>> {log.err} 1>> {log.out}
           """



rule cuffcompare:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/merged.gtf"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp.combined.gtf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/cuffcompare.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/cuffcompare.stderr"
    threads: THREADS
    shell: """
           cuffcompare -r {GTF} -s {GENOME_FA_FILE} -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/merged.gtf 2> {log.err} 1> {log.out}
           """



rule select_cuffcompare_class_code:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp.combined.gtf"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp.combined.selected.gtf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/select_cuffcompare_class_code.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/select_cuffcompare_class_code.stderr"
    threads: THREADS
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/select_cuffcompare_class_code.pl -i {input} -o {output} 2> {log.err} 1> {log.out}
           """



rule iSeeRNA:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp.combined.selected.gtf"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA/cuffcmp.combined.selected.result"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/iSeeRNA.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/iSeeRNA.stderr"
    threads: THREADS
    shell: """
           iSeeRNA -c {HPPRNA_SOFTWARE_FOLDER}/iSeeRNA-1.2.2/conf/{ASSEMBLY_VERSION}.conf -i {input} -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA 2> {log.err} 1> {log.out}
           cd {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA 2>> {log.err} 1>> {log.out}
           make 2>> {log.err} 1>> {log.out}
           """



rule identify_noncoding:
    input:
        result=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA/cuffcmp.combined.selected.result",
        gtf=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/Assembly/cuffcmp.combined.selected.gtf"
    output:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA/cuffcmp.combined.selected.noncoding.gtf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/identify_noncoding.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/identify_noncoding.stderr"
    threads: THREADS
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/identify_noncoding.pl -i {input.result} -j {input.gtf} -o {output} 2> {log.err} 1> {log.out}
           """



rule compile_gtf:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA/cuffcmp.combined.selected.noncoding.gtf"
    output:
        denovo_map=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_detection/iseeRNA/denovo_lncRNA.transcript_to_gene.txt",
        all_map=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/all_lncRNA_mRNA.transcript_to_gene.txt",
        all_gtf=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/all_lncRNA_mRNA.gtf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/compile_gtf.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/compile_gtf.stderr"
    threads: THREADS
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/get_transcript_to_gene_table.pl -i {input} -o {output.denovo_map} 2> {log.err} 1> {log.out}
           cat {output.denovo_map} {MAP_FILE} > {output.all_map}
           cat {input} {GTF} > {output.all_gtf}
           """



rule kallisto_index:
    input:
         fa=HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + "/genome.fa",
         gtf=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/all_lncRNA_mRNA.gtf"
    output:
         fa=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/" + ASSEMBLY_VERSION + ".transcript.fa",
         idx=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/" + ASSEMBLY_VERSION + ".transcript.idx"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/kallisto_index.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/kallisto_index.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Processed_BAM
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/DEG
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/NGS_plot
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Bigwig
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Index
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Summary
           gffread -w {output.fa} -g {input.fa} {input.gtf} 2> {log.err} 1> {log.out}
           kallisto index -i {output.idx} {output.fa} 2>> {log.err} 1>> {log.out}
           """



rule quant:
    input:
          fastq=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq",
          idx=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/" + ASSEMBLY_VERSION + ".transcript.idx"
    output:
          sam=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Processed_BAM/{smp}.sam",
          h5=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{smp}/kallisto/abundance.h5",
          tsv_1=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{smp}/kallisto/abundance.tsv",
          json=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{smp}/kallisto/run_info.json",
          tsv_2=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Summary/{smp}.abundance.tsv"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.quant.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.quant.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{wildcards.smp}
           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{wildcards.smp}/kallisto
           kallisto quant --pseudobam -i {input.idx} -o {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{wildcards.smp}/kallisto -b 100 --single -l {KALLISTO_L} -s {KALLISTO_S} {input.fastq} > {output.sam}
           cp {output.tsv_1} {output.tsv_2}
           """



rule sam2bam_2:
    input:
          WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Processed_BAM/{smp}.sam"
    output:
          WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Processed_BAM/{smp}.bam"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.sam2bam_2.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/{smp}.sam2bam_2.stderr"
    shell: """
           samtools view -S -b -o {output} {input} 2> {log.err} 1> {log.out}
           """



rule matrix:
    input:
          tsv=expand(WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Summary/{sample}.abundance.tsv", sample = SAMPLES_2),
          txt=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Index/all_lncRNA_mRNA.transcript_to_gene.txt"
    output:
          transcript=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/" + PROJECT_NAME + ".lncRNA_denovo.transcript.TPM.txt",
          gene=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/matrix.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/matrix.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/combine_kallisto_results.pl -i {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Summary -j {input.txt} -o {output.transcript} 2> {log.err} 1> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/sleuth_transcript2gene.pl -i {input.txt} -j {output.transcript} -o {output.gene} 2>> {log.err} 1>> {log.out}
           """



';



for (my $i=0;$i<@treatment;$i++) {



print Profile "rule deg_".$control[$i]."_vs_".$treatment[$i].":\n";
print Profile "    input:\n";
print Profile "         map=WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/Index/all_lncRNA_mRNA.transcript_to_gene.txt\",\n";
print Profile "         h5=expand(WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{sample}/kallisto/abundance.h5\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         tsv=expand(WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{sample}/kallisto/abundance.tsv\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "         json=expand(WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/{sample}/kallisto/run_info.json\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i].")\n";
print Profile "    output:\n";
print Profile "         R=WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/\" + PROJECT_NAME + \".lncRNA_denovo.".$control[$i]."_vs_".$treatment[$i].".sleuth.R\",\n";
print Profile "         transcript=WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/\" + PROJECT_NAME + \".lncRNA_denovo.".$control[$i]."_vs_".$treatment[$i].".transcript.DEG.txt\",\n";
print Profile "         gene=WORKING_FOLDER + \"/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/\" + PROJECT_NAME + \".lncRNA_denovo.".$control[$i]."_vs_".$treatment[$i].".gene.DEG.txt\"\n";
print Profile "    log:\n";
print Profile "         out=WORKING_FOLDER + \"/standard_results/lncRNA_denovo/Snakemake_logs/deg_".$control[$i]."_vs_".$treatment[$i].".stdout\",\n";
print Profile "         err=WORKING_FOLDER + \"/standard_results/lncRNA_denovo/Snakemake_logs/deg_".$control[$i]."_vs_".$treatment[$i].".stderr\"\n";
print Profile "    shell: \"\"\"\n";
print Profile "           mkdir -p {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/".$control[$i]."_vs_".$treatment[$i]."\n";



for (my $j=0;$j<@{$replicate{$control[$i]}};$j++) {



print Profile "           cp -r {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/".$replicate{$control[$i]}[$j]." {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/".$control[$i]."_vs_".$treatment[$i]."/\n";



}



for (my $j=0;$j<@{$replicate{$treatment[$i]}};$j++) {



print Profile "           cp -r {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/Quantification/".$replicate{$treatment[$i]}[$j]." {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/DEG/".$control[$i]."_vs_".$treatment[$i]."/\n";



}



print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_comparison_info_lncRNA.pl -p {PROJECT_NAME} -d {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/DEG -s \"{SAMPLE_join_".$control[$i]."_vs_".$treatment[$i]."}\" -g \"{SAMPLE_group_join_".$control[$i]."_vs_".$treatment[$i]."}\" -c \"".$control[$i].":".$treatment[$i].";\"\n";
print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_sleuth_Rscript_lncRNA.pl -p {PROJECT_NAME} -d {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification -c \"".$control[$i].":".$treatment[$i].";\" -o {output.R}\n";
print Profile "           R --slave --vanilla < {output.R}\n";
print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/sleuth_DEG_gene.pl -i {input.map} -j {output.transcript} -o {output.gene}\n";
print Profile "           \"\"\"\n\n\n";



}#for (my $i=0;$i<@treatment;$i++) {



print Profile '
rule clean_tpm:
    input:
         WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Gene_expression_matrix/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.txt"
    output:
         WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/clean_tpm.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/clean_tpm.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/clean_FPKM_STAR_RSEM_EBSeq.pl -i {input} -o {output} 2> {log.err} 1> {log.out}
           """



rule heatmap:
    input:
        WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt"
    output:
        R=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.heatmap.R",
        pdf=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt.1.heatmap.pdf"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/heatmap.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/heatmap.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_heatmap_Rscript.pl -i {PROJECT_NAME}.lncRNA_denovo.gene.TPM.clean.txt -d {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



rule pca:
    input:
         WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt"
    output:
         R=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.PCA.R",
         matrix=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt.1.matrix.txt",
         pdf=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt.1.PCA.pdf",
         rotation=WORKING_FOLDER + "/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster/" + PROJECT_NAME + ".lncRNA_denovo.gene.TPM.clean.txt.1.rotation.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/pca.stdout",
        err=WORKING_FOLDER + "/standard_results/lncRNA_denovo/Snakemake_logs/pca.stderr"
    shell: """
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_PCA_Rscript.pl -i {PROJECT_NAME}.lncRNA_denovo.gene.TPM.clean.txt -d {WORKING_FOLDER}/standard_results/lncRNA_denovo/lncRNA_quantification/Cluster -o {output.R} 2> {log.err} 1> {log.out}
           R --slave --vanilla < {output.R} 2>> {log.err} 1>> {log.out}
           """



';



}#if ($Type eq "Single-End") {



##################################Single-End end###################################









###############################################lncRNA de novo workflow end############################################################



}#if ($Analysis_type eq "novel lncRNA") {



###############################################circRNA workflow start############################################################



if ($Analysis_type eq "circRNA") {



##################################Paired-End start###################################



if ($Type eq "Paired-End") {



print Profile '
rule STAR_index:
    input:
        GENOME_FA_FILE
    output:
        WORKING_FOLDER + "/standard_results/CircRNA/Index/SA",
        WORKING_FOLDER + "/standard_results/CircRNA/Index/SAindex"
    log:
        out=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/STAR_index.stdout",
        err=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/STAR_index.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/Mapping
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/CircTest
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/Index
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/Snakemake_logs
           STAR --runMode genomeGenerate --genomeDir {WORKING_FOLDER}/standard_results/CircRNA/Index --genomeFastaFiles {input} --runThreadN 1 2> {log.err} 1> {log.out}
           """



rule star_joint_mapping:
    input:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq",
        index_1=WORKING_FOLDER + "/standard_results/CircRNA/Index/SA",
        index_2=WORKING_FOLDER + "/standard_results/CircRNA/Index/SAindex"
    output:
        WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{smp}/joint_mapping/{smp}.Aligned.sortedByCoord.out.bam",
        WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{smp}/joint_mapping/{smp}.Chimeric.out.junction",
        WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{smp}/mate1_mapping/{smp}.Chimeric.out.junction",
        WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{smp}/mate2_mapping/{smp}.Chimeric.out.junction"
    log:
        out=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/{smp}.star_joint_mapping.stdout",
        err=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/{smp}.star_joint_mapping.stderr"
    threads: THREADS
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/Mapping/{wildcards.smp}
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/Mapping/{wildcards.smp}/joint_mapping
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/Mapping/{wildcards.smp}/mate1_mapping
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/Mapping/{wildcards.smp}/mate2_mapping
           STAR --runThreadN {threads} --genomeDir {WORKING_FOLDER}/standard_results/CircRNA/Index --outSAMtype BAM SortedByCoordinate --readFilesIn {input.R1} {input.R2} --outFileNamePrefix {WORKING_FOLDER}/standard_results/CircRNA/Mapping/{wildcards.smp}/joint_mapping/{wildcards.smp}. --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 2> {log.err} 1> {log.out}
           STAR --runThreadN {threads} --genomeDir {WORKING_FOLDER}/standard_results/CircRNA/Index --outSAMtype None --readFilesIn {input.R1} --outFileNamePrefix {WORKING_FOLDER}/standard_results/CircRNA/Mapping/{wildcards.smp}/mate1_mapping/{wildcards.smp}. --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 2>> {log.err} 1>> {log.out}
           STAR --runThreadN {threads} --genomeDir {WORKING_FOLDER}/standard_results/CircRNA/Index --outSAMtype None --readFilesIn {input.R2} --outFileNamePrefix {WORKING_FOLDER}/standard_results/CircRNA/Mapping/{wildcards.smp}/mate2_mapping/{wildcards.smp}. --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 2>> {log.err} 1>> {log.out}
           """



rule index_bam:
    input:
        WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{smp}/joint_mapping/{smp}.Aligned.sortedByCoord.out.bam"
    output:
        WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{smp}/joint_mapping/{smp}.Aligned.sortedByCoord.out.bam.bai"
    log:
        out=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/{smp}.index_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/{smp}.index_bam.stderr"
    threads: THREADS
    shell: """
           samtools index {input} 2> {log.err} 1> {log.out}
           """



rule DCC:
    input:
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/joint_mapping/{sample}.Aligned.sortedByCoord.out.bam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/joint_mapping/{sample}.Aligned.sortedByCoord.out.bam.bai", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/joint_mapping/{sample}.Chimeric.out.junction", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/mate1_mapping/{sample}.Chimeric.out.junction", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/mate2_mapping/{sample}.Chimeric.out.junction", sample = SAMPLES_2)
    output:
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/CircCoordinates",
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/CircRNACount",
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/CircSkipJunctions",
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/LinearCount"
    log:
        out=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/DCC.stdout",
        err=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/DCC.stderr"
    threads: THREADS
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/DCC
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_DCC_samplesheet_all.pl -d {WORKING_FOLDER}/standard_results/CircRNA/Mapping -s "{SAMPLE_join}" -n joint -o {WORKING_FOLDER}/standard_results/CircRNA/DCC/samplesheet 2> {log.err} 1> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_DCC_samplesheet_all.pl -d {WORKING_FOLDER}/standard_results/CircRNA/Mapping -s "{SAMPLE_join}" -n mate1 -o {WORKING_FOLDER}/standard_results/CircRNA/DCC/mate1 2>> {log.err} 1>> {log.out}
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_DCC_samplesheet_all.pl -d {WORKING_FOLDER}/standard_results/CircRNA/Mapping -s "{SAMPLE_join}" -n mate2 -o {WORKING_FOLDER}/standard_results/CircRNA/DCC/mate2 2>> {log.err} 1>> {log.out}
           cd {WORKING_FOLDER}/standard_results/CircRNA/DCC 2>> {log.err} 1>> {log.out}
           DCC @samplesheet -mt1 @mate1 -mt2 @mate2 -O {WORKING_FOLDER}/standard_results/CircRNA/DCC -T {threads} -D {DCC_STRAND} -R {REPEAT_GTF} -an {GTF} -Pi -F -M -Nr 5 3 -fg -G -A {GENOME_FA_FILE} 2>> {log.err} 1>> {log.out}
           """



';



for (my $i=0;$i<@treatment;$i++) {



my $num_control=@{$replicate{$control[$i]}};



my $num_treatment=@{$replicate{$treatment[$i]}};



print Profile "rule circRNA_".$control[$i]."_vs_".$treatment[$i].":\n";
print Profile "    input:\n";
print Profile "        expand(WORKING_FOLDER + \"/standard_results/CircRNA/Mapping/{sample}/joint_mapping/{sample}.Aligned.sortedByCoord.out.bam\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "        expand(WORKING_FOLDER + \"/standard_results/CircRNA/Mapping/{sample}/joint_mapping/{sample}.Aligned.sortedByCoord.out.bam.bai\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "        expand(WORKING_FOLDER + \"/standard_results/CircRNA/Mapping/{sample}/joint_mapping/{sample}.Chimeric.out.junction\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "        expand(WORKING_FOLDER + \"/standard_results/CircRNA/Mapping/{sample}/mate1_mapping/{sample}.Chimeric.out.junction\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "        expand(WORKING_FOLDER + \"/standard_results/CircRNA/Mapping/{sample}/mate2_mapping/{sample}.Chimeric.out.junction\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i].")\n";
print Profile "    output:\n";
print Profile "        WORKING_FOLDER + \"/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/CircCoordinates\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/CircRNACount\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/CircSkipJunctions\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/LinearCount\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/\" + PROJECT_NAME + \".STAR_DCC_circTest.".$control[$i]."_vs_".$treatment[$i].".circTest.csv\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/\" + PROJECT_NAME + \".STAR_DCC_circTest.".$control[$i]."_vs_".$treatment[$i].".circTest.txt\"\n";
print Profile "    log:\n";
print Profile "        out=WORKING_FOLDER + \"/standard_results/CircRNA/Snakemake_logs/circRNA_".$control[$i]."_vs_".$treatment[$i].".stdout\",\n";
print Profile "        err=WORKING_FOLDER + \"/standard_results/CircRNA/Snakemake_logs/circRNA_".$control[$i]."_vs_".$treatment[$i].".stderr\"\n";
print Profile "    threads: THREADS\n";
print Profile "    shell: \"\"\"\n";
print Profile "           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."\n";
print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_DCC_samplesheet.pl -d {WORKING_FOLDER}/standard_results/CircRNA/Mapping -c \"{CLASS_".$control[$i]."_2}\" -t \"{CLASS_".$treatment[$i]."_2}\" -n joint -o {WORKING_FOLDER}/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/samplesheet\n";
print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_DCC_samplesheet.pl -d {WORKING_FOLDER}/standard_results/CircRNA/Mapping -c \"{CLASS_".$control[$i]."_2}\" -t \"{CLASS_".$treatment[$i]."_2}\" -n mate1 -o {WORKING_FOLDER}/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/mate1\n";
print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_DCC_samplesheet.pl -d {WORKING_FOLDER}/standard_results/CircRNA/Mapping -c \"{CLASS_".$control[$i]."_2}\" -t \"{CLASS_".$treatment[$i]."_2}\" -n mate2 -o {WORKING_FOLDER}/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/mate2\n";
print Profile "           cd {WORKING_FOLDER}/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."\n";
print Profile "           DCC \@samplesheet -mt1 \@mate1 -mt2 \@mate2 -O {WORKING_FOLDER}/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]." -T {threads} -D {DCC_STRAND} -R {REPEAT_GTF} -an {GTF} -Pi -F -M -Nr 5 3 -fg -G -A {GENOME_FA_FILE}\n";
print Profile "           R --slave --vanilla --args {WORKING_FOLDER}/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]." ".$num_control." ".$num_treatment." ".$control[$i]." ".$treatment[$i]." {PROJECT_NAME} < {HPPRNA_SOFTWARE_FOLDER}/my_R/circTest.R\n";
print Profile "           \"\"\"\n";



}#for (my $i=0;$i<@treatment;$i++) {



}#if ($Type eq "Paired-End") {



##################################Paired-End end###################################




##################################Single-End start###################################



if ($Type eq "Single-End") {



print Profile '
rule STAR_index:
    input:
        GENOME_FA_FILE
    output:
        WORKING_FOLDER + "/standard_results/CircRNA/Index/SA",
        WORKING_FOLDER + "/standard_results/CircRNA/Index/SAindex"
    log:
        out=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/STAR_index.stdout",
        err=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/STAR_index.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/Mapping
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/CircTest
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/Index
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/Snakemake_logs
           STAR --runMode genomeGenerate --genomeDir {WORKING_FOLDER}/standard_results/CircRNA/Index --genomeFastaFiles {input} --runThreadN 1 2> {log.err} 1> {log.out}
           """



rule star_joint_mapping:
    input:
        fastq=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq",
        index_1=WORKING_FOLDER + "/standard_results/CircRNA/Index/SA",
        index_2=WORKING_FOLDER + "/standard_results/CircRNA/Index/SAindex"
    output:
        WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{smp}/single_mapping/{smp}.Aligned.sortedByCoord.out.bam",
        WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{smp}/single_mapping/{smp}.Chimeric.out.junction"
    log:
        out=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/{smp}.star_joint_mapping.stdout",
        err=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/{smp}.star_joint_mapping.stderr"
    threads: THREADS
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/Mapping/{wildcards.smp}
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/Mapping/{wildcards.smp}/single_mapping
           STAR --runThreadN {threads} --genomeDir {WORKING_FOLDER}/standard_results/CircRNA/Index --outSAMtype BAM SortedByCoordinate --readFilesIn {input.fastq} --outFileNamePrefix {WORKING_FOLDER}/standard_results/CircRNA/Mapping/{wildcards.smp}/single_mapping/{wildcards.smp}. --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 2> {log.err} 1> {log.out}
           """



rule index_bam:
    input:
        WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{smp}/single_mapping/{smp}.Aligned.sortedByCoord.out.bam"
    output:
        WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{smp}/single_mapping/{smp}.Aligned.sortedByCoord.out.bam.bai"
    log:
        out=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/{smp}.index_bam.stdout",
        err=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/{smp}.index_bam.stderr"
    threads: THREADS
    shell: """
           samtools index {input} 2> {log.err} 1> {log.out}
           """



rule DCC:
    input:
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/single_mapping/{sample}.Aligned.sortedByCoord.out.bam", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/single_mapping/{sample}.Aligned.sortedByCoord.out.bam.bai", sample = SAMPLES_2),
        expand(WORKING_FOLDER + "/standard_results/CircRNA/Mapping/{sample}/single_mapping/{sample}.Chimeric.out.junction", sample = SAMPLES_2)
    output:
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/CircCoordinates",
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/CircRNACount",
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/CircSkipJunctions",
        WORKING_FOLDER + "/standard_results/CircRNA/DCC/LinearCount"
    log:
        out=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/DCC.stdout",
        err=WORKING_FOLDER + "/standard_results/CircRNA/Snakemake_logs/DCC.stderr"
    threads: THREADS
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/DCC
           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_DCC_samplesheet_all.pl -d {WORKING_FOLDER}/standard_results/CircRNA/Mapping -s "{SAMPLE_join}" -n single -o {WORKING_FOLDER}/standard_results/CircRNA/DCC/samplesheet 2> {log.err} 1> {log.out}
           cd {WORKING_FOLDER}/standard_results/CircRNA/DCC 2>> {log.err} 1>> {log.out}
           DCC @samplesheet -O {WORKING_FOLDER}/standard_results/CircRNA/DCC -T {threads} -D {DCC_STRAND} -R {REPEAT_GTF} -an {GTF} -Pi -F -M -Nr 5 3 -fg -G -A {GENOME_FA_FILE} 2>> {log.err} 1>> {log.out}
           """



';



for (my $i=0;$i<@treatment;$i++) {



my $num_control=@{$replicate{$control[$i]}};



my $num_treatment=@{$replicate{$treatment[$i]}};



print Profile "rule circRNA_".$control[$i]."_vs_".$treatment[$i].":\n";
print Profile "    input:\n";
print Profile "        expand(WORKING_FOLDER + \"/standard_results/CircRNA/Mapping/{sample}/single_mapping/{sample}.Aligned.sortedByCoord.out.bam\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "        expand(WORKING_FOLDER + \"/standard_results/CircRNA/Mapping/{sample}/single_mapping/{sample}.Aligned.sortedByCoord.out.bam.bai\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i]."),\n";
print Profile "        expand(WORKING_FOLDER + \"/standard_results/CircRNA/Mapping/{sample}/single_mapping/{sample}.Chimeric.out.junction\", sample = CLASS_".$control[$i]."_vs_".$treatment[$i].")\n";
print Profile "    output:\n";
print Profile "        WORKING_FOLDER + \"/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/CircCoordinates\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/CircRNACount\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/CircSkipJunctions\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/LinearCount\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/\" + PROJECT_NAME + \".STAR_DCC_circTest.".$control[$i]."_vs_".$treatment[$i].".circTest.csv\",\n";
print Profile "        WORKING_FOLDER + \"/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/\" + PROJECT_NAME + \".STAR_DCC_circTest.".$control[$i]."_vs_".$treatment[$i].".circTest.txt\"\n";
print Profile "    log:\n";
print Profile "        out=WORKING_FOLDER + \"/standard_results/CircRNA/Snakemake_logs/circRNA_".$control[$i]."_vs_".$treatment[$i].".stdout\",\n";
print Profile "        err=WORKING_FOLDER + \"/standard_results/CircRNA/Snakemake_logs/circRNA_".$control[$i]."_vs_".$treatment[$i].".stderr\"\n";
print Profile "    threads: THREADS\n";
print Profile "    shell: \"\"\"\n";
print Profile "           mkdir -p {WORKING_FOLDER}/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."\n";
print Profile "           perl {HPPRNA_SOFTWARE_FOLDER}/my_perl/generate_DCC_samplesheet.pl -d {WORKING_FOLDER}/standard_results/CircRNA/Mapping -c \"{CLASS_".$control[$i]."_2}\" -t \"{CLASS_".$treatment[$i]."_2}\" -n single -o {WORKING_FOLDER}/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."/samplesheet\n";
print Profile "           cd {WORKING_FOLDER}/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]."\n";
print Profile "           DCC \@samplesheet -O {WORKING_FOLDER}/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]." -T {threads} -D {DCC_STRAND} -R {REPEAT_GTF} -an {GTF} -Pi -F -M -Nr 5 3 -fg -G -A {GENOME_FA_FILE}\n";
print Profile "           R --slave --vanilla --args {WORKING_FOLDER}/standard_results/CircRNA/CircTest/".$control[$i]."_vs_".$treatment[$i]." ".$num_control." ".$num_treatment." ".$control[$i]." ".$treatment[$i]." {PROJECT_NAME} < {HPPRNA_SOFTWARE_FOLDER}/my_R/circTest.R\n";
print Profile "           \"\"\"\n";



}#for (my $i=0;$i<@treatment;$i++) {



}#if ($Type eq "Single-End") {



##################################Single-End end###################################



###############################################circRNA workflow end############################################################



}#if ($Analysis_type eq "circRNA") {



###############################################DNA variation start############################################################



############################SNP start#################################



if ($SNP eq "Yes") {



if ($Type eq "Paired-End") {



print Profile '
rule prepare_snp_index:
    input:
        HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + "/genome.fa"
    output:
        fa=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.fa",
        dict=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.dict",
        fai=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.fa.fai",
        index_1=WORKING_FOLDER + "/standard_results/SNP/Index/SA",
        index_2=WORKING_FOLDER + "/standard_results/SNP/Index/SAindex"
    log:
        out=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/prepare_snp_index.stdout",
        err=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/prepare_snp_index.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/Index
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/Snakemake_logs
           cp {input} {output.fa} 2> {log.err} 1> {log.out}
           java -jar {HPPRNA_SOFTWARE_FOLDER}/picard-tools-2.3.0/picard.jar CreateSequenceDictionary R={output.fa} O={output.dict} 2>> {log.err} 1>> {log.out}
           samtools faidx {output.fa} 2>> {log.err} 1>> {log.out}
           STAR --runMode genomeGenerate --genomeDir {WORKING_FOLDER}/standard_results/SNP/Index --genomeFastaFiles {output.fa} --runThreadN 1 2>> {log.err} 1>> {log.out}
           """



rule star_for_snp_1pass:
    input:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq",
        index_1=WORKING_FOLDER + "/standard_results/SNP/Index/SA",
        index_2=WORKING_FOLDER + "/standard_results/SNP/Index/SAindex"
    output:
        WORKING_FOLDER + "/standard_results/SNP/{smp}/1pass/{smp}.Aligned.out.sam",
        WORKING_FOLDER + "/standard_results/SNP/{smp}/1pass/{smp}.SJ.out.tab"
    log:
        out=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.star_for_snp_1pass.stdout",
        err=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.star_for_snp_1pass.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/1pass
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/Index_2pass
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/2pass
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/BAM
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/vcf
           STAR --genomeDir {WORKING_FOLDER}/standard_results/SNP/Index --readFilesIn {input.R1} {input.R2} --outFileNamePrefix {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/1pass/{wildcards.smp}. --runThreadN 1 2> {log.err} 1> {log.out}
           """



rule star_for_snp_index_2pass:
    input:
        fa=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.fa",
        dict=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.dict",
        fai=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.fa.fai",
        tab=WORKING_FOLDER + "/standard_results/SNP/{smp}/1pass/{smp}.SJ.out.tab"
    output:
        WORKING_FOLDER + "/standard_results/SNP/{smp}/Index_2pass/SA",
        WORKING_FOLDER + "/standard_results/SNP/{smp}/Index_2pass/SAindex"
    log:
        out=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.star_for_snp_index_2pass.stdout",
        err=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.star_for_snp_index_2pass.stderr"
    shell: """
           STAR --runMode genomeGenerate --genomeDir {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/Index_2pass --genomeFastaFiles {input.fa} --sjdbFileChrStartEnd {input.tab} --sjdbOverhang 75 --runThreadN 1 2> {log.err} 1> {log.out}
           """



rule star_for_snp_2pass:
    input:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq",
        index_1=WORKING_FOLDER + "/standard_results/SNP/{smp}/Index_2pass/SA",
        index_2=WORKING_FOLDER + "/standard_results/SNP/{smp}/Index_2pass/SAindex"
    output:
        WORKING_FOLDER + "/standard_results/SNP/{smp}/2pass/{smp}.Aligned.out.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.star_for_snp_2pass.stdout",
        err=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.star_for_snp_2pass.stderr"
    shell: """
           STAR --genomeDir {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/Index_2pass --readFilesIn {input.R1} {input.R2} --outFileNamePrefix {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/2pass/{wildcards.smp}. --runThreadN 1 2> {log.err} 1> {log.out}
           """



rule picard:
    input:
        WORKING_FOLDER + "/standard_results/SNP/{smp}/2pass/{smp}.Aligned.out.sam"
    output:
        sort_bam=WORKING_FOLDER + "/standard_results/SNP/{smp}/BAM/{smp}.rg_added_sorted.bam",
        dedupped_dam=WORKING_FOLDER + "/standard_results/SNP/{smp}/BAM/{smp}.dedupped.bam",
        metrics=WORKING_FOLDER + "/standard_results/SNP/{smp}/BAM/{smp}.output.metrics"
    log:
        out=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.picard.stdout",
        err=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.picard.stderr"
    shell: """
           java -jar {HPPRNA_SOFTWARE_FOLDER}/picard-tools-2.3.0/picard.jar AddOrReplaceReadGroups I={input} O={output.sort_bam} SO=coordinate RGID={wildcards.smp} RGLB={wildcards.smp} RGPL=platform RGPU=machine RGSM={wildcards.smp} 2> {log.err} 1> {log.out}
           java -jar {HPPRNA_SOFTWARE_FOLDER}/picard-tools-2.3.0/picard.jar MarkDuplicates I={output.sort_bam} O={output.dedupped_dam} CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={output.metrics} 2>> {log.err} 1>> {log.out}
           """



rule GATK:
    input:
        fa=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.fa",
        bam=WORKING_FOLDER + "/standard_results/SNP/{smp}/BAM/{smp}.dedupped.bam"
    output:
        bam=WORKING_FOLDER + "/standard_results/SNP/{smp}/BAM/{smp}.split.bam",
        vcf_call=WORKING_FOLDER + "/standard_results/SNP/{smp}/vcf/{smp}.call.vcf",
        vcf_filter=WORKING_FOLDER + "/standard_results/SNP/{smp}/vcf/{smp}.filter.vcf"
    log:
        out=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.GATK.stdout",
        err=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.GATK.stderr"
    threads: THREADS
    shell: """
           java -jar {GATK_FOLDER}/GenomeAnalysisTK.jar -T SplitNCigarReads -R {input.fa} -I {input.bam} -o {output.bam} -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 2> {log.err} 1> {log.out}
           java -jar {GATK_FOLDER}/GenomeAnalysisTK.jar -T HaplotypeCaller -R {input.fa} -I {output.bam} -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o {output.vcf_call} 2>> {log.err} 1>> {log.out}
           java -jar {GATK_FOLDER}/GenomeAnalysisTK.jar -T VariantFiltration -R {input.fa} -V {output.vcf_call} -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o {output.vcf_filter} 2>> {log.err} 1>> {log.out}
           """



';



}#if ($Type eq "Paired-End") {



if ($Type eq "Single-End") {



print Profile '
rule prepare_snp_index:
    input:
        HPPRNA_SOFTWARE_FOLDER + "/hppRNA_genome/" + ASSEMBLY_VERSION + "/genome.fa"
    output:
        fa=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.fa",
        dict=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.dict",
        fai=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.fa.fai",
        index_1=WORKING_FOLDER + "/standard_results/SNP/Index/SA",
        index_2=WORKING_FOLDER + "/standard_results/SNP/Index/SAindex"
    log:
        out=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/prepare_snp_index.stdout",
        err=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/prepare_snp_index.stderr"
    threads: THREADS
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/Index
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/Snakemake_logs
           cp {input} {output.fa} 2> {log.err} 1> {log.out}
           java -jar {HPPRNA_SOFTWARE_FOLDER}/picard-tools-2.3.0/picard.jar CreateSequenceDictionary R={output.fa} O={output.dict} 2>> {log.err} 1>> {log.out}
           samtools faidx {output.fa} 2>> {log.err} 1>> {log.out}
           STAR --runMode genomeGenerate --genomeDir {WORKING_FOLDER}/standard_results/SNP/Index --genomeFastaFiles {output.fa} --runThreadN {threads} 2>> {log.err} 1>> {log.out}
           """



rule star_for_snp:
    input:
        fastq=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq",
        fa=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.fa",
        dict=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.dict",
        fai=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.fa.fai",
        index_1=WORKING_FOLDER + "/standard_results/SNP/Index/SA",
        index_2=WORKING_FOLDER + "/standard_results/SNP/Index/SAindex"
    output:
        WORKING_FOLDER + "/standard_results/SNP/{smp}/1pass/{smp}.Aligned.out.sam",
        WORKING_FOLDER + "/standard_results/SNP/{smp}/Index_2pass/SA",
        WORKING_FOLDER + "/standard_results/SNP/{smp}/Index_2pass/SAindex",
        WORKING_FOLDER + "/standard_results/SNP/{smp}/2pass/{smp}.Aligned.out.sam"
    log:
        out=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.star_for_snp.stdout",
        err=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.star_for_snp.stderr"
    threads: THREADS
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/1pass
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/Index_2pass
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/2pass
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/BAM
           mkdir -p {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/vcf
           STAR --genomeDir {WORKING_FOLDER}/standard_results/SNP/Index --readFilesIn {input.fastq} --outFileNamePrefix {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/1pass/{wildcards.smp}. --runThreadN {threads} 2> {log.err} 1> {log.out}
           STAR --runMode genomeGenerate --genomeDir {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/Index_2pass --genomeFastaFiles {input.fa} --sjdbFileChrStartEnd {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/1pass/{wildcards.smp}.SJ.out.tab --sjdbOverhang 75 --runThreadN {threads} 2>> {log.err} 1>> {log.out}
           STAR --genomeDir {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/Index_2pass --readFilesIn {input.fastq} --outFileNamePrefix {WORKING_FOLDER}/standard_results/SNP/{wildcards.smp}/2pass/{wildcards.smp}. --runThreadN {threads} 2>> {log.err} 1>> {log.out}
           """



rule picard:
    input:
        WORKING_FOLDER + "/standard_results/SNP/{smp}/2pass/{smp}.Aligned.out.sam"
    output:
        sort_bam=WORKING_FOLDER + "/standard_results/SNP/{smp}/BAM/{smp}.rg_added_sorted.bam",
        dedupped_dam=WORKING_FOLDER + "/standard_results/SNP/{smp}/BAM/{smp}.dedupped.bam",
        metrics=WORKING_FOLDER + "/standard_results/SNP/{smp}/BAM/{smp}.output.metrics"
    log:
        out=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.picard.stdout",
        err=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.picard.stderr"
    shell: """
           java -jar {HPPRNA_SOFTWARE_FOLDER}/picard-tools-2.3.0/picard.jar AddOrReplaceReadGroups I={input} O={output.sort_bam} SO=coordinate RGID={wildcards.smp} RGLB={wildcards.smp} RGPL=platform RGPU=machine RGSM={wildcards.smp} 2> {log.err} 1> {log.out}
           java -jar {HPPRNA_SOFTWARE_FOLDER}/picard-tools-2.3.0/picard.jar MarkDuplicates I={output.sort_bam} O={output.dedupped_dam} CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={output.metrics} 2>> {log.err} 1>> {log.out}
           """



rule GATK:
    input:
        fa=WORKING_FOLDER + "/standard_results/SNP/Index/" + ASSEMBLY_VERSION + ".genome.fa",
        bam=WORKING_FOLDER + "/standard_results/SNP/{smp}/BAM/{smp}.dedupped.bam"
    output:
        bam=WORKING_FOLDER + "/standard_results/SNP/{smp}/BAM/{smp}.split.bam",
        vcf_call=WORKING_FOLDER + "/standard_results/SNP/{smp}/vcf/{smp}.call.vcf",
        vcf_filter=WORKING_FOLDER + "/standard_results/SNP/{smp}/vcf/{smp}.filter.vcf"
    log:
        out=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.GATK.stdout",
        err=WORKING_FOLDER + "/standard_results/SNP/Snakemake_logs/{smp}.GATK.stderr"
    threads: THREADS
    shell: """
           java -jar {HPPRNA_SOFTWARE_FOLDER}/GATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R {input.fa} -I {input.bam} -o {output.bam} -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 2> {log.err} 1> {log.out}
           java -jar {HPPRNA_SOFTWARE_FOLDER}/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R {input.fa} -I {output.bam} -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o {output.vcf_call} 2>> {log.err} 1>> {log.out}
           java -jar {HPPRNA_SOFTWARE_FOLDER}/GATK/GenomeAnalysisTK.jar -T VariantFiltration -R {input.fa} -V {output.vcf_call} -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o {output.vcf_filter} 2>> {log.err} 1>> {log.out}
           """



';



}#if ($Type eq "Single-End") {



}#if ($SNP eq "Yes") {



############################SNP end#################################



############################Fusion gene start#################################



if ($Fusion eq "Yes") {



if ($Type eq "Paired-End") {



print Profile '
rule prepare_fusion_fastq:
    input:
        R1=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}_2.fastq"
    output:
        R1=WORKING_FOLDER + "/standard_results/Fusion_gene/Raw_data/{smp}/{smp}_1.fastq",
        R2=WORKING_FOLDER + "/standard_results/Fusion_gene/Raw_data/{smp}/{smp}_2.fastq"
    log:
        out=WORKING_FOLDER + "/standard_results/Fusion_gene/Snakemake_logs/{smp}.prepare_fusion_fastq.stdout",
        err=WORKING_FOLDER + "/standard_results/Fusion_gene/Snakemake_logs/{smp}.prepare_fusion_fastq.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene/Raw_data
           mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene/Raw_data/{wildcards.smp}
           mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene/Snakemake_logs
           cp {input.R1} {output.R1} 2> {log.err} 1> {log.out}
           cp {input.R2} {output.R2} 2>> {log.err} 1>> {log.out}
           """



rule fusioncatcher:
    input:
        WORKING_FOLDER + "/standard_results/Fusion_gene/Raw_data/{smp}/{smp}_1.fastq",
        WORKING_FOLDER + "/standard_results/Fusion_gene/Raw_data/{smp}/{smp}_2.fastq"
    output:
        WORKING_FOLDER + "/standard_results/Fusion_gene/Fusion_result/{smp}/summary_candidate_fusions.txt",
        WORKING_FOLDER + "/standard_results/Fusion_gene/Fusion_result/{smp}/final-list_candidate-fusion-genes.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/Fusion_gene/Snakemake_logs/{smp}.fusioncatcher.stdout",
        err=WORKING_FOLDER + "/standard_results/Fusion_gene/Snakemake_logs/{smp}.fusioncatcher.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene/Fusion_result
           mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene/Fusion_result/{wildcards.smp}
           fusioncatcher -d {FUSIONCATCHER_DATA_FOLDER} -i {WORKING_FOLDER}/standard_results/Fusion_gene/Raw_data/{wildcards.smp}/ -o {WORKING_FOLDER}/standard_results/Fusion_gene/Fusion_result/{wildcards.smp}/ 2> {log.err} 1> {log.out}
           """



';



}#if ($Type eq "Paired-End") {



if ($Type eq "Single-End") {



print Profile '
rule prepare_fusion_fastq:
    input:
        WORKING_FOLDER + "/standard_results/Processed_FASTQ/{smp}.adapt.qc.fastq"
    output:
        WORKING_FOLDER + "/standard_results/Fusion_gene/Raw_data/{smp}/{smp}.adapt.qc.fastq"
    log:
        out=WORKING_FOLDER + "/standard_results/Fusion_gene/Snakemake_logs/{smp}.prepare_fusion_fastq.stdout",
        err=WORKING_FOLDER + "/standard_results/Fusion_gene/Snakemake_logs/{smp}.prepare_fusion_fastq.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene/Raw_data
           mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene/Raw_data/{wildcards.smp}
           mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene/Snakemake_logs
           cp {input} {output} 2> {log.err} 1> {log.out}
           """



rule fusioncatcher:
    input:
        WORKING_FOLDER + "/standard_results/Fusion_gene/Raw_data/{smp}/{smp}.adapt.qc.fastq"
    output:
        WORKING_FOLDER + "/standard_results/Fusion_gene/Fusion_result/{smp}/summary_candidate_fusions.txt",
        WORKING_FOLDER + "/standard_results/Fusion_gene/Fusion_result/{smp}/final-list_candidate-fusion-genes.txt"
    log:
        out=WORKING_FOLDER + "/standard_results/Fusion_gene/Snakemake_logs/{smp}.fusioncatcher.stdout",
        err=WORKING_FOLDER + "/standard_results/Fusion_gene/Snakemake_logs/{smp}.fusioncatcher.stderr"
    shell: """
           mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene/Fusion_result
           mkdir -p {WORKING_FOLDER}/standard_results/Fusion_gene/Fusion_result/{wildcards.smp}
           fusioncatcher --single-end -d {FUSIONCATCHER_DATA_FOLDER} -i {WORKING_FOLDER}/standard_results/Fusion_gene/Raw_data/{wildcards.smp}/ -o {WORKING_FOLDER}/standard_results/Fusion_gene/Fusion_result/{wildcards.smp}/ 2> {log.err} 1> {log.out}
           """



';



}#if ($Type eq "Single-End") {



}#if ($Fusion eq "Yes") {



############################Fusion gene end#################################



###############################################DNA variation end############################################################
#####################################################################################Snakemake end##########################################################################################



}



close Profile;


