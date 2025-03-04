#!/bin/bash

#Script continuing from 1_DataProcessing_MLphylogenetics.sh, this one is for SNP calling steps 
#Using Clade1 of genus Acropora data (Scleractinia: Acroporidae)
#Author: Hanaka Mera   mailto:hanaka.mera[@]my.jcu.edu.au
    #Tutorials and papers followed (in this script):
    #Erickson, K.L., Pentico, A., Quattrini, A.M., McFadden, C.S., 2021. New approaches to species delimitation and population structure of anthozoans: Two case studies of octocorals using ultraconserved elements and exons. Mol Ecol Resour 21, 78–92. https://doi.org/10.1111/1755-0998.13241
    #https://gatk.broadinstitute.org/hc/en-us/articles/360035532412-Can-t-use-VQSR-on-non-model-organism-or-small-dataset
    #https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
    #And thanks to many discussions with Lauriane Baraf

#Last edited: July 2023

###################################################################

## SNP calling with GATK4

### 0. General overview of all the commands

#The base idea is that for non-model organisms, we don't have a good reference genome that we can map samples against and extract good quality SNPs. So in the below steps, we first create crude SNPs/INDELs lists from our own crude samples (3_GVCFgenotype1-filter1), recalibrate those lists using the same samples (4_recal1_variantcall2+5_GVCFgenotype2-filter2), then do base score correction on the crude samples in preparation for a 'true' SNP calling (6_recal2_variantcall3) then true SNP call  using the refined SNPs/INDELs lists (7_GVCFgenotype3-filter3). And if the lists are still not refined (i.e. AnalyzeCovariates not converging), you keep doing the recalibration steps until it is refined.

#Packages (if in conda environment, you might want to do conda activate --stack package. Use singularity containers where possible)
    #phyluce-1.7.1
    #bwa-0.7.17
    #samtools-1.12
    #gatk4-4.3.0.0
    #R-4.0.3
    #bcftools-1.16
    #vcftools-0.1.16

##1. Data preparation
    ###Reference index
    ###Generate dictionary
    ###Map reads to reference, sort, mark&remove duplicates, index bam (ARRAY JOB #1 - 1_bwa_mapping_array.pbs)
##2. Initial round of variant calling with HaplotypeCaller (ARRAY JOB #2 - 2_variantcall1_array.pbs)
##3. Combine GVCF, joint genotype, extract SNPs/INDELs from the call set and hard filter (3_GVCFgenotype1-filter1.pbs)
##4. 1st Base recalibration and variant call on recalibrated bams (ARRAY JOB #3 - 4_recal1_variantcall2_array.pbs)
##5. Combine GVCF, joint genotype, extract SNPs/INDELs from the call set and hard filter (5_GVCFgenotype2-filter2.pbs)
##6. 2nd Base recalibration and variant call on recalibrated bams (ARRAY JOB #4 - 6_recal2_variantcall3_array.pbs)
    #At this step, check the before/after plots for convergence with AnalyzeCovariates. If converged, the output from variantcall3 is your final.
##7. Combine GVCF, joint genotype, extract SNPs/INDELs from the call set and hard filter (7_GVCFgenotype3-filter3.pbs)
##8. Get stats and prepare final vcf for downstream analysis

###################################################################

### 1. Data preparation

mkdir -p /scratch/snp_calling/pbsscripts
mkdir -p /scratch/snp_calling/Clade1_refC99/{log,1_dataprep,2_out_variantcall1,3_out_GVCF1filter1,4_out_recal1variantcall2,5_out_GVCF2filter2,6_out_recal2variantcall3,7_out_GVCF3filter3,8_out_fordownstream}

#Find sample with most CONTIGS to use as reference individual. This info can be obtained from the exploded fastas summary stats.
sort -t ',' -k 2,2 -nfr /home/4_match/taxon-filt/results_explodedfastas.txt  #C99
cd /scratch/snp_calling/Clade1_refC99/1_dataprep
#Create config file
printf '[ref]\n''C99\n' > ref.conf
#Run phyluce to create list of all loci present in the reference individual
singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_assembly_get_match_counts \
    --locus-db /home/4_match/search-results/probe.matches.sqlite \
    --taxon-list-config ref.conf \
    --taxon-group ref \
    --output ref-loci.conf \
    --log-path ../log
#Run phyluce to create fasta file of loci present in the reference individual
singularity run /fast/tmp/containers/phyluce-1.7.1.sif phyluce_assembly_get_fastas_from_match_counts \
    --contigs /home/3_correction/correctedcontigs \
    --locus-db /home/4_match/search-results/probe.matches.sqlite \
    --match-count-output ref-loci.conf \
    --output ref-loci.fasta \
    --log-path ../log

#### Generate dictionary

singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx2g' CreateSequenceDictionary \
    --REFERENCE ref-loci.fasta --OUTPUT ref-loci.dict
singularity run /fast/tmp/containers/samtools-1.16.1.sif samtools faidx ref-loci.fasta
#Create a list of directories which have cleaned reads (both READ1 and READ2)
ls -d /home/2_clean/C* > dir_clean.list
ls -d /home/3_subsampled/* >> dir_clean.list


#### Map reads to reference, sort, mark&remove duplicates, index bam
    #ARRAY JOB #1 - 1_bwa_mapping_array.pbs
nano /scratch/snp_calling/pbsscripts/1_bwa_mapping_array.pbs
#####BEGIN PBS
#!/bin/bash
#PBS -j oe
#PBS -m ae
#PBS -N 1_bwa_array
#PBS -l select=1:ncpus=5:mem=50gb
#PBS -l walltime=5:00:00
#PBS -M my@email
#PBS -J 0-179

shopt -s expand_aliases
#Requires: bwa, samtools, GATK4

cd /scratch/snp_calling/Clade1_refC99/1_dataprep
id=$PBS_ARRAY_INDEX
clean_reads=($(cat dir_clean.list))
sample=$(echo ${clean_reads[$id]} | sed 's/.$//; s/.*\///')

#R1/R2 files
r1=$(echo "${clean_reads[$id]}$sample-READ1.fastq.gz")
r2=$(echo "${clean_reads[$id]}$sample-READ2.fastq.gz")

#Map reads with algorithm mem for illumina reads 70bp-1Mb
    #-t numThreads, -B Mismatch penalty (seq error rate), -M Mark shorter split hits as secondary, -R read group header line, ref prefix
singularity run /fast/tmp/containers/bwa-0.7.17.sif bwa mem -t 5 -B 10 \
    -M -R "@RG\tID:$sample\tSM:$sample\tPL:Illumina" ref-loci $r1 $r2 > $sample.ref_pair.sam

#Sort reads
singularity run /fast/tmp/containers/samtools-1.16.1.sif samtools view -@ 5 -bS $sample.ref_pair.sam | \
    singularity run /fast/tmp/containers/samtools-1.16.1.sif samtools sort -@ 5 -m 30000000000 -o $sample.ref_pair_sorted.bam

#Mark and remove duplicates
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk --java-options '-Xmx2g' MarkDuplicates \
    -I $sample.ref_pair_sorted.bam \
    -O $sample.ref_all_dedup.bam \
    -M $sample.ref_all_dedup_metricsfile \
    --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 250 \
    --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT
    #--ASSUME_SORTED true is depreciated in GATK4
    #--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP is HPC dependent. check by ulimit -n and set the number to lower than that.

#Index bam file
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk --java-options '-Xmx2g' BuildBamIndex \
    -I $sample.ref_all_dedup.bam
singularity run /fast/tmp/containers/samtools-1.16.1.sif samtools flagstat \
    -@ 5 $sample.ref_all_dedup.bam > $sample.ref_all_dedup_stats.txt
#Remove intermediate files
rm $sample.ref_pair*
#####END PBS

cd /scratch/snp_calling/Clade1_refC99/log
qsub ../../pbsscripts/1_bwa_mapping_array.pbs

#Get stats
cd /scratch/snp_calling/Clade1_refC99/1_dataprep
#MarkDuplicates stats
echo 'LIBRARY READ_PAIRS_EXAMINED SECONDARY_OR_SUPPLEMENTARY_RDS UNMAPPED_READS UNPAIRED_READ_DUPLICATES READ_PAIR_DUPLICATES READ_PAIR_OPTICAL_DUPLICATES PERCENT_DUPLICATION ESTIMATED_LIBRARY_SIZE' > markduplicates_all_stats.txt
for file in *dedup_metricsfile; do
    sample=$(echo $file | cut -d. -f1)
    #Print the relevant line in the metricsfile with sample name
    stats=$(sed -n '8p' $file | cut -f3-11)
    echo $sample $stats >> markduplicates_all_stats.txt; done

#samtools stats
echo "LIBRARY,Total_QC_passed_reads,secondary,supplementary,duplicates,MAPPED,PAIRED,read1,read2,properly_paired,with_itself_mate_paired,singletons,mate_mapped_different_chr,mate_mapped_different_chrQ5" > samtools_stats.txt
for file in *dedup_stats.txt; do
    sample=$(echo $file | cut -d. -f1)
    stats=$(tr '\n' , < $file)
    #Print the relevant line in the stats.txt file with sample name
    echo $sample','$stats >> samtools_stats.txt; done
    #The duplicates number should be 0 because we removed it in gatk MarkDuplicates

#############################################################

### 2. Initial round of variant calling with HaplotypeCaller
    #ARRAY JOB #2 - 2_variantcall1_array.pbs, 160samples, ~5hrs

nano ../../pbsscripts/2_variantcall1_array.pbs 
#####BEGIN PBS
#!/bin/bash
#PBS -j oe
#PBS -m ae
#PBS -N 2_variantcall1_array
#PBS -l select=1:ncpus=5:mem=50gb
#PBS -l walltime=15:00:00
#PBS -M my@email
#PBS -J 0-179

shopt -s expand_aliases
#Requires: GATK4

cd /scratch/snp_calling/Clade1_refC99/2_out_variantcall1
#Reference taxon
ref='/scratch/snp_calling/Clade1_refC99/1_dataprep/ref-loci.fasta'
#Bams after removing duplicates with gatk
id=$PBS_ARRAY_INDEX
dedup_bams=($(ls /scratch/snp_calling/Clade1_refC99/1_dataprep/*all_dedup.bam))
dedup_bam=${dedup_bams[$id]}
#Create a variable with the dedup_bam. Use path variable to extract the last filename with cut -d
sample=$(echo $dedup_bam | cut -d/ -f7 | cut -d. -f1)
echo $sample

echo "Processing $dedup_bam"
#Execute HaplotypeCaller in GATK for Variant discovery
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk --java-options '-Xmx4g' HaplotypeCaller \
    -R $ref -I $dedup_bam -O ${sample}.g.vcf \
    --contamination-fraction-to-filter 0.0002 \
    --emit-ref-confidence GVCF \
    --min-base-quality-score 20 \
    --phred-scaled-global-read-mismapping-rate 30 \
    --standard-min-confidence-threshold-for-calling 40.0
#####END PBS
qsub ../../pbsscripts/2_variantcall1_array.pbs

#############################################################

### 3. Combine GVCF, joint genotype, extract SNPs/INDELs from the call set and hard filter 
    #3_GVCFgenotype1-filter1.pbs, 177samples ~ 45min

#####BEGIN PBS
#!/bin/bash
#PBS -j oe
#PBS -m ae
#PBS -N 3_GVCFgenotype1-filter1
#PBS -l select=1:ncpus=5:mem=30gb
#PBS -l walltime=2:00:00
#PBS -M my@email

cd /scratch/snp_calling/Clade1_refC99/3_out_GVCF1filter1
#Reference taxon
ref='/scratch/snp_calling/Clade1_refC99/1_dataprep/ref-loci.fasta'

#Get paths for all vcf files and combine into a cohort gvcf
for i in ../2_out_variantcall1/*.g.vcf; do echo " -V $i "; done | tr -d "\n" > gvcf.list
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' CombineGVCFs \
    -R $ref $(cat gvcf.list) -O cohort.g.vcf

#Genotyping with GVCF in the cohort variant files produced by HaplotypeCaller gvcf
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' GenotypeGVCFs \
    -R $ref -V cohort.g.vcf -O genotyped.vcf

#Extract the SNPs from the call set
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SelectVariants \
    -R $ref -V genotyped.vcf -O genotyped_snps.vcf \
    --select-type-to-include SNP
#Extract the indels from the call set
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SelectVariants \
    -R $ref -V genotyped.vcf -O genotyped_indels.vcf \
    --select-type-to-include INDEL

#Hard filter variants according to filter expression: SNPs
    #Designed for hard-filtering using statistics annotations. Better to use VQSR because it uses machine-learning algorithms, but this is only available for large variants and well-curated variant resources e.g. humans, mice. For non-model organisms, we have to use the hard-filtering and combine different criteria (e.g. low quality sequencing reads) and hopefully it doesn’t filter out ‘good’ ones that happened to have bad score for one of the criteria. The output will have a FILTER field to something other than PASS according to the criteria.
#ALSO NOTE: The space in '-filter "xx < 00"' seems to not be a problem for the first 5 filters but is a problem for MQRankSum and ReadPosRankSum. It will still produce warnings because not all variants get MQRS and RPRS scores.
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' VariantFiltration \
    -R $ref -V genotyped_snps.vcf \
    -O genotyped_snps_filtered.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "SOR > 4.0" --filter-name "SOR4" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-filter "MQRankSum<-12.5" --filter-name "MQRS-12.5" \
	-filter "ReadPosRankSum<-8.0" --filter-name "RPRS-8"
#Hard filter variants according to filter expression: INDELs
    #These values are mostly GATK's recommendation
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' VariantFiltration \
    -R $ref -V genotyped_indels.vcf \
    -O genotyped_indels_filtered.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "SOR > 10.0" --filter-name "SOR10" 

#Remove filtered snps from vcf file (this should only keep the ones marked as PASS)
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SelectVariants \
    -R $ref -V genotyped_snps_filtered.vcf \
    -O genotyped_snps_onlyPASS.vcf \
    --exclude-filtered true
#Remove filtered indels from vcf file
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SelectVariants \
    -R $ref -V genotyped_indels_filtered.vcf \
    -O genotyped_indels_onlyPASS.vcf \
    --exclude-filtered true
#####END PBS
qsub ../../pbsscripts/3_GVCFgenotype1-filter1.pbs

#Check outputs (grep -v "^#" to exclude the headers that start with '#')
grep -v "^#" genotyped_snps.vcf | wc -l #186627
grep -v "^#" genotyped_snps_filtered.vcf | wc -l #186627
grep -v "^#" genotyped_snps_onlyPASS.vcf | wc -l #158699
#This does the same thing in gatk:
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk CountVariants -V genotyped_snps_onlyPASS.vcf #158699

#############################################################

### 4. 1st Base recalibration and 2nd variant call on recalibrated bams
    #ARRAY JOB #3 - 4_recal1_variantcall2_array.pbs, 160samples ~5hrs

#####BEGIN PBS
#!/bin/bash
#PBS -j oe
#PBS -m ae
#PBS -N 4_recal1_variantcall2_array
#PBS -l select=1:ncpus=5:mem=50gb
#PBS -l walltime=15:00:00
#PBS -M my@email
#PBS -J 0-179

shopt -s expand_aliases
#Requires: GATK4, R
#Currently singularity `module load R/4.1.2` does not work for this so have to use conda until it becomes unavailable.
module load conda3
source $CONDA_PROF/conda.sh
conda activate R-4.0.3

cd /scratch/snp_calling/Clade1_refC99/4_out_recal1variantcall2
#Reference taxon
ref=/scratch/snp_calling/Clade1_refC99/1_dataprep/ref-loci.fasta
#Bams after removing duplicates with gatk
id=$PBS_ARRAY_INDEX
dedup_bams=($(ls /scratch/snp_calling/Clade1_refC99/1_dataprep/*all_dedup.bam))
dedup_bam=${dedup_bams[$id]}
#Create a variable with the dedup_bam. Use path variable to extract the last filename with cut -d
sample=$(echo $dedup_bam | cut -d/ -f7 | cut -d. -f1)
echo $sample

####Base recalibration 1
echo "Processing $dedup_bam"
#Execute BaseRecalibratorn GATK to get PRE recalibration table and ApplyBQSR to run 1st base recalibration
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx4g' BaseRecalibrator \
    -R $ref -I $dedup_bam -O ${sample}_pre-recal1.table \
    --known-sites ../3_out_GVCF1filter1/genotyped_snps_onlyPASS.vcf \
    --known-sites ../3_out_GVCF1filter1/genotyped_indels_onlyPASS.vcf
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx4g' ApplyBQSR \
    -R $ref -I $dedup_bam --bqsr-recal-file $sample'_pre-recal1.table' \
    -O ${sample}_recal1.bam
#Execute BaseRecalibratorn GATK to get POST recalibration table
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx4g' BaseRecalibrator \
    -R $ref -I ${sample}_recal1.bam -O ${sample}_post-recal1.table \
    --known-sites ../3_out_GVCF1filter1/genotyped_snps_onlyPASS.vcf \
    --known-sites ../3_out_GVCF1filter1/genotyped_indels_onlyPASS.vcf
#Analyze covariates to see whether it converged
    #This requires ggplot2 and gsalib Rpackages to be installed on the same system.
    #wget https://cran.r-project.org/src/contrib/Archive/gsalib/gsalib_2.1.tar.gz
    #Then open R, run `install.packages("/dir/to/gsalib_2.1.tar.gz")`
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx4g' AnalyzeCovariates \
    -before ${sample}_pre-recal1.table \
    -after ${sample}_post-recal1.table \
    -plots ${sample}_AnalyzeCovariates.pdf \
    -csv ${sample}_AnalyzeCovariates.csv
    #If error persists, copy BQSR.R script to local, then edit to 'data <- read.csv(args[1], stringsAsFactors = TRUE)', add the Rscript location to $PATH then run below instead
/home/R/BQSR.R ${sample}_AnalyzeCovariates.csv ${sample}_pre-recal1.table ${sample}_ext_AnalyzeCovariates.pdf

#HaplotypeCaller in GATK for variant discovery on recalibrated bams
echo "Processing ${sample}_recal1.bam"
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx4g' HaplotypeCaller \
    -R $ref -I ${sample}_recal1.bam -O ${sample}_recal1.g.vcf \
    --contamination-fraction-to-filter 0.0002 \
    --emit-ref-confidence GVCF \
    --min-base-quality-score 20 \
    --phred-scaled-global-read-mismapping-rate 30 \
    --standard-min-confidence-threshold-for-calling 40.0
#####END PBS
qsub ../../pbsscripts/4_recal1_variantcall2_array.pbs

#Check $sample_AnalyzeCovariates.pdfs. The before and after plots should look different (and after plots should be better e.g. closer to 0 or linear. see GATK website for how to interpret individual plots). The goal is to do this whole process again so that the before and after plots converge i.e. look the same.

#############################################################

### 5. Combine GVCF, joint genotype, extract SNPs/INDELs from the call set and hard filter 
    #5_GVCFgenotype2-filter2.pbs, 160samples ~ 45min

#####BEGIN PBS
#!/bin/bash
#PBS -j oe
#PBS -m ae
#PBS -N 5_GVCFgenotype2-filter2
#PBS -l select=1:ncpus=5:mem=30gb
#PBS -l walltime=5:00:00
#PBS -M my@email

cd /scratch/snp_calling/Clade1_refC99/5_out_GVCF2filter2
#Reference taxon
ref='/scratch/snp_calling/Clade1_refC99/1_dataprep/ref-loci.fasta'

#Get paths for all vcf files and combine into a cohort gvcf
for i in ../4_out_recal1variantcall2/*_recal1.g.vcf; do echo " -V $i "; done | tr -d "\n" > gvcf_recal1.list
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' CombineGVCFs \
    -R $ref $(cat gvcf_recal1.list) -O cohort_recal1.g.vcf

#Genotyping with GVCF in the cohort variant files produced by HaplotypeCaller gvcf
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' GenotypeGVCFs \
    -R $ref -V cohort_recal1.g.vcf -O genotyped_recal1.vcf

#Extract the SNPs from the call set
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SelectVariants \
    -R $ref -V genotyped_recal1.vcf \
    -O genotyped_recal1_snps.vcf \
    --select-type-to-include SNP
#Extract the indels from the call set
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SelectVariants \
    -R $ref -V genotyped_recal1.vcf \
    -O genotyped_recal1_indels.vcf \
    --select-type-to-include INDEL

#Hard filter variants: SNPs
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' VariantFiltration \
    -R $ref -V genotyped_recal1_snps.vcf \
    -O genotyped_recal1_snps_filtered.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "SOR > 4.0" --filter-name "SOR4" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-filter "MQRankSum<-12.5" --filter-name "MQRS-12.5" \
	-filter "ReadPosRankSum<-8.0" --filter-name "RPRS-8"
#Hard filter variants: INDELs
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' VariantFiltration \
    -R $ref -V genotyped_recal1_indels.vcf \
    -O genotyped_recal1_indels_filtered.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "SOR > 10.0" --filter-name "SOR10" 

#Remove filtered snps from vcf file
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SelectVariants \
    -R $ref -V genotyped_recal1_snps_filtered.vcf \
    -O genotyped_recal1_snps_onlyPASS.vcf \
    --exclude-filtered true
#Remove filtered indels from vcf file
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SelectVariants \
    -R $ref -V genotyped_recal1_indels_filtered.vcf \
    -O genotyped_recal1_indels_onlyPASS.vcf \
    --exclude-filtered true
#####END PBS
qsub ../../pbsscripts/5_GVCFgenotype2-filter2.pbs

#Check outputs (grep -v "^#" to exclude the headers that start with '#')
grep -v "^#" genotyped_recal1_snps.vcf | wc -l #186475
grep -v "^#" genotyped_recal1_snps_onlyPASS.vcf | wc -l #158599
#This does the same thing in gatk:
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk CountVariants -V genotyped_recal1_snps_onlyPASS.vcf #158599

#############################################################

### 6. 2nd Base recalibration and variant call on recalibrated bams 
    #ARRAY JOB #4 - 6_recal2_variantcall3_array.pbs, 160samples ~ 000min

#####BEGIN PBS
#!/bin/bash
#PBS -j oe
#PBS -m ae
#PBS -N 6_recal2_variantcall3_array
#PBS -l select=1:ncpus=5:mem=50gb
#PBS -l walltime=15:00:00
#PBS -M my@email
#PBS -J 0-176

shopt -s expand_aliases
#Requires: GATK4, R

cd /scratch/snp_calling/Clade1_refC99/6_out_recal2variantcall3
#Reference taxon
ref=/scratch/snp_calling/Clade1_refC99/1_dataprep/ref-loci.fasta
#bams after removing duplicates with picard
id=$PBS_ARRAY_INDEX
dedup_bams=($(ls /scratch/snp_calling/Clade1_refC99/1_dataprep/*all_dedup.bam))
dedup_bam=${dedup_bams[$id]}
#Create a variable with the dedup_bam. Use path variable to extract the last filename with cut -d
sample=$(echo $dedup_bam | cut -d/ -f7 | cut -d. -f1)
echo $sample

####Base recalibration 2
echo "Processing $dedup_bam"
#Get pre-recal table with uncalibrated bam but "refined" SNPs and INDELs
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx4g' BaseRecalibrator \
    -R $ref -I $dedup_bam -O ${sample}_pre-recal2.table \
    --known-sites ../5_out_GVCF2filter2/genotyped_recal1_snps_onlyPASS.vcf \
    --known-sites ../5_out_GVCF2filter2/genotyped_recal1_indels_onlyPASS.vcf
#Now apply 'true' base score correction on the crude samples in preparation for a 'true' convergence check
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx4g' ApplyBQSR \
    -R $ref -I $dedup_bam --bqsr-recal-file $sample'_pre-recal2.table' \
    -O ${sample}_recal2.bam
#Get Post-recal table on the newly recalibrated bam to compare in AnalyzeCovariates
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx4g' BaseRecalibrator \
    -R $ref -I ${sample}_recal2.bam -O ${sample}_post-recal2.table \
    --known-sites ../5_out_GVCF2filter2/genotyped_recal1_snps_onlyPASS.vcf \
    --known-sites ../5_out_GVCF2filter2/genotyped_recal1_indels_onlyPASS.vcf
#Analyze covariates to see whether it converged
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx4g' AnalyzeCovariates \
    -before ${sample}_pre-recal2.table \
    -after ${sample}_post-recal2.table \
    -plots ${sample}_AnalyzeCovariates2.pdf \
    -csv ${sample}_AnalyzeCovariates2.csv
#Also see if i made any difference with the first post-recal
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx4g' AnalyzeCovariates \
    -before ../4_out_recal1variantcall2/${sample}_post-recal1.table \
    -after ${sample}_post-recal2.table \
    -plots ${sample}_AnalyzeCovariates3.pdf \
    -csv ${sample}_AnalyzeCovariates3.csv

    #Rerun BQSR.R if it doesn't run with GATK
/home/R/BQSR.R ${sample}_AnalyzeCovariates2.csv ${sample}_pre-recal2.table ${sample}_ext_AnalyzeCovariates2.pdf
/home/R/BQSR.R ${sample}_AnalyzeCovariates3.csv ../4_out_recal1variantcall2/${sample}_post-recal1.table ${sample}_ext_AnalyzeCovariates3.pdf

#HaplotypeCaller in GATK for variant discovery on officially recalibrated bams
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx4g' HaplotypeCaller \
    -R $ref -I ${sample}_recal2.bam -O ${sample}_recal2.g.vcf \
    --contamination-fraction-to-filter 0.0002 \
    --emit-ref-confidence GVCF \
    --min-base-quality-score 20 \
    --phred-scaled-global-read-mismapping-rate 30 \
    --standard-min-confidence-threshold-for-calling 40.0
#####END PBS
qsub ../../pbsscripts/6_recal2_variantcall3_array.pbs

#At this stage, check the before/after plots for convergence with AnalyzeCovariates. If converged, we can do a "true" snp calling in the next step. If not, repeat the genotyping, filtering, ApplyBQSR and HaplotypeCaller all over again. (modify scripts 5_GVCFgenotype2-filter2.pbs and 6_recal2_variantcall3_array.pbs)

#############################################################

### 7. Combine GVCF, joint genotype, extract SNPs/INDELs from the call set and hard filter
    #7_GVCFgenotype3-filter3.pbs, 160samples ~ 40min

#####BEGIN PBS
#!/bin/bash
#PBS -j oe
#PBS -m ae
#PBS -N 7_GVCFgenotype3-filter3
#PBS -l select=1:ncpus=5:mem=30gb
#PBS -l walltime=5:00:00
#PBS -M my@email

cd /scratch/snp_calling/Clade1_refC99/7_out_GVCF3filter3
#Reference taxon
ref='/scratch/snp_calling/Clade1_refC99/1_dataprep/ref-loci.fasta'

#Get paths for all vcf files and combine into a cohort gvcf
for i in ../6_out_recal2variantcall3/*_recal2.g.vcf; do echo " -V $i "; done | tr -d "\n" > gvcf_recal2.list
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' CombineGVCFs \
    -R $ref $(cat gvcf_recal2.list) -O cohort_recal2.g.vcf

#Genotyping with GVCF in the cohort variant files produced by HaplotypeCaller gvcf
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' GenotypeGVCFs \
    -R $ref -V cohort_recal2.g.vcf -O genotyped_recal2.vcf

#Extract the SNPs from the call set
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SelectVariants \
    -R $ref -V genotyped_recal2.vcf \
    -O genotyped_recal2_snps.vcf \
    --select-type-to-include SNP
#Extract the indels from the call set
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SelectVariants \
    -R $ref -V genotyped_recal2.vcf \
    -O genotyped_recal2_indels.vcf \
    --select-type-to-include INDEL

#Hard filter variants: SNPs and also mask INDELs
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' VariantFiltration \
    -R $ref -V genotyped_recal2_snps.vcf \
    -O genotyped_recal2_snps_filtered.vcf \
    --mask genotyped_recal2_indels.vcf --mask-extension 5 --mask-name INDEL \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "SOR > 4.0" --filter-name "SOR4" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum<-12.5" --filter-name "MQRS-12.5" \
    -filter "ReadPosRankSum<-8.0" --filter-name "RPRS-8"

#Remove filtered snps from vcf file
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk \
    --java-options '-Xmx8G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' SelectVariants \
    -R $ref -V genotyped_recal2_snps_filtered.vcf \
    -O genotyped_recal2_snps_onlyPASS.vcf \
    --exclude-filtered true

#END OF SNP CALLING WITH GATK4. 
#The final output 'genotyped_recal2_snps_onlyPASS.vcf' can now be used for filtering, subsetting, then structure analysis, SNAPP etc.

#####END PBS
qsub ../../pbsscripts/7_GVCFgenotype3-filter3.pbs

#Check outputs (grep -v "^#" to exclude the headers that start with '#')
grep -v "^#" genotyped_recal2_snps.vcf | wc -l #182182
grep -v "^#" genotyped_recal2_snps_onlyPASS.vcf | wc -l #107784
singularity run /fast/tmp/containers/gatk4-4.3.0.0.sif gatk CountVariants -V genotyped_recal2_snps_onlyPASS.vcf #same thing as grep
grep -v "^#" genotyped_recal2_snps_onlyPASS.vcf | grep "INDEL" genotyped_recal2_snps_onlyPASS.vcf #should be none

#############################################################

### 8. Get stats and prepare vcf for downstream analysis

#### Compress, index, and get stats with bcftools
singularity run /fast/tmp/containers/tabix-1.13.sif bgzip \
    -c genotyped_recal2_snps_onlyPASS.vcf > ../8_out_fordownstream/genotyped_recal2_snps_onlyPASS.vcf.gz
singularity run /fast/tmp/containers/tabix-1.13.sif tabix \
    -p vcf ../8_out_fordownstream/genotyped_recal2_snps_onlyPASS.vcf.gz
cd ../8_out_fordownstream
singularity run /fast/tmp/containers/bcftools-1.16.sif bcftools view -H genotyped_recal2_snps_onlyPASS.vcf.gz | wc -l #same thing as grep
singularity run /fast/tmp/containers/bcftools-1.16.sif bcftools stats genotyped_recal2_snps_onlyPASS.vcf.gz > stats_snps_recal2.csv

#### Filter with bcftools (parameters following Erickson et al. 2020)
#75% filtering (-e:exclude)
singularity run /fast/tmp/containers/bcftools-1.16.sif bcftools view -e 'AC==0 || AC==AN || F_MISSING>=0.25' -m2 -M2 -O z -o filtered_snps_75.vcf.gz genotyped_recal2_snps_onlyPASS.vcf.gz
singularity run /fast/tmp/containers/bcftools-1.16.sif bcftools +prune -w 1000bp -n 1 -N 1st -o filtered_snps_75.vcf filtered_snps_75.vcf.gz
singularity run /fast/tmp/containers/bcftools-1.16.sif bcftools view -H filtered_snps_75.vcf | wc -l  #refC99:1948
#Equivalent to:
    #singularity run /fast/tmp/containers/vcftools-0.1.16.sif vcftools --vcf ../7_out_GVCF3filter3/genotyped_recal2_snps_onlyPASS.vcf --out filtered_snps_75 --thin 1000 --max-non-ref-af 0.99 --min-alleles 2 --max-alleles 2 --max-missing 0.75 --recode
        #- --thin 1000 to ensure that no two snps are within 1000 bp of one another (a proxy filter to get 1 snp per locus)
        #- --max-non-ref 0.99 to filter out any SNPs that are the same across all samples
        #- --min-alleles 2 & --max-alleles 2 to only include biallelic sites
        #- --max-missing 0.75 parameter to set a threshold for taxon completeness (e.g. 0.75 only includes SNPs for at least 75% of individuals have an allele thus less than 25% of individuals have missing data)
        #- --recode to create a new vcf instead of overwriting existing file
        #- More parameters can be found [here](http://vcftools.sourceforge.net/man_0112a.html)

#95% filtering (-e:exclude)
singularity run /fast/tmp/containers/bcftools-1.16.sif bcftools view -e 'AC==0 || AC==AN || F_MISSING>=0.05' -m2 -M2 -O z -o filtered_snps_95.vcf.gz genotyped_recal2_snps_onlyPASS.vcf.gz
singularity run /fast/tmp/containers/bcftools-1.16.sif bcftools +prune -w 1000bp -n 1 -N 1st -o filtered_snps_95.vcf filtered_snps_95.vcf.gz
singularity run /fast/tmp/containers/bcftools-1.16.sif bcftools view -H filtered_snps_95.vcf | wc -l  #1941

#Subset for subclades,then 95% filtering
##1A
singularity run /fast/tmp/containers/bcftools-1.16.sif bcftools view -S samples_1A.list -O z -o 1A_snps.vcf.gz genotyped_recal2_snps_onlyPASS.vcf.gz
singularity run /fast/tmp/containers/bcftools-1.16.sif bcftools view -e 'AC==0 || AC==AN || F_MISSING>=0.05' -m2 -M2 -O z -o filtered_1A_snps_95.vcf.gz 1A_snps.vcf.gz
singularity run /fast/tmp/containers/bcftools-1.16.sif bcftools +prune -w 1000bp -n 1 -N 1st -o filtered_1A_snps_95.vcf filtered_1A_snps_95.vcf.gz
singularity run /fast/tmp/containers/bcftools-1.16.sif bcftools view -H filtered_1A_snps_95.vcf | wc -l  #1924

#### Format with plink for structure/dapc
singularity run /fast/tmp/containers/plink-1.90b6.26.sif plink1.9 --vcf filtered_noOG_snps_75.vcf \
    --const-fid 0 \
    --allow-extra-chr \
    --recode 'structure' \
    --out structure_snps_vcf75
    #--const-fid: converts sample IDs to within-family IDs while setting all family IDs to a single value (default '0'). Basically it tells Plink to use the same name as those set in the vcf file
    #--allow-extra-chr: allows new contig name format e.g uce-...
#This spits out 3 files, one of which we need with a weird extension
mv structure_snps_vcf75.recode.strct_in structure_snps_vcf75.str
#Continue to edit the structure formatted file output to remove Cnumber's 'C' and the first row because STRUCTURE doesn't like it
sed '1,2d' structure_snps_vcf75.str > structure_snps_vcf75-edited.str #sed 's/C//; 1,2d'

#95
singularity run /fast/tmp/containers/plink-1.90b6.26.sif plink1.9 --vcf filtered_noOG_snps_95.vcf --const-fid 0 --allow-extra-chr --recode 'structure' --out structure_snps_vcf95
mv structure_snps_vcf95.recode.strct_in structure_snps_vcf95.str
sed '1,2d' structure_snps_vcf95.str > structure_snps_vcf95-edited.str


#### BEFORE USING THE VCF/STR OUTPUT, check mean depth of coverage per indiv  and amount of missing data
singularity run /fast/tmp/containers/vcftools-0.1.16.sif vcftools  --vcf refC99_1C_snps_95.vcf --depth --out stats/refC99_1C_snps_95
#proportion missing data per indiv
singularity run /fast/tmp/containers/vcftools-0.1.16.sif vcftools  --vcf refC99_1C_snps_95.vcf --missing-indv --out stats/refC99_1C_snps_95
#other stats see: https://speciationgenomics.github.io/filtering_vcfs/
