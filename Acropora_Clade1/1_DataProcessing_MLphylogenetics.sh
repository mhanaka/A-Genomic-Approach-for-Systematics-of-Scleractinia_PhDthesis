#!/bin/bash

#Script for basic phylogenomics using target capture of ultraconserved elements (UCEs) and phylogenetic reconstruction with maximum likelihood
#Using Clade1 of genus Acropora data (Scleractinia: Acroporidae)
#Author: Hanaka Mera   mailto:hanaka.mera[@]my.jcu.edu.au
    #Tutorials followed (in this script):
    #Phyluce (phylogenomics): https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html
    #Phyluce (genome harvesting): https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-3.html

#Last edited: July 2023

###################################################################

## Summary
##0. Data
##1. Organise clean reads
### 1-A. Raw UCE reads
### 1-B. Publicly available genomes: NCBI
### 1-C. Other sequencing formats e.g. WGS
##2. Assemble
##3. Correct assembled contigs
##4. Match contigs to probe
##5. Extract UCE loci of interest
##6. Align UCE loci
##7. ML trees

###################################################################

## 0. Data

#1. Targeted enrichment and capture of UCEs for Acropora Clade1
    #Combining data from multiple sequencing runs = samples
#2. Publicly available genome that match PSH names
    #A.tenuis, A.yongei = 2 samples
    #We already know A.echinata is not echinata from Bridge et al 2023 so not included
#3. Whole genome sequencing data shared by collaborators
#4. Outgroup

#NOTE:
#'singularity run ... ' or 'module load conda ...' before most commands are high performance cluster / local machine specific. Modify as see fit

#####################################################################

## 1. Organise clean reads and downloaded data

mkdir -p /home/{probes,1_genomes,1_raw,2_clean,2_subsample,3_assembly,3_correction,4_match,5_align,6_tree}
#Download UCE probe set
#https://doi.org/10.1016/j.ympev.2020.106944 under Appendix

### 1-A. Raw UCE reads

#Raw sequence reads were first checked for read counts then cleaned of adapters. 
#Data from before 2023 were processed with Illumiprocessor wrapper program (Faircloth et al. 2012) for trimmomatic (Bolger et al. 2014; see below step1-A-1).
#Data from 2023 and onwards were processed with fastp v0.23.2 (Chen 2023; see below step1A-2).

#### 1-A-1. Trim adapters with Illumiprocessor
illumiprocessor --input /home/1_raw/ --output 2_clean --config illumi.conf --cores 5
    #illumi.conf requires adapter information, see https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html#clean-the-read-data
#Count reads manually
for i in 2_clean/*_R1_*.fastq.gz; do 
    echo $i; gunzip -c $i | wc -l | awk '{print $1/4}'; done
#OR get stats
printf 'sample,reads\n' > 2_clean/results_illumi.txt
for i in $(cat dir_clean.list); do
    sample=$(echo $i | sed 's/\(.*\)\//\1/; s/.*\///')
    echo $sample','"$(grep -s "Input Read Pairs" "$i"stats/*)" >> 2_clean/results_illumi.txt; done

#### 1-A-2. Preprocess with fastp
cd /home/1_raw
for dir in *; do
    cd $dir
    R1=$(ls *_1*)
    R2=$(ls *_2*)
    singularity run /sw/containers/fastp-0.23.2.sif fastp 
        -i $R1 -I $R2 \
        -o ../../2_clean/"$dir"/"$dir"-READ1.fastq.gz \
        -O ../../2_clean/"$dir"/"$dir"-READ2.fastq.gz \
        -j ../../2_clean/"$dir"/"$dir"-READ.json \
        -h ../../2_clean/"$dir"/"$dir"-READ.html \
        --thread 1 \
    cd ..; done
#Get stats
printf 'sample,before_totalreads,after_totalreads\n' > 2_clean/results_fastp.txt
for i in $(ls -d 2_clean/*/); do
    sample=$(echo $i | sed 's/\(.*\)\//\1/; s/.*\///')
    echo $sample','"$(grep -s -A1 '"before_filtering"' $i"$sample"-READ.json | \
        tail -n1)""$(grep -s -A1 '"after_filtering"' $i"$sample"-READ.json | \
        tail -n1)" >> 2_clean/results_fastp.txt; done

### 1-B. Publicly available genomes: NCBI

#Download, convert .fna/.fasta to .2bit
cd /home/1_genomes
#Get full list of currently available genome sequences
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
    #For RefSeq use: wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
#Then change the genus/species names to what you're after and check that it is the correct ones
grep -E 'Acropora.*' assembly_summary_genbank.txt | cut -f 1,8,12,15,16,19 > ftp_download.list
#Create links to wget for each dataset
grep -E 'Acropora.*' assembly_summary_genbank.txt | cut -f 20 | sed 'h; s/^.*\///; H; g; s/\n/\//' > ftp_folder.list
#For loop through to wget each .fna and unzip
while read line; do
    sample=$(echo $line | sed 's/.*\///')
    wget $line'_genomic.fna.gz'
    gunzip $sample*; done < ftp_folder.list
#Create new directory for each species and move. Note: sqlite (phyluce match contigs to probe step) doesn't like period in file names so had to change the last dot to underscore.
sed -i 's/ /_/g' download_list.txt 
while read id taxon chromo short; do
  id_=$(echo $id | sed 's/\./_/')
  sample=$(echo ${taxon}_${id_})
  mkdir $sample
  mv $id* ${sample}/; done < download_list.txt
#Delete intermediate files
rm ftp_*

#Others e.g.
wget http://aten.reefgenomics.org/download/aten_final_0.11.fasta.gz
gunzip -c aten_final_0.11.fasta.gz > Acropora_tenuis_ReefGenomics_v0_11.fna

#Once you have all the data, convert to .2bit and run genome harvesting steps with Phyluce:

#Convert .fna to .2bit
cd /home/1_genome/
for dir in Acropora_*; do
  cd ${dir}
  singularity run /fast/tmp/containers/ucsc-software-20221018.sif faToTwoBit *.fna ${dir}'.2bit'
  cd ..; done

#Extract .fasta that match UCE loci from genome and treat as corrected contigs (after step ##)
#Create sqlite database
ls -1d archives/Acropora_* | sed 's/archives\///; s/\// /' | tr -d '\n'  
    #Then copy paste under 'scaffoldlist' below, in 1 line, spece in between taxa
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_probe_run_multiple_lastzs_sqlite \
    --db genomes.sqlite \
    --output genome-lastz \
    --scaffoldlist Acropora_tenuis_GCA_014633955_1 Acropora_yongei_GCA_014634225_1 \
    --genome-base-path ./ \
    --probefile /home/probes/Cowman_etal_APPENDIX_C-hexa-v2-final-probes.fasta \
    --cores 1
#Make a configuration file
printf [scaffolds]'\n' > genomes.conf
ls -1d Acropora_* | sed 's/archives\///; h; s/$/:\/home\/1_genomes\//; G; s/\n//; G; s/\n/\//; s/$/\.2bit/' >> genomes.conf
#SO that it looks like this:
[scaffolds]
Acropora_tenuis_GCA_014633955_1:/home/1_genome/Acropora_tenuis_GCA_014633955_1/Acropora_tenuis_GCA_014633955_1.2bit
Acropora_yongei_GCA_014634225_1:/home/1_genome/Acropora_yongei_GCA_014634225_1/Acropora_yongei_GCA_014634225_1.2bit
#Then run:
phyluce_probe_slice_sequence_from_genomes \
    --lastz genome-lastz \
    --conf genomes.conf \
    --flank 500 \
    --name-pattern "Cowman_etal_APPENDIX_C-hexa-v2-final-probes.fasta_v_{}.lastz.clean" \
    --output genome-fasta
#Rename files because everything gets lowercased :(, then move to respective folders
cd genome-fasta
for fa in *; do
  sample=$(echo $fa | sed '1 s/a/A/; s/gca/GCA/')
  dir=$(echo $sample | sed 's/\.fasta//')
  mv $fa ../${dir}/${sample}; done
rm -r genome-*/
#Use the fasta output as corrected contigs

### 1-C. Other sequencing formats e.g. WGS

#Preprocess with fastp the same way as step1-A-2. If LARGE data e.g. whole genome, subsample reads to reduce data size.
cd /scratch/2_clean
for dir in /scratch/2_clean/*; do
    sample=$(echo $dir | sed 's/.*\///')
    echo $sample
    singularity run /sw/containers/seqtk-1.3.sif seqtk sample -s 999 "$dir"/"$sample"-READ1.fastq.gz 5000000 | gzip > /scratch/2_subsampled/"$sample"/"$sample"-READ1.fastq.gz
    singularity run /sw/containers/seqtk-1.3.sif seqtk sample -s 999 "$dir"/"$sample"-READ2.fastq.gz 5000000 | gzip > /scratch/2_subsampled/"$sample"/"$sample"-READ2.fastq.gz
done

###################################################################

## 2. Assemble

#Assemble reads with SPAdes
cd 3_assembly
nano spades.pbs
#####PBS START
#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N spades
#PBS -l select=1:ncpus=40:mem=250gb
#PBS -l walltime=30:00:00
#PBS -M my@email
#PBS -J 0-99

shopt -s expand_aliases

cd /home/2_clean
id=${PBS_ARRAY_INDEX}
clean_reads=($(ls -d *))
sample=${clean_reads[$id]}
r1=$(echo "/home/2_clean/$sample/${sample}-READ1.fastq.gz")
r2=$(echo "/home/2_clean/$sample/${sample}-READ2.fastq.gz")

mkdir /tmp/${sample}_spades_tmp
singularity run /sw/containers/spades-3.15.5.sif spades \
    -o /tmp/personal/${sample}_spades \
    -1 $r1 -2 $r2 --tmp-dir /tmp/${sample}_spades_tmp \
    --careful --threads 40 --cov-cutoff 2 --checkpoints last
rsync -av /tmp/personal/${sample}_spades 
    /home/3_assembly/${sample}/ \
    --exclude /tmp/personal/${sample}_spades_tmp
rm -r /tmp/${sample}_spades_tmp
rm -r /tmp/personal/${sample}_spades
#####PBS END
qsub spades.pbs

#Symlink to resulting contigs.fasta file and run QC
mkdir contigs
for i in /home/3_assembly/*/*_spades/contigs.fasta; do
    sample=$(echo $i | sed 's/\..*/\.contigs\.fasta/; s/.*\///')
    ln -s $i /home/3_assembly/contigs/"$sample".contigs.fasta
done
printf 'sample,contigs,total bp,mean length,95 CI length,min length,max length,median length,contigs >1kb''\n' > 3_assembly/results_assembly.txt
for i in /home/3_assembly/contigs/*; do
    singularity run /sw/containers/phyluce-1.7.1.sif phyluce_assembly_get_fasta_lengths --input $i --csv >> 3_assembly/results_assembly.txt
done

###################################################################

## 3. Correct assembled contigs

cd 3_correction
### Create config file for mapping
printf reads\:'\n' > mapping.conf
#Below should be in a format as:
    #[space]'samplename': /path/to/clean/reads/ (where READ1 and READ2 are)
for file in /home/3_assembly/contigs/*contigs.fasta; do
    i=$(echo $file | sed 's/\.contigs\.fasta//; s/.*\///')
    printf '  '"$i"\:' ''/home/2_clean/'"$i"'/''\n' >> mapping.conf; done
printf '\n'contigs:'\n' >> mapping.conf
#Below should be in a format as:
    #[space]'samplename': /path/to/assembled/contigs.fasta
for file in /home/3_assembly/contigs/*contigs.fasta; do
    i=$(echo $file | sed 's/\.contigs\.fasta//; s/.*\///')
    printf '  '$i\:' '$file'\n' >> mapping.conf; done

### Run mapping phyluce_workflow
nano mapping.pbs
#####BEGIN PBS
#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N wf_map
#PBS -l select=1:ncpus=5:mem=100gb
#PBS -l walltime=30:00:00
#PBS -M my@email

#Set JAVA memory
export SINGULARITYENV__JAVA_OPTS="-Xms800m -Xmx80g"

cd /home/3_correction
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_workflow \
    --config mapping.conf --output mapping_contigs \
    --workflow mapping --cores 5
#####END PBS
qsub mapping.pbs

#Copy the coverage stats
cp mapping_contigs/coverage/all-taxon.summary.csv .

### Create config file for correction
printf bams\:'\n' > correction.conf
#Change directory to the output of previous step (mapping/mapped_reads)
cd mapping_contigs/mapped_reads
for file in *.fxm.sorted.md.bam; do
    Cnum=$(echo $file | sed 's/\.fxm\.sorted\.md\.bam//')
    printf '  '$Cnum\:' '$PWD'/'$file'\n' >> /home/3_correction/correction.conf; done
printf '\n'contigs:'\n' >> /home/3_correction/correction.conf
cd /home
for file in /home/3_assembly/contigs/*contigs.fasta; do
    i=$(echo $file | sed 's/\.contigs\.fasta//; s/.*\///')
    printf '  '$i\:' '$file'\n' >> 3_correction/correction.conf; done

### Run correction phyluce_workflow
cd 3_correction
nano correction.pbs
#####BEGIN PBS
#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N wf_corr
#PBS -l select=1:ncpus=5:mem=100gb
#PBS -l walltime=10:00:00
#PBS -M my@email

#Set JAVA memory
export SINGULARITYENV__JAVA_OPTS="-Xms800m -Xmx80g"

cd /home/3_correction
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_workflow \
    --config correction.conf --output correction_contigs \
    --workflow correction --cores 4
#####END PBS
qsub correction.pbs

#Symlink to corrected contigs (and genome harvested ones)
mkdir correctedcontigs
cd /home
for file in 3_correction/correction_contigs/consensus/*; do 
    sample=$(echo $file | sed 's/\..*/\.contigs\.fasta/; s/.*\///')
    ln -s "$PWD"/"$file" "$PWD"/3_correction/correctedcontigs/"$sample"; done
for file in /home/1_genome/Acropora_*; do
    sample=$(echo $file | sed 's/.*\///')
    ln -s "$file"/"$sample".fasta "$PWD"/3_correction/correctedcontigs/"$sample".fasta; done

#Run QC
cd /home/3_correction
printf 'sample,contigs,total bp,mean length,95 CI length,min length,max length,median length,contigs >1kb''\n' > results_correctedcontigsQC.txt
for i in correctedcontigs/*.fasta; do
    singularity run /sw/containers/phyluce-1.7.1.sif phyluce_assembly_get_fasta_lengths --input $i --csv >> results_correctedcontigsQC.txt; done
    
###################################################################

## 4. Match contigs to probe

cd 4_match
mkdir log
nano match.pbs
#####BEGIN PBS
#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N match
#PBS -l select=1:ncpus=2:mem=100gb
#PBS -l walltime=30:00:00
#PBS -M my@email

cd /home
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_assembly_match_contigs_to_probes \
    --contigs 3_correction/correctedcontigs \
    --probes /home/probes/Cowman_etal_APPENDIX_C-hexa-v2-final-probes.fasta \
    --min-identity 70 \
    --min-coverage 70 \
    --output 4_match/search-results
#####END PBS

#Get match stats
printf 'sample,unique_contigs,percent,contigs,dupe match,removed_multicontig,removed_multiUCE''\n' > results_searchQC.txt
cut -d '-' -f 6 log/phyluce_assembly_match_contigs_to_probes.log | \
    sed 's/ dupe probe matches//; s/ UCE loci removed for matching multiple contigs//; s/ contigs removed for matching multiple UCE loci//; s/\: /,/' | sed -e 's/ (/,/; s/) uniques of /,/; s/ contigs, /,/; s/ //g' | \
    head -n -4 | sed 1,17d >> results_searchQC.txt

###################################################################

## 5. Extract UCE loci of interest

#Create matrix configuration file
printf '[all]''\n' > taxon-set.conf
for i in search-results/*.lastz; do 
    sample=$(echo $i | sed 's/.*\///; s/\..*//'); 
    printf $sample'\n' >> taxon-set.conf; done
mkdir -p taxon-all/log
cd taxon-all
#Get match counts
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_assembly_get_match_counts \
    --locus-db ../search-results/probe.matches.sqlite \
    --taxon-list-config ../taxon-set.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output taxa-incomplete.conf \
    --log-path log
#Get monolithic fasta files
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../../3_correction/correctedcontigs \
    --locus-db ../search-results/probe.matches.sqlite \
    --match-count-output taxa-incomplete.conf \
    --output taxa-incomplete.fasta \
    --incomplete-matrix taxa-incomplete.incomplete \
    --log-path log
#Explode fastas and get summary stats
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_assembly_explode_get_fastas_file \
    --input taxa-incomplete.fasta \
    --output exploded-fastas \
    --by-taxon
printf 'sample,contigs,total bp,mean_length,95_CI_length,min_length,max_length,median_length,contigs_>1kb''\n' > results_explodedfastas.txt
for i in exploded-fastas/*.fasta; do
    singularity run /sw/containers/phyluce-1.7.1.sif phyluce_assembly_get_fasta_lengths --input $i --csv >> results_explodedfastas.txt;
done

######################################################################

## 6. Align UCE loci (edge trim)

cd 5_align
mkdir log

#####BEGIN PBS
#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N align
#PBS -l select=1:ncpus=5:mem=20gb
#PBS -l walltime=80:00:00
#PBS -M my@email

cd /home/5_align

#Align
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_seqcap_align \
    --input ../4_match/taxon-all/taxa-incomplete.fasta \
    --output mafft-nexus-edge-trimmed \
    --taxa 180 \
    --aligner mafft \
    --cores 5 \
    --incomplete-matrix \
    --log-path log
#Summary stats
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_align_summary_data \
    --alignments mafft-nexus-edge-trimmed \
    --cores 5 \
    --log-path log \
    --show-taxon-counts 
#Remove locus names
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_remove_locus_name_from_files \
    --alignments mafft-nexus-edge-trimmed \
    --output mafft-nexus-edge-trimmed-clean \
    --cores 5 \
    --log-path log
#Get % complete matrix (change percent value and output name)
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-edge-trimmed-clean \
    --taxa 180 \
    --percent 0.50 \
    --output mafft-nexus-edge-trimmed-clean-50p \
    --cores 5 \
    --log-path log
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-edge-trimmed-clean \
    --taxa 180 \
    --percent 0.70 \
    --output mafft-nexus-edge-trimmed-clean-70p \
    --cores 5 \
    --log-path log
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_only_loci_with_min_taxa \
    --alignments mafft-nexus-edge-trimmed-clean \
    --taxa 180 \
    --percent 0.80 \
    --output mafft-nexus-edge-trimmed-clean-80p \
    --cores 5 \
    --log-path log
#Get parsimonious informative site stats
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_informative_sites \
    --alignments mafft-nexus-edge-trimmed-clean-50p \
    --output infosites-e50.csv \
    --cores 5 \
    --log-path log
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_informative_sites \
    --alignments mafft-nexus-edge-trimmed-clean-70p \
    --output infosites-e70.csv \
    --cores 5 \
    --log-path log
singularity run /sw/containers/phyluce-1.7.1.sif phyluce_align_get_informative_sites \
    --alignments mafft-nexus-edge-trimmed-clean-80p \
    --output infosites-e80.csv \
    --cores 5 \
    --log-path log
#####END PBS

#alignment stats: length ± sampleSD, PIsite% 
    #see https://stackoverflow.com/questions/18786073/compute-average-and-standard-deviation-with-awk
cut -d, -f2,3 infosites-e50.csv | sed 1d | awk -F ',' '{ sum_length += $1; sumsq_length += ($1)^2; sum_pis += ($2/$1); sumsq_pis += ($2/$1)^2 } END { if (NR > 0) print "mean length: " sum_length / NR, "± " sqrt((sumsq_length-sum_length^2/NR)/(NR-1)), "SD bp; mean percent PIS: " sum_pis / NR*100, "± " sqrt((sumsq_pis-sum_pis^2/NR)/(NR-1))*100, "SD %" }'

######################################################################

## 7. ML trees
cd 6_tree
nano iqtree.pbs
#####BEGIN PBS
#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N ML_e50_cl1
#PBS -l select=1:ncpus=5:mem=100gb
#PBS -l walltime=200:00:00
#PBS -M my@email

alignment='/home/5_align/mafft-nexus-edge-trimmed-clean-50p'
run='e50.1832loci'

cd /home/6_tree

#Infer a concatenation-based species tree with 1000 UFbootstrap and SH-aLRT test, and an edge-linked partition model
singularity run /sw/containers/iqtree-2.2.2.2.sif iqtree2 \
    -p $alignment --prefix concat.$run -B 1000 -alrt 1000 \
    -T AUTO --threads-max 5 -m MFP+MERGE -rcluster 10 -alninfo --cptime 300
#Compute sCF using likelihood with v2.2.2
singularity run /sw/containers/iqtree-2.2.2.2.sif iqtree2 \
    -te concat.${run}.treefile -p $alignment --scfl 100 \
    --prefix siteconcord.$run -T AUTO --threads-max 5 --cf-verbose --cf-quartet

mkdir $run
mv {concat,loci,geneconcord,siteconcord}."$run"* "$run"/
#####END PBS
qsub iqtree.pbs

#Check PIsites
sed -n '/SEQUENCE ALIGNMENT/,/Column meanings/p' siteconcord.*.iqtree | head -n-2 | tail -n+7 | \
    awk '{i++; pis=$6/$4; total+=pis; mean=total/i; dev=(pis-mean)^2; var+=dev} END {print "[mean]",mean*100,"%\n[stdev]",sqrt(var/i)*100,"%"}'
