#!/bin/bash

#Script for ML phylogenetics with mitochondria data
#Using target capture of ultraconserved elements (UCEs) off-target reads
#Author: Hanaka Mera   mailto:hanaka.mera[@]my.jcu.edu.au
    #Tutorials followed (in this script):
    #jasonleebrown repo: https://github.com/jasonleebrown/UCE_phyluce_pipeline/blob/master/README.md#mitogenome-pipeline

#Last edited: July 2025

##################################################################

## 1. Prepare input files

### Decide on reference and download previous studies' data

#Downloaded one Montipora capitata from https://www.ncbi.nlm.nih.gov/nuccore/OX421459.1
OX421459.1_Montipora_capitata.fasta

#Published dataset: van Oppen et al. 2004
    #Genbank acc. AY313547-AY313605
    #Downlaod using NCBI's edirect program suit
/Users/programs/edirect/epost -db nucleotide -input vanOppen2004.list | /Users/programs/edirect/efetch -format fasta > vanOppen2004_mtCR.fasta
sed -i -e 's/ mitochondrial control region, partial sequence/\|mtDNA/; s/ /_/g' vanOppen2004_mtCR.fasta

#HPC
mkdir mito/{1_ref,2_clean,3_bwamem,4_angsd,5_align,6_tree}
nano samples.list
    #All UCE sample numbers in one column

#then rsync all dataset from local
rsync -LavP local/*.fasta mito/1_ref/

for i in $(cat samples.list); do
    fullpath=$(grep $i/ cleaned_reads_fullpath_to_directory.list)
    mkdir 2_clean/$i
    ln -s "$fullpath"*READ1* $PWD/2_clean/"$i"/"$i"-READ1.fastq.gz
    ln -s "$fullpath"*READ2* $PWD/2_clean/"$i"/"$i"-READ2.fastq.gz
done

### Index reference
cd 1_ref
singularity run /fast/tmp/containers/bwa-0.7.17.sif bwa index OX421459.1_Montipora_capitata.fasta

######################################################################

### 2. Run bwa-mem samtools

#Modified from bams_loopBWA-MEM-mtDNA.sh in jasonleebrown github repo

nano bwa.pbs
#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N bwa
#PBS -l select=1:ncpus=2:mem=10gb
#PBS -l walltime=1:00:00
#PBS -M hanaka.mera@my.jcu.edu.au
#PBS -J 0-99

cd /home/mito
id=$PBS_ARRAY_INDEX
clean_reads=($(ls 2_clean))
sample=$(echo ${clean_reads[$id]})

cd 2_clean/$sample
echo $sample
singularity run /fast/tmp/containers/bwa-0.7.17.sif bwa mem \
    ../../1_ref/OX421459.1_Montipora_capitata.fasta \
    *READ1.fastq.gz *READ2.fastq.gz -B 2 > $sample.sam
singularity run /fast/tmp/containers/samtools-1.20.sif samtools view \
    -uS $sample.sam | singularity run /fast/tmp/containers/samtools-1.20.sif samtools sort -o $sample.sorted.bam
rm $sample.sam
mv $sample.sorted.bam ../../3_bwamem/$sample.bam
cd ../..

####PBS END
qsub bwa.pbs

#qsub as many times needed for the number of samples, just change PBS -J 0-99 to 100-199 etc.

######################################################################

### 3. Run angsd

#Modified from angsd_Dofasta4_iupac0.2_minDepth_2_mtDNA.sh

for i in 3_bwamem/*; do 
    sample=$(echo $i | sed 's/.*\///; s/\.bam//')
    echo $sample
    singularity run /fast/tmp/containers/angsd-0.940.sif angsd \
    -i $i -doFasta 2 -doCounts 1 \
    -minQ 10 -minMapQ 20 -remove_bads 1 -uniqueOnly 1 -setMinDepth 1 \
    -out 4_angsd/$sample

    gunzip 4_angsd/$sample.fa.gz
    sed -i "s/>.*/>$sample |mtDNA/" 4_angsd/$sample.fa
    mv 4_angsd/$sample.fa 4_angsd/$sample.fasta
    rm 4_angsd/$sample.arg
done
    #-doFasta 2: Generate a fasta for a BAM file; 2=use the most common (non N) base (needs -doCounts 1) 
    #-doCounts 1
    #-minQ 10: min base quality 10
    #-minMapQ 20: min mapping quality 20
    #-remove_bads 1: remove reads with flag above 255
    #-uniqueOnly 1: remove reads with multiple hits
    #-setMinDepth 1: discard site if total sequencing depth is below 1

# Get stats
    #angsd adds "N" where there is no data. So:
#This prints total characters per fasta file (including Ns)
awk '/>.*/{if (x)print x;print;x="";next}{x=(!x)?$0:x$0;}END{print x;}' C100.fasta | tail -n+2 | awk '{print length}'
#This is without the Ns
awk '/>.*/{if (x)print x;print;x="";next}{x=(!x)?$0:x$0;}END{print x;}' C100.fasta | tail -n+2 | sed 's/[N]//g' | awk '{print length}'

printf 'sample,total_bp,total_bp_without_missing''\n' > results_angsd.txt
for i in *.fasta; do
    sample=$(echo $i | sed 's/\.fasta//')
    total=$(awk '/>.*/{if (x)print x;print;x="";next}{x=(!x)?$0:x$0;}END{print x;}' $i | tail -n+2 | awk '{print length}')
    ATGC=$(awk '/>.*/{if (x)print x;print;x="";next}{x=(!x)?$0:x$0;}END{print x;}' $i | tail -n+2 | sed 's/[N]//g' | awk '{print length}')
    printf $sample','$total','$ATGC'\n' >> results_angsd.txt; done

######################################################################

### 4. Align

cd 5_align

#concatenate all fastas
cat ../4_angsd/*.fasta > concatenated.fasta
cat concatenated.fasta vanOppen2004_mtCR.fasta > UCE_vanOppen2004_mtCR_allconcat.fasta

nano align.pbs
#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N align
#PBS -l select=1:ncpus=5:mem=50gb
#PBS -l walltime=30:00:00
#PBS -M hanaka.mera@my.jcu.edu.au

cd /home/mito/5_align
singularity run /fast/tmp/containers/mafft-7.505.sif mafft \
    --maxiterate 1000 --thread 4 \
    UCE_vanOppen2004_mtCR_allconcat.fasta > aligned.nex

singularity run /fast/tmp/containers/gblocks-0.91b.sif Gblocks \
    aligned.nex -t=d -b2=267 -b3=10 -b4=10 -b5=a -p=n
mv aligned.nex-gb aligned_gblocks.nex
#####PBS END
    #Original alignment: 141444 positions
    #Gblocks alignment:  16762 positions (11 %) in 299 selected block(s)

#local
rsync -LavP HPC:/home/mito/5_align/aligned_gblocks.nex .
cp aligned_gblocks.nex aligned_gblocks_trimmed.nex
#open in aliview, trim 6679~end, 1~6174
rsync -LavP aligned_gblocks_trimmed.nex HPC:/home/mito/5_align/.

#Remove samples with only ambiguous bases or gaps
    #This happens because I manually removed most of the data to match the regions from previous studies
sed '/^[^>]/ s/[atgc]//g' aligned_gblocks_trimmed.nex | awk '$0 ~ />.*/ {print $0} length>252 {print $0}' | grep -B1 "nnn\|---" --no-group-separator | grep ">" > tmp.remove
grep -f tmp.remove aligned_gblocks_trimmed.nex -A1 --no-group-separator > tmp2.remove
grep -vf tmp2.remove aligned_gblocks_trimmed.nex > aligned_gblocks_trimmed_cleaned2.nex
rm tmp*

# Get alignment stats
mkdir stats
cp aligned_gblocks_trimmed_cleaned2.nex stats/aligned_gblocks_trimmed_cleaned2.fasta
cd stats
sed -i -e '/>C.*/{h; s/ |mtDNA/_/; G; s/\n//; s/>//2}' aligned_gblocks_trimmed_cleaned.fasta 
sed -i -e '/>K.*/{h; s/ |mtDNA/_/; G; s/\n//; s/>//2}' aligned_gblocks_trimmed_cleaned.fasta 
singularity run /fast/tmp/containers/phyluce-1.7.3.sif phyluce_assembly_explode_get_fastas_file \
    --input aligned_gblocks_trimmed_cleaned.fasta \
    --output exploded-fastas \
    --by-taxon
printf 'sample,total_bp,total_bp_without_missing''\n' > results_alignment.txt
for i in exploded-fastas/*.fasta; do
    sample=$(echo $i | sed 's/.*\///; s/\.fasta//')
    total=$(awk '/>.*/{if (x)print x;print;x="";next}{x=(!x)?$0:x$0;}END{print x;}' $i | tail -n+2 | awk '{print length}')
    atgc=$(awk '/>.*/{if (x)print x;print;x="";next}{x=(!x)?$0:x$0;}END{print x;}' $i | tail -n+2 | sed 's/[n]//g' | awk '{print length}')
    printf $sample','$total','$atgc'\n' >> results_alignment.txt; done


######################################################################

### 5. Run iqtree
cd 6_tree

nano iqtree_single.pbs
#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N ML_mito
#PBS -l select=1:ncpus=5:mem=10gb
#PBS -l walltime=5:00:00
#PBS -M hanaka.mera@my.jcu.edu.au

alignment='/home/mito/5_align/aligned_gblocks_trimmed_cleaned.nex'

cd /home/mito/6_tree

#Infer a concatenation-based species tree with 1000 UFbootstrap and SH-aLRT test, and an edge-linked partition model
singularity run /fast/tmp/containers/iqtree-2.2.2.2.sif iqtree2 \
    -s $alignment --prefix concat -B 1000 -alrt 1000 \
    -T AUTO --threads-max 5 -m MFP+MERGE -rcluster 10 --cptime 300
####PBS END
qsub iqtree_single.pbs
