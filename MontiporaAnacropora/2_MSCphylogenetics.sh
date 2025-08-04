#!/bin/bash

#Script for phylogenetics with multi-species coalescent using target capture of ultraconserved elements (UCEs)
#Using genus Montipora and Anacropora data (Scleractinia: Acroporidae)
#Author: Hanaka Mera   mailto:hanaka.mera[@]my.jcu.edu.au
    #Tutorials followed (in this script):
    #ASTRAL: https://github.com/smirarab/ASTRAL/tree/master
    #ASTRAL tutorial: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md
    #Sayyari and Mirarab 2016 https://doi.org/10.1093/molbev/msw079

#Last edited: July 2025

##################################################################

### 1. Run iqtree gene trees for each alignment

#Use output after phyluce_align_get_only_loci_with_min_taxa, steps after phyluce_align_seqcap_align

nano iqtree_for_astral_0-99.pbs
#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N ML_loci
#PBS -l select=1:ncpus=5:mem=50gb
#PBS -l walltime=50:00:00
#PBS -M hanaka.mera@my.jcu.edu.au
#PBS -J 0-99

shopt -s expand_aliases
id=$PBS_ARRAY_INDEX

run='e70.newMon.remEx.1000loci'
alignments='/home/UCE/newMontipora/5_align/taxon-set/mafft-nexus-edge-trimmed-clean-70p'
alignments_list=($(ls $alignments))
aln=${alignments_list[$id]}
aln_uce=$(echo $aln | sed 's/.*\///; s/\..*//')

mkdir -p /tmp/"$run"_$id
outdir="/tmp/${run}_${id}/$aln_uce"

#Infer gene/locus trees
    #note that -m MFP is not listed because this is default behaviour since iqtree v1.5.4
singularity run /fast/tmp/containers/iqtree-2.2.2.2.sif iqtree2 \
    -s $alignments/"$aln" --prefix $outdir -B 1000 -T AUTO --threads-max 5 --cptime 300

rsync -Lav /tmp/"$run"_"$id"/ /scratch/ASTRAL/1_iqtree_"$run"/"$aln_uce"
rm -r /tmp/"$run"_$id
#####PBS END #####

qsub iqtree_for_astral_0-99.pbs

**You will have to run the array jobs as many times as the number of alignments**
for i in {1..9}; do 
    pbs=$(echo iqtree_for_astral_0-99.pbs | sed "s/0-99/${i}00-${i}99/")
    sed "s/0-99/${i}00-${i}99/" iqtree_for_astral_0-99.pbs > $pbs
    qsub $pbs; done
rm iqtree_for_astral_*00*

##################################################################

### 8-3. Collapse low bootstraps (ASTRAL-III tutorial did <10)

#Chose BS<33 following Sayyari and Mirarab 2016 Molecular Biology and Evolution
for i in 1_iqtree_e70.newMon.remEx.1000loci/*/*.treefile; do
    output=$(echo $i | sed 's/.*\///; s/\.treefile//')
    echo $output
    singularity run /fast/tmp/containers/newick_utils-1.6.sif nw_ed $i 'i & b<=33' o > 2_collapsed_e70.newMon.remEx.1000loci/$output.BS33collapsed.tre; done
#Concatenate trees into one file
cat 2_collapsed_e70.newMon.remEx.1000loci/*.tre > 3_astral_e70.newMon.remEx.1000loci/e70.newMon.remEx.1000loci_BS33collapsed.tre

##################################################################

### 8-4. Run ASTRAL

#Fast enough to run on login node if fewer loci, but better to submit pbs job

nano astral.pbs
#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N sp_e50.newMon.remEx
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -l walltime=50:00:00
#PBS -M hanaka.mera@my.jcu.edu.au

cd /scratch/ASTRAL/3_astral_e70.newMon.remEx.1000loci

singularity run /fast/tmp/containers/astral-5.7.8.sif java -jar /opt/Astral/astral.5.7.8.jar \
    -i e70.newMon.remEx.1000loci_BS33collapsed.tre -t 2 \
    -o e70.newMon.remEx.1000loci_BS33collapsed_astral.tre \
    2> e70.newMon.remEx.1000loci_BS33collapsed_astral.log

####PBS END

mv 3_* /home/UCE/newMontipora/6_tree/ASTRAL/

#Run the add-bl.py script from https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md to add arbitrary terminal branch lengths just for plotting sake

#local
#rsync the astral tree output
./add-bl.py postfix e70.newMon.remEx.1000loci_BS33collapsed_astral.tre
