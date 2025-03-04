#!/bin/bash

#Script continuing from 1_DataProcessing_MLphylogenetics.sh and 2_SNPcalling.sh, this one is for running SPEEDEMON analysis with BEAST and SNAPPER
#Using Clade1 of genus Acropora data (Scleractinia: Acroporidae)
#Author: Hanaka Mera   mailto:hanaka.mera[@]my.jcu.edu.au
    #Tutorials and papers followed (in this script):
    #Douglas, J., Bouckaert, R., 2022. Quantitatively defining species boundaries with more efficiency and more biological realism. Commun Biol 5, 755. https://doi.org/10.1038/s42003-022-03723-z
    #https://github.com/rbouckaert/speedemon/tree/master/tutorial
#Last edited: July 2023

###################################################################

## SNAPPER + SPEEDEMON

### 1. Format using local PGDSpider for SNAPP (vcf to nexus)

#Use the subsetted vcf (relevant subclade) before it got formatted for anything
/scratch/snp_calling/Clade1_refC99/*_filtered_snps_95.vcf

#local
#Download PGDSpider2 and unzip: http://www.cmpg.unibe.ch/software/PGDSpider/
java -Xmx1024m -Xms512m -jar /Users/jc313394/bioinformatics/programs/PGDSpider_2.1.1.5/PGDSpider2.jar 

#First click on 'Create/Edit SPID file' (or open vcf-to-nexus.spid)
 #VCF-ParserQuestions
    #VCF to NEXUS
    #diplod, No, leave optionals blank, non-polymorphic No, population definitions No
 #NEXUS-WriterQuestions
    #Convert to binary format YES
 #Save and Apply
#Input file > format: vcf; select refC99_1*_snps_75.vcf
#Output file > format: nexus; select same vcf file, change the extension to .nex
#CONVERT

#Make a population file to add to the nexus output (it'll make the setting in BEAUti slightly easier)
grep -w -f samples_1A.list ../metadata_Cl1.csv | cut -d, -f1,5 | sed 's/,/\t/' > pop_1A.txt
#Rename taxon to 'PSH_Cnum'
while read Cnum taxon; do sed -i -e "s/$Cnum/${taxon}_$Cnum/" refC99_1A_snps_95.nex; done < pop_1A.txt
rm *.nex-e

#############################################################

### 2. Create .xml file in BEAUti

#Change taxon name to 'samplename'_sp so that each individual becomes its own species
#Change the prior to Yule Skyline Collapse
#Change epsilon to 1.0e-4 (or something else relevant to the taxa), save, then change to different number and save as a separate .xml. This allows us to test different epsilons
#Change mcmc chainLength="1000000" preBurnin="5000" storeEvery="1000"
#Under Operators tab (if not visible, change preference in Beauti), change GammaMover to 10.0, RateMixer to 10.0

#############################################################

### 3. Run BEAST

cd /scratch/SPEEDEMON
#Don't forget to add/update the BEAST packages
singularity run /fast/tmp/containers/beast2-2.7.4.sif packagemanager -list
singularity run /fast/tmp/containers/beast2-2.7.4.sif packagemanager -add speedemon

nano speedemon.pbs
#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N speed_1A_run1
#PBS -l select=1:ncpus=5:mem=30gb
#PBS -l walltime=200:00:00
#PBS -M my@email

xml='refC99_1A_8.5e-3'
clade='1A'
run='run1'

cd /scratch/SPEEDEMON/"$clade"/"$run"
singularity run /fast/tmp/containers/beast2-2.7.4.sif beast -threads 5 ../"$xml".xml #-resume
#### END PBS
qsub speedemon.pbs

#############################################################

### 4. Post-processing

#local
#If ran over multiple runs for ESS, check individual logs in Tracer, then combine logs and trees
for run in run*; do echo -log $run/*e-2*.trees  | tr -s '\n' ' '; done
/Applications/BEAST\ 2.7.4/bin/logcombiner -b 0 -log run1/refC99_1A_2.0e-2.trees -log run2/refC99_1A_2.0e-2.trees -log run3/refC99_1A_2.0e-2.trees -log run4/refC99_1A_2.0e-2.trees -log run5/refC99_1A_2.0e-2.trees -log run6/refC99_1A_2.0e-2.trees -o refC99_1A_2.0e-2_run1-6.trees

#Run ClusterTreeSetAnalyser to get SPEEDEMON results
/Applications/BEAST\ 2.7.4/bin/applauncher ClusterTreeSetAnalyser -trees refC99_1A_2.0e-2_run1-6.trees -out refC99_1A_2.0e-2_species.txt -epsilon 2.0e-2 -burnin 10
#Run TreeAnnotator to get .tre file
/Applications/BEAST\ 2.7.4/bin/treeannotator -burnin 10 refC99_1A_2.0e-2_run1-6.trees refC99_1A_2.0e-2_run1-6.tre
