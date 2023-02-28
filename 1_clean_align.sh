#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=rad_align
#SBATCH --partition nocona
#SBATCH --nodes=1 
#SBATCH --ntasks=8 
#SBATCH --time=48:00:00
#SBATCH --array=1-2

# set up variables for this project
workdir=/lustre/scratch/emilwrig/aoudad/final_analyses
reference_genome=/lustre/scratch/emilwrig/aoudad/genome2_domesticgoat/GCF_001704415.2_ARS1.2_genomic.fna
array_input=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/info/list.txt | tail -n1 )

module load gcc/10.1.0 stacks/2.55

# filter the data with the STACKS process radtags script
process_radtags -1 ${workdir}/00_fastq/${array_input}_R1_001.fastq.gz -2 ${workdir}/00_fastq/${array_input}_R2_001.fastq.gz -o ${workdir}/01_cleaned -e ecoRI -c -q 

# trim the first 5 bp off (cut site)
trimmomatic SE -threads 1 -phred33 ${workdir}/01_cleaned/${array_input}_R1_001.1.fq.gz ${workdir}/01_cleaned/${array_input}_R1.fastq.gz HEADCROP:5
trimmomatic SE -threads 1 -phred33 ${workdir}/01_cleaned/${array_input}_R2_001.2.fq.gz ${workdir}/01_cleaned/${array_input}_R2.fastq.gz HEADCROP:5

# load new gcc and bwa because bwa requires a different gcc version?
module load gcc/9.2.0 bwa/0.7.17

# align to reference using bwa
bwa mem -t 8 ${reference_genome} ${workdir}/01_cleaned/${array_input}_R1.fastq.gz ${workdir}/01_cleaned/${array_input}_R2.fastq.gz > ${workdir}/02_bam/${array_input}.sam

# switch modules again
module load gcc/10.1.0 samtools

# convert sam to bam
samtools view -b -S -o ${workdir}/02_bam/${array_input}.bam ${workdir}/02_bam/${array_input}.sam

# sort the bam file
samtools sort ${workdir}/02_bam/${array_input}.bam > ${workdir}/02_bam/${array_input}_sorted.bam

# remove the unneeded sam and bam and fastq files
rm ${workdir}/02_bam/${array_input}.sam
rm ${workdir}/02_bam/${array_input}.bam
rm ${workdir}/01_cleaned/${array_input}_R1_001.1.fq.gz
rm ${workdir}/01_cleaned/${array_input}_R2_001.2.fq.gz
rm ${workdir}/01_cleaned/${array_input}_R1_001.rem.1.fq.gz
rm ${workdir}/01_cleaned/${array_input}_R2_001.rem.2.fq.gz

