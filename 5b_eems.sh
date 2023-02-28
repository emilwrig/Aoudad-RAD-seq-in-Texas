#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=bam
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G

source activate eems

module load gnu/5.4.0

/home/jmanthey/eems/runeems_snps/src/runeems_snps --params aoudad_eems.params --seed 123
