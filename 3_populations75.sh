#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=POPS75
#SBATCH --output=%x.o%j
#SBATCH --error=%x.o%j
#SBATCH --partition nocona
#SBATCH --nodes=1
#SBATCH --ntasks=12

module load gcc/10.1.0 stacks/2.55

DIR="/lustre/scratch/emilwrig/aoudad/final_analyses/populations/"
RUNNAME="populations_75"
POPMAP="/lustre/scratch/emilwrig/aoudad/final_analyses/info/popmap.txt"
SORTED="/lustre/scratch/emilwrig/aoudad/final_analyses/03_stacks"
mkdir $DIR
cd $DIR
mkdir $RUNNAME

populations -P $SORTED --popmap $POPMAP -O ./$RUNNAME  -R 0.75 -t 12 --vcf --structure --hwe 

echo "Done"

