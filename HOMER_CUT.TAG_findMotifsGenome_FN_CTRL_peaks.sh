#!/bin/bash -l
#SBATCH -J HOMER_CUT&TAG_findMotifsGenome
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=deborah.gerard@uni.lu
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH -p batch
#SBATCH --qos=normal

# Export PATH of HOMER
export PATH=$PATH:/home/users/dgerard/Environment/HOMER/bin/

# Check HOMER version
perl $HOME/Environment/HOMER/configureHomer.pl -list

# It is HOMER v5.1

# Install hg 38
#perl $HOME/Environment/HOMER/configureHomer.pl -install hg38

# Run the findMotifsGenome command
findMotifsGenome.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_Fn.bed $SCRATCH/HOMER_hg38.fa $SCRATCH/HOMER/Sophia/results/FNpeaks -size 1000

findMotifsGenome.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_CTRL.bed $SCRATCH/HOMER_hg38.fa $SCRATCH/HOMER/Sophia/results/CTRLpeaks -size 1000
