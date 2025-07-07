#!/bin/bash -l
#SBATCH -J HOMER_CUT&TAG_findMotifsGenome_3reps
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=deborah.gerard@uni.lu
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH -p batch
#SBATCH --qos=normal

# Export PATH of HOMER
export PATH=$PATH:/home/users/dgerard/Environment/HOMER/bin/

# Check HOMER version
perl $HOME/Environment/HOMER/configureHomer.pl -list

# It is HOMER v5.1

# Run the findMotifsGenome command on the significant differential peaks
findMotifsGenome.pl $SCRATCH/HOMER/Sophia/HOMER_peak_file_3repsPerCondition.bed $SCRATCH/HOMER_hg38.fa $SCRATCH/HOMER/Sophia/results/diffPeaks -size 1000

# Run the findMotifsGenome command on the significant peaks up in the Fn condition
findMotifsGenome.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_Fn_3repsPerCondition.bed $SCRATCH/HOMER_hg38.fa $SCRATCH/HOMER/Sophia/results/FN -size 1000

# Run the findMotifsGenome command on the significant peaks up in the control condition
findMotifsGenome.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_CTRL_3repsPerCondition.bed $SCRATCH/HOMER_hg38.fa $SCRATCH/HOMER/Sophia/results/CTRL -size 1000

