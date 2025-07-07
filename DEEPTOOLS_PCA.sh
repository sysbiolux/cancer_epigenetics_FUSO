#!/bin/bash -l
#SBATCH -J deeptools_PCA
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=deborah.gerard@uni.lu
#SBATCH -N 1
#SBATCH --time=06:00:00
#SBATCH -p batch
#SBATCH --exclusive
#SBATCH --qos=normal

# Activate the python environment for deeptools
source $HOME/Environment/deeptools/bin/activate

# Check deeptools version
deeptools --version

# Run first multiBigwigSummary command on bigwig files
multiBigwigSummary bins -b $SCRATCH/deeptools/bigwig/*.bw \
-o $SCRATCH/deeptools/results/AF_FN_matrix.npz \
-p 24

# Run then the plotPCA command
plotPCA -in $SCRATCH/deeptools/results/AF_FN_matrix.npz \
-o $SCRATCH/deeptools/results/AF_FN_CUTandTAG_PCA.pdf \
--transpose \
--outFileNameData $SCRATCH/deeptools/results/AF_FN_CUTandTAG_PCA.tab

deactivate