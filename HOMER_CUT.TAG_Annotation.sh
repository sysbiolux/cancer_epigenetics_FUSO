#!/bin/bash -l
#SBATCH -J HOMER_CUT&TAG_AnnotatePeaks
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=deborah.gerard@uni.lu
#SBATCH -N 1
#SBATCH --time=03:00:00
#SBATCH -p batch
#SBATCH --exclusive
#SBATCH --qos=normal

# Export PATH of HOMER
export PATH=$PATH:/home/users/dgerard/Environment/HOMER/bin/

# Check HOMER version
perl $HOME/Environment/HOMER/configureHomer.pl -list

# It is HOMER v5.1

# Run the annotatePeaks command
## On the differential peaks
annotatePeaks.pl $SCRATCH/HOMER/Sophia/HOMER_peak_file.bed \
hg38 \
-m $SCRATCH/HOMER/Sophia/diffPeaks/all.diffPeaks.known.motifs.motif > $SCRATCH/HOMER/Sophia/diffPeaks/diff_peaks_motif_anno.txt

## On the peaks UP in Fuso
annotatePeaks.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_Fn.bed \
hg38 \
-m $SCRATCH/HOMER/Sophia/FN/all.FN.known.motifs.motif > $SCRATCH/HOMER/Sophia/FN/FN_peaks_motif_anno.txt

## On the peaks UP in CTRL
annotatePeaks.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_CTRL.bed \
hg38 \
-m $SCRATCH/HOMER/Sophia/CTRL/all.CTRL.known.motifs.motif > $SCRATCH/HOMER/Sophia/CTRL/CTRL_peaks_motif_anno.txt

