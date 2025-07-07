---
title: "HOMER_TFmotifEnrichment"
author: Deborah Gérard^[University of Luxembourg - FSTM - DLSM - Systems Biology group - Epigenetics team]
date: "07 July, 2025"
output: 
  html_document:
    keep_md: true
    df_print: paged
    toc: true
    toc_depth: 6
  pdf_document: default
editor_options:
  chunk_output_type: console
---

This is the script for running [HOMER](https://www.sciencedirect.com/science/article/pii/S1097276510003667?via%3Dihub) on Sophia's data. The goal is to find transcription factor motifs that are over-represented between FUSO treated cells and controls from CUT&RUN data.  

#### 1. Extract the coordinates of the differential peaks that have been called by Aurélien.  
Use a FDR strictly less than 0.05
Load necessary libraries

``` r
library("tidyverse")
library("here")
```

Load data. Those are from CUT&RUN experiment and represent differential peaks between FUSO and treated cells.

``` r
# Load data
dat = read_delim(here("data",
                      "2025-02-11_cut_n_tag_peaks_RNA-seq_Fn_vs_Ctr.tsv"),
                 delim = "\t",
                 col_names = TRUE)

# Pull out the coordinates of peaks with a adjusted p-value less than 0.05
dat %>% 
  dplyr::filter(mcols.FDR < 0.05) %>% 
  dplyr::select(mcols.peakID, seqnames, start, end, strand) %>% 
  distinct(mcols.peakID, .keep_all = TRUE) %>%   # 2 peaks identifiers are duplicated (merged_70676 and merged_70691)
  write_delim(.,
              here("results",
                   "HOMER_peak_file.txt"),
                   delim = "\t",
                   col_names = FALSE)

# Make a BED file as well
dat %>% 
  dplyr::filter(mcols.FDR < 0.05) %>% 
  dplyr::mutate(start = start - 1,
                X.add = "") %>% # BED files are 0-based and need an extra empty column
  dplyr::select(seqnames, start, end, mcols.peakID, X.add, strand) %>% 
  distinct(mcols.peakID, .keep_all = TRUE) %>%   # 2 peaks identifiers are duplicated (merged_70676 and merged_70691)
  write_delim(.,
              here("results",
                   "HOMER_peak_file.bed"),
                   delim = "\t",
                   col_names = FALSE)
```

Download HOMER on HPC (iris)

``` bash
# Connect
ssh iris-cluster

# Book some ressources
si

# Download the installation file
cd $HOME/Environment
mkdir HOMER
cd HOMER
wget http://homer.ucsd.edu/homer/configureHomer.pl 

# And install
perl ./configureHomer.pl -install
export PATH=$PATH:/home/users/dgerard/Environment/HOMER/bin/

# Make a directory where to store the results
mkdir -p $SCRATCH/HOMER/Sophia

# Make a results directory
cd $SCRATCH/HOMER/Sophia
mkdir results
```

Transfer the peak file to HPC

``` bash
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_peak_file.txt iris-cluster:/scratch/users/dgerard/HOMER/Sophia/

# Do the same for BED
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_peak_file.bed iris-cluster:/scratch/users/dgerard/HOMER/Sophia/
```

And the script

``` bash
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/scripts/HOMER_CUT.TAG_findMotifsGenome.sh iris-cluster:/scratch/users/dgerard/HOMER/Sophia/
```

Copy the human genome version from HOMER into SCRATCH

``` bash
cp $HOME/Environment/HOMER/data/genomes/hg38/genome.fa $SCRATCH/HOMER_hg38.fa
```

##### 1.1 Use HOMER to perform motifs enrichment analysis.
Launch the script in passive mode

``` bash
sbatch $SCRATCH/HOMER/Sophia/HOMER_CUT.TAG_findMotifsGenome.sh
```

Display the script

``` bash
cat /Volumes/deborah.gerard/Documents/Sophia_Croce/scripts/HOMER_CUT.TAG_findMotifsGenome.sh
```

```
## #!/bin/bash -l
## #SBATCH -J HOMER_CUT&TAG_findMotifsGenome
## #SBATCH --mail-type=begin,end,fail
## #SBATCH --mail-user=deborah.gerard@uni.lu
## #SBATCH -N 1
## #SBATCH --time=48:00:00
## #SBATCH -p batch
## #SBATCH --qos=normal
## 
## # Export PATH of HOMER
## export PATH=$PATH:/home/users/dgerard/Environment/HOMER/bin/
## 
## # Check HOMER version
## perl $HOME/Environment/HOMER/configureHomer.pl -list
## 
## # It is HOMER v5.1
## 
## # Install hg 38
## #perl $HOME/Environment/HOMER/configureHomer.pl -install hg38
## 
## # Run the findMotifsGenome command
## findMotifsGenome.pl $SCRATCH/HOMER/Sophia/HOMER_peak_file.bed $SCRATCH/HOMER_hg38.fa $SCRATCH/HOMER/Sophia/results -size 1000
```

Download the results

``` bash
rsync -avzPhu iris-cluster:/scratch/users/dgerard/HOMER/Sophia/results ~/Desktop
```

Repeat the TF motifs analysis but this time using peaks that more or less present in a condition as separate list
Extract differential peaks that are more present in the Fn condition. Save as BED file

``` r
dat %>% 
  dplyr::filter(mcols.FDR < 0.05,
                mcols.logFC >= 1.0) %>% 
  dplyr::mutate(start = start - 1,
                X.add = "") %>% # BED files are 0-based and need an extra empty column
  dplyr::select(seqnames, start, end, mcols.peakID, X.add, strand) %>% 
  distinct(mcols.peakID, .keep_all = TRUE) %>%   # 2 peaks identifiers are duplicated (merged_70676 and merged_70691)
  write_delim(.,
              here("results/HOMER_MotifEnrichment_AF_FNpeaks",
                   "HOMER_peak_in_AF_Fn.bed"),
                   delim = "\t",
                   col_names = FALSE)
```

Extract peaks that are more present in the control condition. Save as BED file

``` r
dat %>% 
  dplyr::filter(mcols.FDR < 0.05,
                mcols.logFC <= -1.0) %>% 
  dplyr::mutate(start = start - 1,
                X.add = "") %>% # BED files are 0-based and need an extra empty column
  dplyr::select(seqnames, start, end, mcols.peakID, X.add, strand) %>% 
  distinct(mcols.peakID, .keep_all = TRUE) %>%   # 2 peaks identifiers are duplicated (merged_70676 and merged_70691)
  write_delim(.,
              here("results/HOMER_MotifEnrichment_AF_CTRLpeaks",
                   "HOMER_peak_in_AF_CTRL.bed"),
                   delim = "\t",
                   col_names = FALSE)
```

Create 2 new directories where to store the results 

``` bash
# Connect to HPC
ssh iris-cluster

# For peaks only represented in the FN condition
mkdir $SCRATCH/HOMER/Sophia/results/FNpeaks

# For peaks only represented in the CTRL condition
mkdir $SCRATCH/HOMER/Sophia/results/CTRLpeaks
```

Transfer the 2 new BED files to HPC. Same for the launcher script

``` bash
# For the peaks in the FN condition
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_AF_FNpeaks/HOMER_peak_in_AF_Fn.bed iris-cluster:/scratch/users/dgerard/HOMER/Sophia/

# For the peaks in the CTRL condition
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_AF_CTRLpeaks/HOMER_peak_in_AF_CTRL.bed iris-cluster:/scratch/users/dgerard/HOMER/Sophia/

# Script
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/scripts/HOMER_CUT.TAG_findMotifsGenome_FN_CTRL_peaks.sh iris-cluster:/scratch/users/dgerard/HOMER/Sophia/
```

Launch the script

``` bash
sbatch $SCRATCH/HOMER/Sophia/HOMER_CUT.TAG_findMotifsGenome_FN_CTRL_peaks.sh
```

Display the script

``` bash
cat /Volumes/deborah.gerard/Documents/Sophia_Croce/scripts/HOMER_CUT.TAG_findMotifsGenome_FN_CTRL_peaks.sh
```

```
## #!/bin/bash -l
## #SBATCH -J HOMER_CUT&TAG_findMotifsGenome
## #SBATCH --mail-type=begin,end,fail
## #SBATCH --mail-user=deborah.gerard@uni.lu
## #SBATCH -N 1
## #SBATCH --time=48:00:00
## #SBATCH -p batch
## #SBATCH --qos=normal
## 
## # Export PATH of HOMER
## export PATH=$PATH:/home/users/dgerard/Environment/HOMER/bin/
## 
## # Check HOMER version
## perl $HOME/Environment/HOMER/configureHomer.pl -list
## 
## # It is HOMER v5.1
## 
## # Install hg 38
## #perl $HOME/Environment/HOMER/configureHomer.pl -install hg38
## 
## # Run the findMotifsGenome command
## findMotifsGenome.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_Fn.bed $SCRATCH/HOMER_hg38.fa $SCRATCH/HOMER/Sophia/results/FNpeaks -size 1000
## 
## findMotifsGenome.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_CTRL.bed $SCRATCH/HOMER_hg38.fa $SCRATCH/HOMER/Sophia/results/CTRLpeaks -size 1000
```

The results are back. Download them locally and then on ATLAS

``` bash
# The results for the FN peaks
rsync -avzPhu iris-cluster:/scratch/users/dgerard/HOMER/Sophia/results/FNpeaks/ ~/Desktop/HOMER_MotifEnrichment_AF_FNpeaks/

# The results for the CTRL peaks
rsync -avzPhu iris-cluster:/scratch/users/dgerard/HOMER/Sophia/results/CTRLpeaks/ ~/Desktop/HOMER_MotifEnrichment_AF_CTRLpeaks/
```

#### 2. Annotate (meaning in this case finding the positions of the different motifs) the peaks with the motifs

``` bash
# Remake the directory
mkdir -p $SCRATCH/HOMER/Sophia

# Trasnfer back the BED files
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_DiffPeaks/HOMER_peak_file.bed iris-cluster:/scratch/users/dgerard/HOMER/Sophia/

rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_AF_FNpeaks/HOMER_peak_in_AF_Fn.bed iris-cluster:/scratch/users/dgerard/HOMER/Sophia/

rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_AF_CTRLpeaks/HOMER_peak_in_AF_CTRL.bed iris-cluster:/scratch/users/dgerard/HOMER/Sophia/

# Transfer back the motifs files of each condition
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_DiffPeaks/homerMotifs.all.motifs iris-cluster:/scratch/users/dgerard/HOMER/Sophia/diffPeaks/

rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_AF_FNpeaks/homerMotifs.all.motifs iris-cluster:/scratch/users/dgerard/HOMER/Sophia/FN/

rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_AF_CTRLpeaks/homerMotifs.all.motifs iris-cluster:/scratch/users/dgerard/HOMER/Sophia/CTRL/
```

Transfer the script from ATLAS to HPC and run it

``` bash
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/scripts/HOMER_CUT.TAG_Annotation.sh iris-cluster:/scratch/users/dgerard/

# Run
sbatch $SCRATCH/HOMER_CUT.TAG_Annotation.sh
```

Display the script

``` bash
cat /Volumes/deborah.gerard/Documents/Sophia_Croce/scripts/HOMER_CUT.TAG_Annotation.sh
```

```
## #!/bin/bash -l
## #SBATCH -J HOMER_CUT&TAG_AnnotatePeaks
## #SBATCH --mail-type=begin,end,fail
## #SBATCH --mail-user=deborah.gerard@uni.lu
## #SBATCH -N 1
## #SBATCH --time=03:00:00
## #SBATCH -p batch
## #SBATCH --exclusive
## #SBATCH --qos=normal
## 
## # Export PATH of HOMER
## export PATH=$PATH:/home/users/dgerard/Environment/HOMER/bin/
## 
## # Check HOMER version
## perl $HOME/Environment/HOMER/configureHomer.pl -list
## 
## # It is HOMER v5.1
## 
## # Run the annotatePeaks command
## ## On the differential peaks
## annotatePeaks.pl $SCRATCH/HOMER/Sophia/HOMER_peak_file.bed \
## hg38 \
## -m $SCRATCH/HOMER/Sophia/diffPeaks/all.diffPeaks.known.motifs.motif > $SCRATCH/HOMER/Sophia/diffPeaks/diff_peaks_motif_anno.txt
## 
## ## On the peaks UP in Fuso
## annotatePeaks.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_Fn.bed \
## hg38 \
## -m $SCRATCH/HOMER/Sophia/FN/all.FN.known.motifs.motif > $SCRATCH/HOMER/Sophia/FN/FN_peaks_motif_anno.txt
## 
## ## On the peaks UP in CTRL
## annotatePeaks.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_CTRL.bed \
## hg38 \
## -m $SCRATCH/HOMER/Sophia/CTRL/all.CTRL.known.motifs.motif > $SCRATCH/HOMER/Sophia/CTRL/CTRL_peaks_motif_anno.txt
```

Transfer the results back

``` bash
rsync -avzPhu --exclude homerMotifs.all.motifs iris-cluster:/scratch/users/dgerard/HOMER/Sophia/diffPeaks ~/Desktop/

rsync -avzPhu --exclude homerMotifs.all.motifs iris-cluster:/scratch/users/dgerard/HOMER/Sophia/FN ~/Desktop/

rsync -avzPhu --exclude homerMotifs.all.motifs iris-cluster:/scratch/users/dgerard/HOMER/Sophia/CTRL ~/Desktop/
```

That run was actually performed with the list of de novo motifs. Concatenate all known motifs together and rerun

``` bash
# For the differential peaks
cd /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_DiffPeaks/knownResults/
cat *.motif > all.diffPeaks.known.motifs.motif

# For the peaks in FN condition
cd /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_AF_FNpeaks/knownResults/
cat *.motif > all.FN.known.motifs.motif

# For the peaks in the CTRL condition
cd /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_AF_CTRLpeaks/knownResults/
cat *.motif > all.CTRL.known.motifs.motif

# Transfer back the motifs files of each condition
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_DiffPeaks/knownResults/all.diffPeaks.known.motifs.motif iris-cluster:/scratch/users/dgerard/HOMER/Sophia/diffPeaks/

rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_AF_FNpeaks/knownResults/all.FN.known.motifs.motif iris-cluster:/scratch/users/dgerard/HOMER/Sophia/FN/

rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/results/HOMER_MotifEnrichment_AF_CTRLpeaks/knownResults/all.CTRL.known.motifs.motif iris-cluster:/scratch/users/dgerard/HOMER/Sophia/CTRL/
```

Transfer the script from ATLAS to HPC and run it

``` bash
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/scripts/HOMER_CUT.TAG_Annotation.sh iris-cluster:/scratch/users/dgerard/

# Run
sbatch $SCRATCH/HOMER_CUT.TAG_Annotation.sh
```

Display the script

``` bash
cat /Volumes/deborah.gerard/Documents/Sophia_Croce/scripts/HOMER_CUT.TAG_Annotation.sh
```

```
## #!/bin/bash -l
## #SBATCH -J HOMER_CUT&TAG_AnnotatePeaks
## #SBATCH --mail-type=begin,end,fail
## #SBATCH --mail-user=deborah.gerard@uni.lu
## #SBATCH -N 1
## #SBATCH --time=03:00:00
## #SBATCH -p batch
## #SBATCH --exclusive
## #SBATCH --qos=normal
## 
## # Export PATH of HOMER
## export PATH=$PATH:/home/users/dgerard/Environment/HOMER/bin/
## 
## # Check HOMER version
## perl $HOME/Environment/HOMER/configureHomer.pl -list
## 
## # It is HOMER v5.1
## 
## # Run the annotatePeaks command
## ## On the differential peaks
## annotatePeaks.pl $SCRATCH/HOMER/Sophia/HOMER_peak_file.bed \
## hg38 \
## -m $SCRATCH/HOMER/Sophia/diffPeaks/all.diffPeaks.known.motifs.motif > $SCRATCH/HOMER/Sophia/diffPeaks/diff_peaks_motif_anno.txt
## 
## ## On the peaks UP in Fuso
## annotatePeaks.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_Fn.bed \
## hg38 \
## -m $SCRATCH/HOMER/Sophia/FN/all.FN.known.motifs.motif > $SCRATCH/HOMER/Sophia/FN/FN_peaks_motif_anno.txt
## 
## ## On the peaks UP in CTRL
## annotatePeaks.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_CTRL.bed \
## hg38 \
## -m $SCRATCH/HOMER/Sophia/CTRL/all.CTRL.known.motifs.motif > $SCRATCH/HOMER/Sophia/CTRL/CTRL_peaks_motif_anno.txt
```

Transfer the results back

``` bash
rsync -avzPhu --exclude all.diffPeaks.known.motifs.motif iris-cluster:/scratch/users/dgerard/HOMER/Sophia/diffPeaks ~/Desktop/

rsync -avzPhu --exclude all.FN.known.motifs.motif iris-cluster:/scratch/users/dgerard/HOMER/Sophia/FN ~/Desktop/

rsync -avzPhu --exclude all.CTRL.known.motifs.motif iris-cluster:/scratch/users/dgerard/HOMER/Sophia/CTRL ~/Desktop/
```

#### 3. Make BED files of peaks having ELF3 motif

``` r
# For the Fn condition
FN.anno.p = read_delim(here("results/HOMER_MotifEnrichment_AF_FNpeaks/Peak_annotation/",
                            "FN_peaks_motif_anno.txt"),
                       delim = "\t",
                       col_names = TRUE)

# Select the ELF3 motif
FN.anno.p %>% 
  rename(PeakID = `PeakID (cmd=annotatePeaks.pl /scratch/users/dgerard//HOMER/Sophia/HOMER_peak_in_AF_Fn.bed hg38 -m /scratch/users/dgerard//HOMER/Sophia/FN/all.FN.known.motifs.motif)`) %>% 
  select(Chr, Start, End, PeakID,matches("ELF3")) %>% 
  rename(ELF3.motif = `ELF3(ETS)/PDAC-ELF3-ChIP-Seq(GSE64557)/Homer Distance From Peak(sequence,strand,conservation)`) %>% 
  drop_na() %>% 
  select(-ELF3.motif) %>% 
  write_delim(here("data/BED",
                   "FN_UP_peaks_ELF3_motifs.bed"),
              col_names = FALSE,
              delim = "\t")

# Same for CTRL
CTRL.anno.p = read_delim(here("results/HOMER_MotifEnrichment_AF_CTRLpeaks/Peak_annotation/",
                            "CTRL_peaks_motif_anno.txt"),
                       delim = "\t",
                       col_names = TRUE)

# Select the ELF3 motif
CTRL.anno.p %>% 
  rename(PeakID = `PeakID (cmd=annotatePeaks.pl /scratch/users/dgerard//HOMER/Sophia/HOMER_peak_in_AF_CTRL.bed hg38 -m /scratch/users/dgerard//HOMER/Sophia/CTRL/all.CTRL.known.motifs.motif)`) %>% 
  select(Chr, Start, End, PeakID,matches("ELF3")) %>% 
  rename(ELF3.motif = `ELF3(ETS)/PDAC-ELF3-ChIP-Seq(GSE64557)/Homer Distance From Peak(sequence,strand,conservation)`) %>% 
  drop_na() %>% 
  select(-ELF3.motif) %>% 
  write_delim(here("data/BED",
                   "CTRL_UP_peaks_ELF3_motifs.bed"),
              col_names = FALSE,
              delim = "\t")
```

#### 4. Make a PCA plot using deeptools

``` bash
# On HPC (Aion)
ssh aion-cluster

# Book some ressources to install deeptols
salloc -p interactive -N 1 -c 2 --qos debug -C batch -t 01:00:00

# Make separate Python environment
# Load modules to use Python3
module load env/development/2024a
module load lang/Python/3.12.3-GCCcore-13.3.0

# Create the environment in my home directory with ther other environments
python3 -m venv ~/Environment/deeptools

# Activate the newly created environment
source $HOME/Environment/deeptools/bin/activate

# Install via pip
pip install deeptools

# Check the version of deeptools
deeptools --version
```
The version of deeptools is 3.5.6

Run first `multiBigwigSummary` command on bigwig files and then `plotPCA`

``` bash
# Make a directory where to store the bigwig
mkdir -p $SCRATCH/deeptools/bigwig

# Copy the bigwigs to HPC
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/data/BIGWIG/*.bw aion-cluster:/scratch/users/dgerard/deeptools/bigwig/

# Transfer the launcher script to HPC
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia/scripts/DEEPTOOLS_PCA.sh aion-cluster:/scratch/users/dgerard/deeptools/

# Run the launcher script
sbatch $SCRATCH/deeptools/DEEPTOOLS_PCA.sh
```

Display the laucher scripts

``` bash
cat /Volumes/deborah.gerard/Documents/Sophia_Croce/scripts/DEEPTOOLS_PCA.sh
```

```
## #!/bin/bash -l
## #SBATCH -J deeptools_PCA
## #SBATCH --mail-type=begin,end,fail
## #SBATCH --mail-user=deborah.gerard@uni.lu
## #SBATCH -N 1
## #SBATCH --time=06:00:00
## #SBATCH -p batch
## #SBATCH --exclusive
## #SBATCH --qos=normal
## 
## # Activate the python environment for deeptools
## source $HOME/Environment/deeptools/bin/activate
## 
## # Check deeptools version
## deeptools --version
## 
## # Run first multiBigwigSummary command on bigwig files
## multiBigwigSummary bins -b $SCRATCH/deeptools/bigwig/*.bw \
## -o $SCRATCH/deeptools/results/AF_FN_matrix.npz \
## -p 24
## 
## # Run then the plotPCA command
## plotPCA -in $SCRATCH/deeptools/results/AF_FN_matrix.npz \
## -o $SCRATCH/deeptools/results/AF_FN_CUTandTAG_PCA.pdf \
## --transpose \
## --outFileNameData $SCRATCH/deeptools/results/AF_FN_CUTandTAG_PCA.tab
## 
## deactivate
```

Download the results back

``` bash
rsync -avzPhu aion-cluster:/scratch/users/dgerard/deeptools/results/ ~/Desktop/deeptools/
```

#### 5. Results of the PCA plot
From the PCA plot, one can see that the pair AF_ctr_2-1/AF_Fn_2-1 does not cluster well. It has been removed and the differential peak analysis performed by Aurélien was rerun without that pair.  
Now, extract the differential peaks (with a FDR < 0.05) again

``` r
# Load new data
dat = read_delim(here("data",
                      "2025-07-07_cut_n_tag_peaks_RNA-seq_Fn_vs_Ctr_3reps.tsv"),
                 delim = "\t",
                 col_names = TRUE)

# Make a BED file 
dat %>% 
  dplyr::filter(mcols.FDR < 0.05) %>% 
  dplyr::mutate(start = start - 1,
                X.add = "") %>% # BED files are 0-based and need an extra empty column
  dplyr::select(seqnames, start, end, mcols.peakID, X.add, strand) %>% 
  distinct(mcols.peakID, .keep_all = TRUE) %>%
  write_delim(.,
              here("results/HOMER_MotifEnrichment_DiffPeaks/RunWith3RepsPerCondition",
                   "HOMER_peak_file_3repsPerCondition.bed"),
                   delim = "\t",
                   col_names = FALSE)
```
**Conclusion**: with 3 replicates per condition, we have 2535 significant differential peaks using a FDR < 0.05 while there were 3607 using all replicates.  
Extract also the peaks that are "up" in the Fn condition and the ones that are "up" in the control condition

``` r
# up in the Fn condition (FDR < 0.05 and logFC >= 1.0)
dat %>% 
  dplyr::filter(mcols.FDR < 0.05,
                mcols.logFC >= 1.0) %>% 
  dplyr::mutate(start = start - 1,
                X.add = "") %>% # BED files are 0-based and need an extra empty column
  dplyr::select(seqnames, start, end, mcols.peakID, X.add, strand) %>% 
  distinct(mcols.peakID, .keep_all = TRUE) %>%
  write_delim(.,
              here("results/HOMER_MotifEnrichment_AF_FNpeaks/RunWith3RepsPerCondition",
                   "HOMER_peak_in_AF_Fn_3repsPerCondition.bed"),
                   delim = "\t",
                   col_names = FALSE)

# up in the control condition (FDR < 0.05 and logFC >= 1.0)
dat %>% 
  dplyr::filter(mcols.FDR < 0.05,
                mcols.logFC <= -1.0) %>% 
  dplyr::mutate(start = start - 1,
                X.add = "") %>% # BED files are 0-based and need an extra empty column
  dplyr::select(seqnames, start, end, mcols.peakID, X.add, strand) %>% 
  distinct(mcols.peakID, .keep_all = TRUE) %>%
  write_delim(.,
              here("results/HOMER_MotifEnrichment_AF_CTRLpeaks/RunWith3RepsPerCondition",
                   "HOMER_peak_in_AF_CTRL_3repsPerCondition.bed"),
                   delim = "\t",
                   col_names = FALSE)
```
**Conclusion**: with 3 replicates per condition, we have 1444 significant peaks using a FDR < 0.05 and logFC >= 1.0 in the FN condition using 3 replicates while there were 1841 using all replicates. For the control condition, we have 779 significant peaks using a FDR < 0.05 and logFC >= 1.0 in the control condition using 3 replicates while there were 949 using all replicates.  

Transfer the script and the BED files to Iris (HPC)

``` bash
# Script
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia_Croce/scripts/HOMER_CUT.TAG_findMotifsGenome_3reps.sh iris-cluster:/scratch/users/dgerard/

# BED files
rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia_Croce/results/HOMER_MotifEnrichment_DiffPeaks/RunWith3RepsPerCondition/HOMER_peak_file_3repsPerCondition.bed iris-cluster:/scratch/users/dgerard/HOMER/Sophia/

rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia_Croce/results/HOMER_MotifEnrichment_AF_FNpeaks/RunWith3RepsPerCondition/HOMER_peak_in_AF_Fn_3repsPerCondition.bed iris-cluster:/scratch/users/dgerard/HOMER/Sophia/

rsync -avzPhu /Volumes/deborah.gerard/Documents/Sophia_Croce/results/HOMER_MotifEnrichment_AF_CTRLpeaks/RunWith3RepsPerCondition/HOMER_peak_in_AF_CTRL_3repsPerCondition.bed iris-cluster:/scratch/users/dgerard/HOMER/Sophia/
```

Launch the script

``` bash
sbatch $SCRATCH/HOMER_CUT.TAG_findMotifsGenome_3reps.sh
```

Display the script

``` bash
cat /Volumes/deborah.gerard/Documents/Sophia_Croce/scripts/HOMER_CUT.TAG_findMotifsGenome_3reps.sh
```

```
## #!/bin/bash -l
## #SBATCH -J HOMER_CUT&TAG_findMotifsGenome_3reps
## #SBATCH --mail-type=begin,end,fail
## #SBATCH --mail-user=deborah.gerard@uni.lu
## #SBATCH -N 1
## #SBATCH --exclusive
## #SBATCH --time=04:00:00
## #SBATCH -p batch
## #SBATCH --qos=normal
## 
## # Export PATH of HOMER
## export PATH=$PATH:/home/users/dgerard/Environment/HOMER/bin/
## 
## # Check HOMER version
## perl $HOME/Environment/HOMER/configureHomer.pl -list
## 
## # It is HOMER v5.1
## 
## # Run the findMotifsGenome command on the significant differential peaks
## findMotifsGenome.pl $SCRATCH/HOMER/Sophia/HOMER_peak_file_3repsPerCondition.bed $SCRATCH/HOMER_hg38.fa $SCRATCH/HOMER/Sophia/results/diffPeaks -size 1000
## 
## # Run the findMotifsGenome command on the significant peaks up in the Fn condition
## findMotifsGenome.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_Fn_3repsPerCondition.bed $SCRATCH/HOMER_hg38.fa $SCRATCH/HOMER/Sophia/results/FN -size 1000
## 
## # Run the findMotifsGenome command on the significant peaks up in the control condition
## findMotifsGenome.pl $SCRATCH/HOMER/Sophia/HOMER_peak_in_AF_CTRL_3repsPerCondition.bed $SCRATCH/HOMER_hg38.fa $SCRATCH/HOMER/Sophia/results/CTRL -size 1000
```

