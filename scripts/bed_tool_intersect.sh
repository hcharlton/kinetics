#!/bin/bash
#SBATCH --account mutationalscanning
#SBATCH -c 32
#SBATCH --mem 256g

source $(conda info --base)/etc/profile.d/conda.sh
conda activate kinetics_gwf
bedtools intersect -a ../data/ob006_true_mutations_cleaned_noheader.bed -b ~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics/ob006_kinetics_diploid.bam -wb > ob006_intersected_reads.bed