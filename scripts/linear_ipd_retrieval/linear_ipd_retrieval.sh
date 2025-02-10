#!/bin/bash
#SBATCH --account mutationalscanning
#SBATCH -c 32
#SBATCH --mem 256g

source $(conda info --base)/etc/profile.d/conda.sh
conda activate kinetics_gwf
python fetch_kinetics_2.py