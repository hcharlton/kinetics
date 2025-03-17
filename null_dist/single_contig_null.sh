#!/bin/bash
#SBATCH --account mutationalscanning
#SBATCH --time 10:00:00
#SBATCH -c 32
#SBATCH --mem 512g

source $(conda info --base)/etc/profile.d/conda.sh
conda activate kinetics_gwf
python retrieve_null_batch.py