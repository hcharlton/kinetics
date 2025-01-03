#!/bin/bash
#SBATCH --account mutationalscanning
#SBATCH -c 16
#SBATCH --mem 32g
#SBATCH --time 10:00:00

# activate mamba environment
source activate kinetics

# run the script
python fetch_kinetics.py