import os
import polars as pl
import pysam
from tqdm import tqdm
import numpy as np
import argparse


# Define the command line arguments
parser = argparse.ArgumentParser(description='Get the IPD values (forward, reverse) for a given read_name and location in h1tg000001l')
parser.add_argument('read_name', type=str, help='Name of the read to process')
parser.add_argument('position', type=int, help='Position in the read to process')

# Parse the arguments
args = parser.parse_args()


# path for kinetics data folder of ob006
path = '~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics'
# load the kinetics data for ob006 with pysam
kinetics = pysam.AlignmentFile(os.path.expanduser(f'{path}/ob006_kinetics_diploid.bam'), 'rb')

def get_kinetic_data(read_name, position_in_read):
    contig = 'h1tg000001l'
    for read in kinetics.fetch(contig=contig):
        if read.query_name == read_name:
            cigar = read.cigarstring
            is_reverse = read.is_reverse
            if is_reverse:
                ipd_fwd = list(read.get_tag('ri')[position_in_read-4:position_in_read+5])
                ipd_rev = list(read.get_tag('fi')[-(position_in_read+1)-4:-(position_in_read+1)+5])
            else:
                ipd_fwd = list(read.get_tag('fi')[position_in_read-4:position_in_read+5])
                ipd_rev = list(read.get_tag('ri')[-(position_in_read+1)-4:-(position_in_read+1)+5])
            reference_start = read.reference_start
            return ipd_fwd, ipd_rev, reference_start, is_reverse

    return None, None, None, None

# Get the kinetic data for the specified row index
read_name, position = args.read_name, args.position
ipd_fwd, ipd_rev, reference_start, is_reverse = get_kinetic_data(read_name, position)

if ipd_fwd:
    print(f'Retrieved fwd IPD:{ipd_fwd}, rev IPD:{ipd_rev}, reference_start:{reference_start}, is_reverse:{is_reverse}')
else:
    print(f'No IPD data found')



            # if is_reverse:
            #     ipd_fwd = read.get_tag('ri')[position_in_read]
            #     ipd_rev = read.get_tag('fi')[-(position_in_read+1)]
            # else:
            #     ipd_fwd = read.get_tag('fi')[position_in_read]
            #     ipd_rev = read.get_tag('ri')[-(position_in_read+1)]