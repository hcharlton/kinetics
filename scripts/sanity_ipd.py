import os
import polars as pl
import pysam
from tqdm import tqdm
import numpy as np
import argparse


# Define the command line arguments
parser = argparse.ArgumentParser(description='Get the IPD values (forward, reverse) for a given row in the ob006_true_mutations.bed file')
parser.add_argument('row_index', type=int, help='Index of the row to process')
parser.add_argument('read_name', type=str, help='Name of the read to process')

# Parse the arguments
args = parser.parse_args()

# these are repeated in rows and need to be discarded (why are they repeated???)
ignore = ['contig',
        'start',
        'end',
        'ref',
        'alt',
        'context',
        'oriented_ref',
        'oriented_alt',
        'oriented_context',
        'mutation_qual',
        'strand',
        'reads_with_alt',
        'coverage',
        'sm',
        'sx',
        'read_name',
        'avg_read_lengths',
        'position_in_read']

# load the df with observed mutations
df = pl.read_csv("ob006_true_mutations.bed", separator="\t", null_values = ignore).drop_nulls()


# path for kinetics data folder of ob006
path = '~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics'
# load the kinetics data for ob006 with pysam
kinetics = pysam.AlignmentFile(os.path.expanduser(f'{path}/ob006_kinetics_diploid.bam'), 'rb')

def get_kinetic_data(df, i):
    row = df.row(i)
    contig = row[df.columns.index('contig')]
    start = int(row[df.columns.index('start')])
    end = int(row[df.columns.index('end')])
    position_in_read = int(row[df.columns.index('position_in_read')])
    read_name = row[df.columns.index('read_name')]

    for read in kinetics.fetch(contig=contig):
        if read.query_name == read_name:
            ipd_fi = read.get_tag('fi')[position_in_read]
            ipd_ri = read.get_tag('ri')[position_in_read]
            return ipd_fi, ipd_ri, contig, end

    return None, None, None, None

# Get the kinetic data for the specified row index
row_index = args.row_index
ipd_fi, ipd_ri, contig, end = get_kinetic_data(df, row_index)

if ipd_fi:
    print(f'Retrieved fi:{ipd_fi}, ri:{ipd_ri}\nVerify at contig:{contig}, end:{end}')
else:
    print(f'No IPD data found for row {row_index}')