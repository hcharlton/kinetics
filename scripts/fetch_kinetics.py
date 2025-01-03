import os
import polars as pl
import pysam
from tqdm import tqdm
import numpy as np

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

    for read in kinetics.fetch(contig=contig, start=start, end=end):
        if read.query_name == read_name:
            ipd = list(read.get_tag('ri')[position_in_read-4:position_in_read+5])
            return ipd

    return None


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

# Initialize list to store the new column values
ipd_list = []

# Iterate over the DataFrame and get the kinetic data
for i in tqdm(range(10)):
    ipd = get_kinetic_data(df, i)
    ipd_list.append(ipd)

# update the dataframe
# df = df.with_columns(
#     pl.Series('ipd', ipd_list)
# )

# save the DataFrame to a csv file
# df.write_csv("ipd_var_df_8mer_all.csv")

# save the ipd_list to a csv file

pl.DataFrame(ipd_list).write_csv("ipd_8context_all.csv", include_header=False)
