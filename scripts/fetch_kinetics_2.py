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
df = pl.read_csv("../data/ob006_true_mutations_uuid.bed", separator="\t", null_values = ignore).drop_nulls()


# path for kinetics data folder of ob006
path = '~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics'
# load the kinetics data for ob006 with pysam


def get_kinetic_data(df, i):
    kinetics = pysam.AlignmentFile(os.path.expanduser(f'{path}/ob006_kinetics_diploid.bam'), 'rb')
    context = 4
    row = df.row(i)
    strand = row[df.columns.index('strand')]
    contig = row[df.columns.index('contig')]
    start = int(row[df.columns.index('start')])
    end = int(row[df.columns.index('end')])
    position_in_read = int(row[df.columns.index('position_in_read')])
    read_name = row[df.columns.index('read_name')]

    for read in kinetics.fetch(contig=contig, start=start, end=end):
        if read.query_name == read_name:
            if read.is_reverse:
                ipd_fwd = list(read.get_tag('ri')[position_in_read-context:position_in_read+context+1])
                ipd_rev = list(read.get_tag('fi')[-(position_in_read+1)-context:-(position_in_read+1)+context+1])
                base_pairs = list(read.query_sequence[(position_in_read-context):(position_in_read+context+1)])
            else:
                ipd_fwd = list(read.get_tag('fi')[position_in_read-context:position_in_read+context+1])
                ipd_rev = list(read.get_tag('ri')[-(position_in_read+1)-context:-(position_in_read+1)+context+1])
                base_pairs = list(read.query_sequence[(position_in_read-context):(position_in_read+context+1)])
            return ipd_fwd, ipd_rev, base_pairs
    return None, None, None


 

# Initialize list to store the new column values
ipd_rev_list = []
ipd_fwd_list = []
unique_id_list = []
base_pairs_list = []

# Iterate over the DataFrame and get the kinetic data
for i in tqdm(range(1000)):
    unique_id = df['unique_id'][i]
    unique_id_list.append(unique_id)
    ipd_fwd, ipd_rev, base_pairs = get_kinetic_data(df, i)
    if ipd_fwd and ipd_rev:
        ipd_fwd_list.append(ipd_fwd)
        ipd_rev_list.append(ipd_rev)
        base_pairs_list.append(base_pairs)
    else:
        ipd_fwd_list.append(None)
        ipd_rev_list.append(None)
        base_pairs_list.append(None)

# update the dataframe
# df = df.with_columns(
#     pl.Series('ipd', ipd_list)
# )

# save the DataFrame to a csv file
# df.write_csv("ipd_var_df_8mer_all.csv")

# save the ipd_list to a csv file
# print(len(base_pairs_list))
pl.DataFrame({'unique_id':unique_id_list, 'base_pairs': base_pairs_list, 'ipd_fwd':ipd_fwd_list, 'ipd_rev':ipd_rev_list}).write_parquet("updated_ipd_list_1000.parquet")
