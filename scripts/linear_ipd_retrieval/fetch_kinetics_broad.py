import os
import polars as pl
import pysam
from tqdm import tqdm


# load the bed file with observed mutations into a DF
df = pl.read_csv("../../data/ob006_true_mutations_uuid.bed", separator="\t").drop_nulls()
print(df.head())

# path for kinetics data folder of ob006
path = '~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics'
# load the kinetics data for ob006 with pysam
context = 16
def get_kinetic_data(df, i, context):
    kinetics = pysam.AlignmentFile(os.path.expanduser(f'{path}/ob006_kinetics_diploid.bam'), 'rb')
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
for i in tqdm(range(len(df))):
    unique_id = df['unique_id'][i]
    unique_id_list.append(unique_id)
    try:
        ipd_fwd, ipd_rev, base_pairs = get_kinetic_data(df, i)
    except BaseException:
        ipd_fwd = None
        ipd_rev = None
        base_pairs = None
    ipd_fwd_list.append(ipd_fwd)
    ipd_rev_list.append(ipd_rev)
    base_pairs_list.append(base_pairs)  

    # if ipd_fwd and ipd_rev:
    #     ipd_fwd_list.append(ipd_fwd)
    #     ipd_rev_list.append(ipd_rev)
    #     base_pairs_list.append(base_pairs)
    # else:
    #     ipd_fwd_list.append(None)
    #     ipd_rev_list.append(None)
    #     base_pairs_list.append(None)


def create_header(context):
    return (
        ['unique_id', 'base_pairs']
        + [f'ipd_fwd_a{i}' for i in range(context - 1, 0, -1)]
        + ['ipd_fwd_event']
        + [f'ipd_fwd_p{i}' for i in range(1, context)]
        + [f'ipd_rev_a{i}' for i in range(context - 1, 0, -1)]
        + ['ipd_rev_event']
        + [f'ipd_rev_p{i}' for i in range(1, context)]
    )

header = create_header(context)

# make np.array of the lists
data = [unique_id_list, base_pairs_list]

