import os
import polars as pl
import pysam
from tqdm import tqdm
import logging
import configparser

### load config ###
config_path = '~/mutationalscanning/Workspaces/chcharlton/kinetics/scripts/linear_ipd_retrieval/test.ini'

def load_config(config_path):
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found at: {config_path}")
    
    config = configparser.ConfigParser()
    config.read(config_path)

    # paths existence
    if 'Paths' not in config:
        raise ValueError('Paths sections not found in config')
    # path to bam file existence
    if 'bam_filepath' not in config['Paths']:
        raise ValueError('no bam_file specified in config')
    # path bed file existence
    if 'bed_filepath' not in config['Paths']:
        raise ValueError('no bed_filepath specified in config')
    # constants existence
    if 'Constants' not in config:
        raise ValueError('no Constants section in config')
    # context value existence
    if 'context' not in config['Constants']:
        raise ValueError('no context specified in config')
    
    return config 

config = load_config(os.path.expanduser(config_path))

bam_filepath = os.path.expanduser(config['Paths']['bam_filepath'])
bed_filepath = os.path.expanduser(config['Paths']['bed_filepath']) 
results_filepath = os.path.expanduser(config['Paths']['output_dir'])
context = int(config['Constants']['context'])




# enable logging 
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_bed_file(bed_filepath):
    """Loads a bed file into a Polars dataframe and validates it"""
    logging.info(f"Loading BED file: {bed_filepath}")  # Log the *expanded* path
    if not os.path.exists(bed_filepath): 
        raise FileNotFoundError(f"BED file not found at: {bed_filepath}")
    try:
        df = pl.read_csv(bed_filepath, separator="\t").drop_nulls()
        required_columns = ['strand', 'contig', 'start', 'end', 'position_in_read', 'read_name', 'unique_id']
        for column in required_columns:
            assert column in df.columns, f"Column {column} not found in BED file"
        return df
    except Exception as e:
        raise

# load the bed file with observed mutations into a DF
df = load_bed_file(bed_filepath)



# load the kinetics data for ob006 with pysam
def get_kinetic_data(row: tuple, header: list, bam: pysam.AlignmentFile, context: int) -> tuple[list | None, list | None, list | None]:
    """
    retrieves selected kinetics and nt data for a single ead from an open pysam 
    alignment file
    Args:
    row: a tuple representing a single row of a Polars DataFrame
    bam: an *open* pysam.AlignmentFile object
    context: the number of bases to include on either side of the mutation (total size = 2*context + 1)
    Returns:
    a tuple: (ipd_fwd, ipd_rev, base_pairs) where:
    ipd_fwd: a list of IPD values for the forward strand
    ipd_rev: a list of IPD values for the reverse strand
    base_pairs: a list of base pairs centered on the mutation
    """
    unique_id = row[header.index('unique_id')]
    contig = row[header.index('contig')]
    start = int(row[header.index('start')])
    end = int(row[header.index('end')])
    position_in_read = int(row[header.index('position_in_read')])
    read_name = row[header.index('read_name')]
    null_data = None, None, None, None, None, None

    try:
        for read in bam.fetch(contig=contig, start=start, end=end):

            fwd_start_index = position_in_read - context
            fwd_end_index = position_in_read + context + 1
            rev_start_index = -(position_in_read + 1) - context
            rev_end_index = -(position_in_read + 1) + context + 1

            if read.query_name == read_name:
                if read.is_reverse: # note that these were determined by T&E + IGV Verification
                    fwd_tag = 'ri'
                    rev_tag = 'fi'
                else:
                    fwd_tag = 'fi'
                    rev_tag = 'ri'
                try:
                    ipd_fwd = list(read.get_tag(fwd_tag)[fwd_start_index:fwd_end_index])
                    ipd_rev = list(read.get_tag(rev_tag)[rev_start_index:rev_end_index])
                    fn = read.get_tag('fn')
                    rn = read.get_tag('rn')
                    base_pairs = list(read.query_sequence[fwd_start_index:fwd_end_index])

                    expected_length = 2 * context + 1
                    if not (len(ipd_fwd) == len(ipd_rev) == len(base_pairs) == expected_length):
                        logging.warning(f"Unexpected IPD/base length. Read: {read_name}, Expected: {expected_length}, "
                                        f"Got: fwd:{len(ipd_fwd)}, rev:{len(ipd_rev)}, bases:{len(base_pairs)}")
                        return null_data

                    return unique_id, ipd_fwd, ipd_rev, base_pairs, fn, rn
                except KeyError as e:
                    logging.warning(f"Read {read_name} is missing tag {e}. Skipping.")
                    return null_data
                except IndexError as e:
                    logging.warning(f"Index out of range for read {read_name}: {e}. Skipping.")
                    return null_data
    except ValueError as e:
        logging.error(f"Error during bam.fetch (likely BAM index issue): {e}")
        return null_data

    return null_data

 
def process_kinetics(df: pl.DataFrame, bam_filepath: str, context: int) -> pl.DataFrame:
    """Retrieve and store the kinetics data for all reads in the DataFrame."""
    if not os.path.exists(bam_filepath):
        raise FileNotFoundError(f"BAM file not found: {bam_filepath}")

    # try:
    with pysam.AlignmentFile(bam_filepath, 'rb') as bam:
        # Prepare the data for apply.  Extract only necessary columns *before* the apply.
        cols = ['unique_id', 'contig','start', 'end', 'position_in_read', 'read_name']
        subset_df = df.select(cols)
        # map rows applies the get_kinetics data function to all the r
        results = subset_df.map_rows(lambda row: get_kinetic_data(row=row, header=cols, bam=bam, context=context))
        results.columns = ['unique_id', 'ipd_fwd', 'ipd_rev', 'base_pairs', 'fn', 'rn']
        # Join with the original DataFrame based on row index.  Create index columns first.
        df = df.with_row_index('row_index')

    return results



def main():
    # try:
    df = load_bed_file(bed_filepath)
    logging.info(f"Loaded BED file: {bed_filepath}")

    final_df = process_kinetics(df, bam_filepath, context)
    logging.info("retrieval complete")
    output_filepath = os.path.join(results_filepath, os.path.basename(bed_filepath).split("_")[0]) + '.parquet'
    final_df.write_parquet(output_filepath)

    # except Exception as e:
    #     logging.error(f"An error occurred: {e}


if __name__ == "__main__":
    main()
