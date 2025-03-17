import pysam, logging, os
from tqdm import tqdm
import polars as pl

def load_bed_file(bed_path: str) -> list[tuple[str, int, int]]:
    regions = []
    with open(bed_path, "r") as bed:
        for line in bed:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue
            contig = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            regions.append((contig, start, end))
    return regions

def process_read_for_null(read, window_size: int = 64) -> list[dict]:
    # define the fwd and reverse ipd and pulse width tags depending on read orientation
    fwd_ipd_tag, rev_ipd_tag = (("ri", "fi") if read.is_reverse else ("fi", "ri"))
    fwd_pw_tag, rev_pw_tag = (("rp", "fp") if read.is_reverse else ("fp", "rp"))
    seq = read.query_sequence
    qual = list(read.query_qualities)
    try:
        fn = int(read.get_tag('fn'))
        rn = int(read.get_tag('rn'))
        sm_all = list(read.get_tag('sm'))
        sx_all = list(read.get_tag('sx'))
        fwd_ipd_all = list(read.get_tag(fwd_ipd_tag))
        rev_ipd_all = list(read.get_tag(rev_ipd_tag))
        fwd_pw_all = list(read.get_tag(fwd_pw_tag))
        rev_pw_all = list(read.get_tag(rev_pw_tag))
        # NB: ri is indexed from tail (see PacBio docs for fi/ri)
        # https://pacbiofileformats.readthedocs.io/en/13.0/BAM.html
        rev_ipd_all.reverse()
        rev_pw_all.reverse()
    except Exception as e:
        logging.warning(f"Error getting kinetics tags for {read.query_name}: {e}")
        return []
    if len(seq) < window_size:
        logging.info(f"Read {read.query_name} is shorter than window size ({window_size}).")
        return []
    # Get reference positions for each base in the read.
    # using full_length=True to account for indels
    ref_positions = read.get_reference_positions(full_length=True)
    rows = []
    num_windows = len(seq) // window_size
    # skip the first and last two windows since the beginnings and ends of reads are unreliable
    for i in range(2, num_windows-2):
        win_start = i * window_size
        win_end = win_start + window_size
        if win_end > len(seq):
            continue
        # define the center of the window:
        # use the three bases starting at index win_start+31 (i.e. indices 31,32,33)
        center_start = win_start + 31
        center_end = center_start + 3
        # calculate window center
        read_center = center_start + 1
        if read_center >= len(ref_positions):
            continue
        # store center context
        center_ref = ref_positions[read_center]
        center_seq = seq[center_start:center_end]
        center_qual = qual[center_start:center_end]
        center_ipd_fwd = fwd_ipd_all[center_start:center_end]
        center_ipd_rev = rev_ipd_all[center_start:center_end]
        center_pw_fwd = fwd_pw_all[center_start:center_end]
        center_pw_rev = rev_pw_all[center_start:center_end]
        # store whole windows
        window_seq = list(seq[win_start:win_end])
        window_qual = list(qual[win_start:win_end])
        window_ipd_fwd = fwd_ipd_all[win_start:win_end]
        window_ipd_rev = rev_ipd_all[win_start:win_end]
        window_pw_fwd = fwd_pw_all[win_start:win_end]
        window_pw_rev = rev_pw_all[win_start:win_end]
        window_sx = sx_all[win_start:win_end]
        window_sm = sm_all[win_start:win_end]
        if len(center_seq) != 3 or len(center_ipd_fwd) != 3 or len(window_ipd_fwd) != window_size:
            continue
        row = {
            "contig": read.reference_name,
            "read_name": read.query_name,
            "fn": fn,
            "rn": rn,
            "ref_center": center_ref,
            "read_center": read_center,
            "center_seq": center_seq,
            "center_qual": center_qual,
            "center_ipd_fwd_0": center_ipd_fwd[0],
            "center_ipd_fwd_1": center_ipd_fwd[1],
            "center_ipd_fwd_2": center_ipd_fwd[2],
            "center_ipd_rev_0": center_ipd_rev[0],
            "center_ipd_rev_1": center_ipd_rev[1],
            "center_ipd_rev_2": center_ipd_rev[2],
            "center_pw_fwd_0": center_pw_fwd[0],
            "center_pw_fwd_1": center_pw_fwd[1],
            "center_pw_fwd_2": center_pw_fwd[2],
            "center_pw_rev_0": center_pw_rev[0],
            "center_pw_rev_1": center_pw_rev[1],
            "center_pw_rev_2": center_pw_rev[2],
            "window_sm": window_sm,
            "window_sx": window_sx,
            "window_seq": window_seq,
            "window_qual": window_qual,
            "window_ipd_fwd": window_ipd_fwd,
            "window_ipd_rev": window_ipd_rev,
            "window_pw_fwd": window_pw_fwd,
            "window_pw_rev": window_pw_rev
        }
        rows.append(row)
    return rows

def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    bam_filepath = os.path.expanduser("~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics/ob006_kinetics_diploid.bam")
    bed_filepath = "ob006_true_mutations_padded_complement.bed"  # BED file with contig, start, end per region
    test_mode = True  # Set to True for testing (process fewer regions)
    all_rows = []
    regions = load_bed_file(bed_filepath)
    with pysam.AlignmentFile(bam_filepath, "rb") as bam:
        for contig, start, end in regions:
            for read in tqdm(bam.fetch(contig, start, end), desc=f"{contig}:{start}-{end}"):
                if read.is_unmapped:
                    continue
                logging.info(f"Processing read {read.query_name} aligned at {read.reference_start} in region {contig}:{start}-{end}")
                rows = process_read_for_null(read, window_size=64)
                all_rows.extend(rows)

            if test_mode:
                break
    if all_rows:
        df = pl.DataFrame(all_rows)
        output_parquet = "./output/null_distribution_test4.parquet"
        os.makedirs(os.path.dirname(output_parquet), exist_ok=True)
        df.write_parquet(output_parquet)
        logging.info(f"Output written to {output_parquet}")
    else:
        logging.info("No null data generated.")

if __name__ == "__main__":
    main()
