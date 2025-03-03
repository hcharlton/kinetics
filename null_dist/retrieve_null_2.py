import pysam, logging, os
from tqdm import tqdm
from collections import defaultdict
import polars as pl

def process_read_for_null(read, window_size: int = 64) -> list[dict]:
    # kinetics_tag = "ri" if read.is_reverse else "fi"
    fwd_tag, rev_tag = (("ri", "fi") if read.is_reverse else ("fi", "ri"))
    seq = read.query_sequence
    try:
        fwd_ipd_all = list(read.get_tag(fwd_tag))
        rev_ipd_all = list(read.get_tag(rev_tag))
        rev_ipd_all.reverse()
        print(len(fwd_ipd_all),len(rev_ipd_all))
    except Exception as e:
        logging.warning(f"Error getting kinetics tags for {read.query_name}: {e}")
        return []
    if len(seq) < window_size:
        logging.info(f"Read {read.query_name} is shorter than window size ({window_size}).")
        return []
    # Get reference positions for each base in the read.
    # Note: This assumes a simple alignment (e.g. no indels affecting position mapping).
    ref_positions = read.get_reference_positions()
    rows = []
    num_windows = len(seq) // window_size
    for i in range(num_windows):
        win_start = i * window_size
        win_end = win_start + window_size
        if win_end > len(seq):
            continue
        # Define the center of the window:
        # We'll use the three bases starting at index win_start+31 (i.e. indices 31,32,33)
        center_start = win_start + 31
        center_end = center_start + 3
        # The read center is taken as the middle of the center triplet.
        read_center = center_start + 1
        if read_center >= len(ref_positions):
            continue
        ref_center = ref_positions[read_center]
        center_seq = seq[center_start:center_end]
        center_ipd_fwd = fwd_ipd_all[center_start:center_end]
        center_ipd_rev = rev_ipd_all[center_start:center_end]
        window_seq = list(seq[win_start:win_end])
        window_ipd_fwd = fwd_ipd_all[win_start:win_end]
        window_ipd_rev = rev_ipd_all[win_start:win_end]
        if len(center_seq) != 3 or len(center_ipd_fwd) != 3 or len(window_ipd_fwd) != window_size:
            continue
        row = {
            "contig": read.reference_name,
            "read_name": read.query_name,
            "ref_center": ref_center,
            "read_center": read_center,
            "center_seq": center_seq,
            "center_ipd_fwd_0": center_ipd_fwd[0],
            "center_ipd_fwd_1": center_ipd_fwd[1],
            "center_ipd_fwd_2": center_ipd_fwd[2],
            "center_ipd_rev_0": center_ipd_rev[0],
            "center_ipd_rev_1": center_ipd_rev[1],
            "center_ipd_rev_2": center_ipd_rev[2],
            "window_seq": window_seq,
            "window_ipd_fwd": window_ipd_fwd,
            "window_ipd_rev": window_ipd_rev
        }
        rows.append(row)
    return rows

def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    bam_filepath = "~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics/ob006_kinetics_diploid.bam"
    test_contig = "h1tg000002l"  # for now this is just one contig
    test_mode = True           # Set to False to process the whole file.
    all_rows = []
    with pysam.AlignmentFile(os.path.expanduser(bam_filepath), "rb") as bam:
        for read in tqdm(bam.fetch(test_contig)):
            if read.is_unmapped:
                continue
            logging.info(f"Read {read.query_name} aligned at {read.reference_start}")
            rows = process_read_for_null(read, window_size=64)
            all_rows.extend(rows)
            if test_mode:
                break
    if all_rows:
        df = pl.DataFrame(all_rows)
        output_parquet = "./output/null_distribution_test2.parquet"
        df.write_parquet(output_parquet)
        logging.info(f"Output written to {output_parquet}")
    else:
        logging.info("No null data generated.")

if __name__ == "__main__":
    main()
