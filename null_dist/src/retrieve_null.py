import pysam, logging, os
from collections import defaultdict
from tqdm import tqdm
import polars as pl

def process_read_for_null(read, window_size: int = 64) -> dict[str, list]:
    kinetics_tag = "ri" if read.is_reverse else "fi"
    seq = read.query_sequence
    try:
        ipd_all = list(read.get_tag(kinetics_tag))

    except Exception as e:
        logging.warning(f"Error getting kinetics tag {kinetics_tag} for {read.query_name}: {e}")
        return {}
    if len(seq) < window_size:
        logging.info(f"Read {read.query_name} is shorter than window size ({window_size}).")
        return {}
    null_data = defaultdict(list)
    num_windows = len(seq) // window_size
    for i in range(num_windows):
        win_start = i * window_size
        win_end = win_start + window_size
        if win_end > len(seq):
            continue
        # Use the center 3-nt (positions 31-33) as a grouping label.
        center_start = win_start + 31
        center_end = center_start + 3
        if center_end > win_end:
            continue
        center_context = seq[center_start:center_end]
        window_seq = seq[win_start:win_end]
        window_ipd = ipd_all[win_start:win_end]
        if len(window_seq) == window_size and len(window_ipd) == window_size:
            null_data[center_context].append((window_seq, window_ipd))
    return null_data

def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    bam_filepath = "~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics/ob006_kinetics_diploid.bam"    
    test_contig = "h1tg000020l"         
    test_mode = True             
    overall_null = defaultdict(list)
    with pysam.AlignmentFile(os.path.expanduser(bam_filepath), "rb") as bam:
        for read in tqdm(bam.fetch(test_contig)):
            if read.is_unmapped:
                continue
            data = process_read_for_null(read, window_size=64)
            for context, window_data in data.items():
                overall_null[context].extend(window_data)
            if test_mode:
                break
    # Flatten the null distribution into rows for file output.
    rows = []
    # Each row will include the center context, the full window sequence, and the kinetics values for that window.
    for context, window_list in overall_null.items():
        for (window_seq, window_ipd) in window_list:
            rows.append({
                "center_context": context,
                "window_seq": window_seq,
                "ipd_window": window_ipd
            })
    if rows:
        df = pl.DataFrame(rows)
        output_parquet = "./output/null_distribution_full_window_test.parquet"
        df.write_parquet(output_parquet)
        logging.info(f"Output written to {output_parquet}")
    else:
        logging.info("No null data generated.")

if __name__ == "__main__":
    main()
