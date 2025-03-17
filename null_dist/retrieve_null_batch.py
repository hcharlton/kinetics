import os
from typing import Optional, List
import pysam
import polars as pl
from tqdm import tqdm
from line_profiler import profile

WINDOW_SIZE: int = 64
CENTER_OFFSET: int = (WINDOW_SIZE // 2) - 1
CENTER_SIZE: int = 3
START_WINDOW: int = 2
BATCH_SIZE: int = 100000

SCHEMA = {
    "contig": pl.Utf8,
    "read_name": pl.Utf8,
    "read_is_reverse": pl.Boolean,
    "fn": pl.Int64,
    "rn": pl.Int64,
    "ref_center": pl.Int64,
    "read_center": pl.Int64,
    "center_seq": pl.Utf8,
    "center_qual": pl.List(pl.UInt8),
    "center_ipd_fwd_0": pl.UInt16,
    "center_ipd_fwd_1": pl.UInt16,
    "center_ipd_fwd_2": pl.UInt16,
    "center_ipd_rev_0": pl.UInt16,
    "center_ipd_rev_1": pl.UInt16,
    "center_ipd_rev_2": pl.UInt16,
    "center_pw_fwd_0": pl.UInt8,
    "center_pw_fwd_1": pl.UInt8,
    "center_pw_fwd_2": pl.UInt8,
    "center_pw_rev_0": pl.UInt8,
    "center_pw_rev_1": pl.UInt8,
    "center_pw_rev_2": pl.UInt8,
    "window_sm": pl.List(pl.UInt16),
    "window_sx": pl.List(pl.UInt16),
    "window_seq": pl.List(pl.Utf8),
    "window_qual": pl.List(pl.UInt8),
    "window_ipd_fwd": pl.List(pl.UInt16),
    "window_ipd_rev": pl.List(pl.UInt16),
    "window_pw_fwd": pl.List(pl.UInt8),
    "window_pw_rev": pl.List(pl.UInt8),
}

@profile
def process_read(read: pysam.AlignedSegment, window_size: int = WINDOW_SIZE) -> Optional[pl.DataFrame]:
    try:
        tags = {
            "fn": read.get_tag("fn"),
            "rn": read.get_tag("rn"),
            "sm": read.get_tag("sm"),
            "sx": read.get_tag("sx"),
            "qual": read.query_qualities,
            "seq": read.query_sequence,
        }
        if read.is_reverse:
            fwd_ipd_tag, rev_ipd_tag = "ri", "fi"
            fwd_pw_tag, rev_pw_tag = "rp", "fp"
        else:
            fwd_ipd_tag, rev_ipd_tag = "fi", "ri"
            fwd_pw_tag, rev_pw_tag = "fp", "rp"
        for tag in [fwd_ipd_tag, rev_ipd_tag, fwd_pw_tag, rev_pw_tag]:
            tags[tag] = read.get_tag(tag)
        if read.is_reverse:
            tags[rev_ipd_tag] = tags[rev_ipd_tag][::-1]
            tags[rev_pw_tag] = tags[rev_pw_tag][::-1]
    except Exception as e:
        print(f"Error getting tags for {read.query_name}: {e}")
        return None

    expected_len = len(tags["seq"])
    required_tags = ["qual", "sm", "sx", fwd_ipd_tag, rev_ipd_tag, fwd_pw_tag, rev_pw_tag]
    if any(len(tags[tag]) != expected_len for tag in required_tags):
        # print(f"Inconsistent tag lengths in read {read.query_name}")
        return None
    if expected_len < window_size:
        return None

    ref_positions = read.get_reference_positions(full_length=True)
    num_windows = expected_len // window_size
    rows: List[dict] = []

    tags_fwd_ipd_tag = tags[fwd_ipd_tag]
    tags_rev_ipd_tag = tags[rev_ipd_tag]
    tags_fwd_pw_tag = tags[fwd_pw_tag]
    tags_rev_pw_tag = tags[rev_pw_tag]

    for i in range(START_WINDOW, num_windows - START_WINDOW):
        win_start = i * window_size
        win_end = win_start + window_size
        center_start = win_start + CENTER_OFFSET
        center_end = center_start + CENTER_SIZE
        read_center = center_start + (CENTER_SIZE // 2)

        #  or any(pos is None for pos in ref_positions[win_start:win_end])
        if read_center >= len(ref_positions):
            continue

        center_seq = tags["seq"][center_start:center_end]
        if len(center_seq) != CENTER_SIZE:
            continue

        window_seq = list(tags["seq"][win_start:win_end])
        if len(window_seq) != window_size:
            continue

        row = {
            "contig": read.reference_name,
            "read_name": read.query_name,
            "read_is_reverse": read.is_reverse,
            "fn": tags["fn"],
            "rn": tags["rn"],
            "ref_center": ref_positions[read_center],
            "read_center": read_center,
            "center_seq": center_seq,
            "center_qual": tags["qual"][center_start:center_end],
            "center_ipd_fwd_0": tags_fwd_ipd_tag[center_start],
            "center_ipd_fwd_1": tags_fwd_ipd_tag[center_start + 1],
            "center_ipd_fwd_2": tags_fwd_ipd_tag[center_start + 2],
            "center_ipd_rev_0": tags_rev_ipd_tag[center_start],
            "center_ipd_rev_1": tags_rev_ipd_tag[center_start + 1],
            "center_ipd_rev_2": tags_rev_ipd_tag[center_start + 2],
            "center_pw_fwd_0": tags_fwd_pw_tag[center_start],
            "center_pw_fwd_1": tags_fwd_pw_tag[center_start + 1],
            "center_pw_fwd_2": tags_fwd_pw_tag[center_start + 2],
            "center_pw_rev_0": tags_rev_pw_tag[center_start],
            "center_pw_rev_1": tags_rev_pw_tag[center_start + 1],
            "center_pw_rev_2": tags_rev_pw_tag[center_start + 2],
            "window_sm": tags["sm"][win_start:win_end],
            "window_sx": tags["sx"][win_start:win_end],
            "window_seq": window_seq,
            "window_qual": tags["qual"][win_start:win_end],
            "window_ipd_fwd": tags_fwd_ipd_tag[win_start:win_end],
            "window_ipd_rev": tags_rev_ipd_tag[win_start:win_end],
            "window_pw_fwd": tags_fwd_pw_tag[win_start:win_end],
            "window_pw_rev": tags_rev_pw_tag[win_start:win_end],
        }
        rows.append(row) # comment out for testing
    return rows
    #return pl.DataFrame(rows, schema=SCHEMA)

@profile
def write_batch(rows, batch_counter: int, output_dir: str) -> None:
    df = pl.DataFrame(rows, schema=SCHEMA)
    output_parquet = os.path.join(output_dir, f"batch_{batch_counter:05d}.parquet")
    df.write_parquet(output_parquet) 
    print(f"Batch {batch_counter} written to {output_parquet}")

def main() -> None:
    bam_filepath = os.path.expanduser(
        "~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics/ob006_kinetics_diploid.bam"
    )
    test_contig = "h1tg000002l"
    output_dir = "./batch_test"
    os.makedirs(output_dir, exist_ok=True)

    batch_counter = 0
    read_counter = 0
    #current_batch: List[pl.DataFrame] = []
    current_batch = []
    with pysam.AlignmentFile(bam_filepath, "rb") as bam:
        for read in tqdm(bam.fetch(test_contig), desc="Processing Reads", unit="read"):
            if read.is_unmapped:
                continue
            rows = process_read(read)
            current_batch.extend(rows or [])
            #current_batch.append(df)
            read_counter += 1
            if read_counter >= BATCH_SIZE:
                write_batch(current_batch, batch_counter, output_dir)
                batch_counter += 1
                read_counter = 0
                current_batch = []
    if current_batch:
        write_batch(current_batch, batch_counter, output_dir)

if __name__ == "__main__":
    main()
