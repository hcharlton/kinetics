import pysam
import pytest
import polars as pl
from retrieve_null_batch import process_read, SCHEMA

@pytest.fixture
def create_test_bam(tmp_path: str) -> str:
    bam_path = tmp_path / "test.bam"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        aln = pysam.AlignedSegment()
        aln.query_name = "read_1"
        aln.query_sequence = "ACTG" * 20
        aln.flag = 0
        aln.reference_id = 0
        aln.reference_start = 100
        aln.mapping_quality = 60
        aln.cigar = [(0, 80)]
        aln.query_qualities = pysam.qualitystring_to_array("I" * 80)
        aln.set_tag("fn", 1)
        aln.set_tag("rn", 2)
        aln.set_tag("sm", [10] * 80)
        aln.set_tag("sx", [20] * 80)
        aln.set_tag("fi", [1] * 80)
        aln.set_tag("ri", [2] * 80)
        aln.set_tag("fp", [3] * 80)
        aln.set_tag("rp", [4] * 80)
        bam.write(aln)
    pysam.index(str(bam_path))
    return str(bam_path)

def test_process_read(create_test_bam: str):
    with pysam.AlignmentFile(create_test_bam, "rb") as bam:
        reads = list(bam.fetch("chr1"))
        assert len(reads) == 1
        rows = process_read(reads[0])
        assert rows is not None
        df = pl.DataFrame(rows, schema=SCHEMA)

        assert df.shape[0] > 0
        assert df["read_name"][0] == "read_1"
        assert df["contig"][0] == "chr1"
        assert df["fn"][0] == 1
        assert df["rn"][0] == 2
        assert df["center_seq"][0] == "TGA"
        assert all(v == 10 for v in df["window_sm"][0])
        assert all(v == 20 for v in df["window_sx"][0])
        assert all(v == 1 for v in df["window_ipd_fwd"][0])
        assert all(v == 2 for v in df["window_ipd_rev"][0])
