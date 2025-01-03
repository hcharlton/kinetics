#!/usr/bin/env python
import pysam  # type: ignore
import sys
import os


def get_ipd_data(bam_file, contig, start, end, read_name, position_in_read):
    kinetics = pysam.AlignmentFile(bam_file, 'rb')
    for read in kinetics.fetch(contig=contig, start=start, end=end):
        if read.query_name == read_name:
            if read.is_reverse:
                ipd_fwd = list(read.get_tag('ri')[position_in_read-10:position_in_read+11])
                ipd_rev = list(read.get_tag('fi')[position_in_read-10:position_in_read+11])
            else:
                ipd_fwd = list(read.get_tag('fi')[position_in_read-10:position_in_read+11])
                ipd_rev = list(read.get_tag('ri')[position_in_read-10:position_in_read+11])
            return ipd_fwd, ipd_rev
    return None, None

def main(bed_line, bam_file):
    header_names = ['unique_id',
                    'contig',
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
    row = dict(zip(header_names, bed_line))
    # print(row)
    contig = row['contig']
    start = int(row['start'])
    end = int(row['end'])
    read_name = row['read_name']
    unique_id = row['unique_id']
    position_in_read = int(row['position_in_read'])
    
    ipd_fwd, ipd_rev = get_ipd_data(bam_file, contig, start, end, read_name, position_in_read)
    if ipd_fwd and ipd_rev:
        print(f"{unique_id}\t{ipd_fwd}\t{ipd_rev}")
    else:
        print("NULL")
        # with open(output_file, 'w') as f:
        #     f.write(f"{unique_id}\t{ipd_fwd}\t{ipd_rev}")



# if __name__ == "__main__":
#     bed_line = sys.argv[1]
#     # bam_file = sys.argv[2]
#     bam_file = os.path.expanduser('~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics/ob006_kinetics_diploid.bam')
#     main(bed_line, bam_file)

if __name__ == "__main__":
    import ast
    bed_line = ast.literal_eval(sys.argv[1])  # Convert string to list
    bam_file = os.path.expanduser('~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics/ob006_kinetics_diploid.bam')
    main(bed_line, bam_file)