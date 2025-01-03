def move_bed(input_bed_line):
    """a template to take a part of the bed line and copy it to a file"""
    inputs = []
    output_dir = 'output'
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{input_bed_line[:36]}.o')
    outputs = [output_file]
    options = {}
    spec = '''
        echo {} > {}'''.format(input_bed_line,output_file)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# gwf.target_from_template(
#     name='move_a_line',
#     template=move_bed(
#         input_bed_line=' '.join(bed_line)
#     )
# )

# header_names = ['unique_id',
#                     'contig',
#                     'start',
#                     'end',
#                     'ref',
#                     'alt',
#                     'context',
#                     'oriented_ref',
#                     'oriented_alt',
#                     'oriented_context',
#                     'mutation_qual',
#                     'strand',
#                     'reads_with_alt',
#                     'coverage',
#                     'sm',
#                     'sx',
#                     'read_name',
#                     'avg_read_lengths',
#                     'position_in_read']

# print(bed_line)


# bed_file = '../../data/ob006_true_mutations_noheaders.bed'

# SHORT BED FILE FOR TESTING
bed_file = './ob006_10.bed'

#FULL BAM FILE
bam_file = '~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics/ob006_kinetics_diploid.bam'

with open(bed_file) as f:
    bed_lines = f.readlines()
# template for the ipd retrieval script
def get_ipd(bed_line):
    inputs = [bed_line]
    print(bed_line)
    output_dir = 'ipd_output'
    os.makedirs(output_dir, exist_ok=True)
    outfile = os.path.join(output_dir, f"{os.path.basename(bam_file)}_{bed_line[:36]}.ipd")
    print(outfile)
    outputs = [outfile]
    options = {'cores': 1, 'memory': "8g", 'walltime':"00:00:30"}
    spec = '''
            source $(conda info --base)/etc/profile.d/conda.sh
            mamba activate kinetics_gwf
            echo '{bed_line}' > out_file
            '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def copy_bed_line(bed_line):
    inputs = bed_line
    output_dir = 'test_output'
    out_file = os.path.join(output_dir, 'test.o')
    outputs = [out_file]
    options = {'cores': 1, 'memory': "8g", 'walltime':"00:00:30"}
    spec = '''echo {bed_line} > out_file'''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

gwf.map(copy_bed_line, bed_lines)
# gwf.target_from_template(
#     name='test4',
#     template=spl526135it_bed_file(
#         input_file=bed_file,
#         output_dir='split_bed_files/'
#     )
# )
# gwf.target_from_template(
#     name='test1',
#     template=get_ipd(bed_lines[0])
# )

#             ./bin/retrieve_ipd.py {bed_line} {bam_file} > {outfile}

# print(bed_lines[0])



bed_line = content[0]
print(len(content[0]))
def copy_bed_line(bed_line):
    inputs = {'bed_line': bed_line}
    print(bed_line)
    output_dir = 'test_output'
    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, f'{bed_line[0]}.o')
    outputs = [out_file]
    options = {'cores': 1, 'memory': "8g", 'walltime':"00:00:30"}
    spec = '''echo {bed_line[4]} > out_file'''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

gwf.map(copy_bed_line, content)