from gwf import Workflow, AnonymousTarget  # type: ignore
import os
import json 

gwf = Workflow(defaults={'cores': 1, 'memory': '2g', 'walltime':"00:00:30", 'queue':"normal", 'account':"mutationalscanning"})

# SHORT BED FILE FOR TESTING
bed_file = './ob006_10.bed'

#FULL BAM FILE
bam_file = os.path.expanduser('~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics/ob006_kinetics_diploid.bam')

content = []
with open(bed_file)as f:
    for line in f:
        content.append(line.strip().split())

bed_line = content[1]


def get_ipd(input_bed_line):
    """a template to take a part of the bed line and copy it to a file"""
    inputs = []
    output_dir = 'output'
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{input_bed_line[0]}.o')
    outputs = [output_file]
    options = {}
    input_bed_line_json = json.dumps(input_bed_line) 
    spec = f'''
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate kinetics
        python bin/retrieve_ipd.py '{input_bed_line_json}' > {output_file}'''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def combine_kinetics(paths, output_path):
    inputs = {'paths': paths}
    outputs = {'zipped_file': output_path}
    options = {}
    spec = """zip ..."""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


for index, bed_line in enumerate(content):
    gwf.target_from_template(f'ipd_job.{index}', get_ipd(bed_line))



