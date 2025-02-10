from gwf import Workflow, AnonymousTarget  # type: ignore
import os
import json 
from pathlib import Path


gwf = Workflow(defaults={'cores': 1, 'memory': '2g', 'walltime':"00:00:30", 'queue':"normal", 'account':"mutationalscanning"})

# SHORT BED FILE FOR TESTING
# bed_file = './ob006_10.bed'

# Full BED file
bed_file = '../../data/ob006_true_mutations_noheaders.bed'

#FULL BAM FILE
bam_file = os.path.expanduser('~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics/ob006_kinetics_diploid.bam')

content = []
with open(bed_file)as f:
    for line in f:
        content.append(line.strip().split())


def get_ipd(input_bed_line):
    """A template to retrieve kinetics from a bam file given a bed_line"""
    inputs = []
    output_dir = 'output'
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{input_bed_line[0]}.json')
    outputs = [output_file]
    options = {}
    input_bed_line_json = json.dumps(input_bed_line) 
    spec = f'''
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate kinetics
        python bin/retrieve_ipd.py '{input_bed_line_json}' > {output_file}'''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def make_ndjson(json_dir):
    """A template for combining a directory of json files into an ndjson"""
    inputs = [json_dir]
    ndjson_file = os.path.join(json_dir, f'{Path(bed_file).stem}.ndjson')
    outputs = [ndjson_file]
    options = {'cores': 4, 'memory': '16g'}
    spec = f'''
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate kinetics
        find {json_dir} -name "*.json" -print0 | xargs -0 cat -s >> {ndjson_file}
        '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# generate json kinetics files for each bed line
for index, bed_line in enumerate(content):
    gwf.target_from_template(f'ipd_job.{index}', get_ipd(bed_line))

# combine them into a ndjson
gwf.target_from_template('ndjson_maker', make_ndjson(json_dir='output'))

        # if [[{cleanup}]]
        #     find {json_dir} -name "*.json" -delete

