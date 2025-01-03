from gwf import Workflow, AnonymousTarget
import os

gwf = Workflow(defaults={'cores': 1, 'queue':"normal",'memory': "16g", 'walltime':"00:00:60", 'account':"mutationalscanning"})


def move(inputfile, outputfile):
    """A template for moving files."""
    inputs = [inputfile]
    outputs = [outputfile]
    options = {
        'cores': 1,
        'memory': '2g',
    }

    spec = '''
    cat {} > {}
    '''.format(inputfile, outputfile)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


gwf.target_from_template(
    name='test1',
    template=move(
        inputfile='input.txt',
        outputfile = 'other_input.txt'
    )
)


