#!/usr/bin/env nextflow

def bed_file = file(params.bed_filepath)

process RetrieveIPD {
    conda 'pysam'
    input:
        val index
    output:
        path "${index}.o", emit: files
    script:
    """
    bed_line=\$(sed '${index}q;d' ${bed_file})
    retrieve_ipd.py "\$bed_line" ${params.bam_filepath} > ${index}.o
    """
}
process CombineOutputs {
    input:
        path files
    output:
        path 'combined_output.csv'
    publishDir "output_files", mode: 'copy'
    script:
    def cat_files = files.collect{ "cat $it >> combined_output.csv" }.join("\n")
    """
    ${cat_files}
    """
}

workflow {
    // def index = channel.of(1..bed_file.readLines().size())
    def index = channel.of(1..500)
    RetrieveIPD(index)
    CombineOutputs(RetrieveIPD.out.files.collect())
}





















// // process AggregateResults {
// //     input:
// //     path ipd_results

// //     output:
// //     path 'final_results.txt'

// //     script:
// //     """
// //     cat output/* > final_results.txt
// //     """
// // }