#!/usr/bin/env nextflow


params.str = 'Hello world!'

process splitLetters {
    conda 'pandas seaborn'
    output:
    path 'chunk_*'

    """
    printf '${params.str}' | split -b 6 - chunk_
    """
}

process convertToUpper {
    input:
    path x

    output:
    stdout

    """
    rev $x
    """
}

workflow {
    splitLetters | flatten | convertToUpper | view { v -> v.trim() }
}