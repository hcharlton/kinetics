


## make a genome file for the lengths of each contig
This command extracts the lenghts and names of each of the contigs and stores it in a txt file

### recipe
`samtools idxstats yourfile.bam | cut -f1,2 > genome.txt`

### actual (verified)
`samtools idxstats ~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics/ob006_kinetics_diploid.bam | cut -f1,2 > ob006_genome_index.txt`

## pad the observed muatations BED file so we ignore a 16 nt context
This takes a BED file and pads the intervals. Since the BED file with the observed mutations only lists the region as the single nucleotide location, I want to pad the region. Previous work demonstrated that around 16 nucleotides away from the the observed mutation was enough to return to background dynamics of the kinetics. 
### recipe
`bedtools slop -i input.bed -g ob006_genome_index.txt -b 16 > padded.bed`
### test
`bedtools slop -i ~/mutationalscanning/Workspaces/chcharlton/kinetics/data/ob006_10_noheader.bed -g ob006_genome_index.txt -b 16 > ob006_10_noheader_padded.bed`
### actual

## generate a complement of the padded BED

### recipe
`bedtools complement -i padded_merged.bed -g genome.txt > complement.bed`
### test
`bedtools complement -i ob006_10_noheader_padded.bed -g ob006_genome_index.txt > ob006_test_complement.bed`
### actual

## gather context for kinetics 
This includes ipd_fwd, ipd_rev, pw_fwd, pw_rev, fn, rn, 



`python retrieve_null.py ob006_test_complement_head.bed ~/mutationalscanning/DerivedData/bam/HiFi/human/ob006/kinetics/ob006_kinetics_diploid.bam ~/mutationalscanning/DerivedData/fasta/HiFi/human/ob006/diploid_assembly/ob006_diploid.fa ./null_dist.ndjson`









m84108_240530_140750_s2/73663713/ccs
m84108_240530_160720_s3/14616969/ccs









