import pysam

human_id = 'ob2v'
position_id_ob2v = ['220627', '220704', '220705', '220703', '220718', '220719']
reference = pysam.FastaFile('../../../DerivedData/fasta/HiFi/human/{human_id}/haploid_assembly/{human_id}_haploid.fa'.format(human_id = human_id))
haploid_bamfile = pysam.AlignmentFile('../../../DerivedData/bam/HiFi/human/{human_id}/kinetics/{human_id}_kinetics_haploid.bam'.format(human_id = human_id), "rb")
#diploid_bamfile = pysam.AlignmentFile('../bam/human/{human_id}/pacbio/kinetics/{human_id}_kinetics_diploid.bam'.format(human_id = human_id), "rb")

haploid_contig_name = 'ptg000001l'
context = 'TCT'

def reverse_complement(seq):
    reverse_complement = ''
    watson_crick_pairs = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    reverse = seq[::-1]
    for base in reverse:
            reverse_complement += watson_crick_pairs[base]
    return reverse_complement

def extract_read_names_and_coordinates(haploid_contig_name, context):
    mutation_coordinates = []
    mutation_read_names = []
    diploid_contig = []  
    with open('../../pepo/de_novo_assembly/results/variant_calling/human/ob2v/pacbio/haploid/post_blast_ob2v_true_mutations.txt', 'r') as df:
        for i, line in enumerate(df):
            if line.split('\t')[5] == haploid_contig_name and line.split('\t')[12] == context:# line.split('\t')[6] == 'C' and line.split('\t')[7] == 'T':
                diploid_contig.append(line.split('\t')[0])
                mutation_read_names.append(line.split('\t')[4])
                mutation_coordinates.append(int(line.split('\t')[7]))
        return mutation_coordinates, mutation_read_names, diploid_contig

def mutation_kinetics(number_of_flanking_bases):
    i = 0
    mutation_coordinates = extract_read_names_and_coordinates(haploid_contig_name, context)[0]
    mutation_read_names = extract_read_names_and_coordinates(haploid_contig_name, context)[1]
    diploid_contig = extract_read_names_and_coordinates(haploid_contig_name, context)[2]
    for coordinate in mutation_coordinates:
        for pileupcolumn in haploid_bamfile.pileup(haploid_contig_name, coordinate-1, coordinate, truncate=True):
            if '+' in mutation_read_names[i]:
                i += 1
                continue
            read_names = [pileupread.alignment.query_name for pileupread in pileupcolumn.pileups if pileupread.query_position != None]
            strand = [pileupread.alignment.is_reverse for pileupread in pileupcolumn.pileups if pileupread.query_position != None]
            if strand[read_names.index(mutation_read_names[i])] == True: # In this line the statement is 'True' if the read is from the negative strand
                ipd = [pileupread.alignment.get_tag('ri') for pileupread in pileupcolumn.pileups if pileupread.query_position != None]
                strand = '-'
            else:
                ipd = [pileupread.alignment.get_tag('fi') for pileupread in pileupcolumn.pileups if pileupread.query_position != None]
                strand = '+'
            read_start_position = [pileupread.alignment.reference_start + 1 for pileupread in pileupcolumn.pileups if pileupread.query_position != None]
            print(pileupcolumn.pos+1, 
                  #(*list(ipd[read_names.index(mutation_read_names[i])][pileupcolumn.get_query_positions()[read_names.index(mutation_read_names[i])]-number_of_flanking_bases:pileupcolumn.get_query_positions()[read_names.index(mutation_read_names[i])]+number_of_flanking_bases+1])), 
                  strand, 
                  #mutation_read_names[i], 
                  #diploid_contig[i],
                  'mutation_kinetics', 
                  sep = '\t')
            i += 1
    return

def find_next_context(haploid_contig_name, context, coordinate):
    i = 0
    #mutation_coordinates = extract_read_names_and_coordinates(haploid_contig_name, context)[0]
    #for coordinate in mutation_coordinates:
    fasta_string = reference.fetch(start = coordinate-1, region = haploid_contig_name)
    while True:
        trinucleotide = fasta_string[i:3+i]
        if trinucleotide == context or trinucleotide == reverse_complement(context):
            first_positon_in_context = coordinate+i
            return first_positon_in_context
        else:
            i += 1
            trinucleotide == fasta_string[i:3+i]

def reference_kinetics(coordinate):
    count = 1
    first_coordinate_of_next_context = find_next_context(haploid_contig_name, context, coordinate)
    last_coordinate_of_next_context = first_coordinate_of_next_context + 3
    first_column = []
    second_column = []
    third_column = []
    for pileupcolumn in haploid_bamfile.pileup(haploid_contig_name, start = first_coordinate_of_next_context - 1, end = last_coordinate_of_next_context, truncate=True):
        ipd = []
        read_names = [pileupread.alignment.query_name for pileupread in pileupcolumn.pileups]
        for pileupread in pileupcolumn.pileups:
            if pileupread.alignment.is_reverse == True and pileupread.query_position != None: # In this line the statement is 'True' if the read is from the negative strand
                ipd.append((pileupread.alignment.query_name, pileupread.alignment.get_tag('ri')))
            if pileupread.alignment.is_reverse != True and pileupread.query_position != None:
               ipd.append((pileupread.alignment.query_name, pileupread.alignment.get_tag('fi')))
        if len(ipd) != len(read_names):
            print('NA', 'NA', 'NA', 'NA', 'reference_kinetics' ,sep = '\t')
            return
        for i in range(len(ipd)): # The coverage is not alway equal to the number of reads at the postion since some bases have a quality score below the threshold
            if count == 1 and len(ipd[i][1]) != 0: # The second condition is added, because some reads do not have kinetics data
                first_column.append(ipd[i][1][pileupcolumn.get_query_positions()[i]]) # To check if the kinetic value correspond to the read name make a tuple: first_column.append((ipd[i][0], ipd[i][1][pileupcolumn.get_query_positions()[i]]))  ¨
            if count == 2 and len(ipd[i][1]) != 0:
                #print(ipd[i][1])
                #print(pileupcolumn.pos+1)
                #print(pileupcolumn.get_query_positions()[i])
                #print(len(ipd[i][1]))
                #print(ipd[i][1][pileupcolumn.get_query_positions()[i]])
                second_column.append(ipd[i][1][pileupcolumn.get_query_positions()[i]])
                second_column_pos = pileupcolumn.pos+1
            if count == 3 and len(ipd[i][1]) != 0:
                third_column.append(ipd[i][1][pileupcolumn.get_query_positions()[i]])
        count += 1
        if count == 4:
            print(second_column_pos, sum(first_column) / len(second_column), sum(second_column) / len(second_column), sum(third_column) / len(third_column), 'reference_kinetics', sep = '\t')


mutation_kinetics(6)
for coordinate in extract_read_names_and_coordinates(haploid_contig_name, context)[0]:
    output = reference_kinetics(coordinate)


#print(extract_read_names_and_coordinates(haploid_contig_name, context)[0])
#reference_kinetics(94391394)
