#!/usr/bin/env python3

import argparse, sys, os, glob, numpy
from collections import Counter, OrderedDict

def main():

    #
    # parser object and command line input arguments
    #
    parser = argparse.ArgumentParser(description="Search open reading frames (ORFs) in a sequence.")
    parser.add_argument("input", help="Enter file/s to analyze (for path to\
                        directory enable batch mode option.)", nargs='+')
    parser.add_argument("-b", "--batch", help="Enables batch mode for processing\
                        all '.fasta' files in a directory. Path to directory must\
                        be specified as input.", action="store_true")
    parser.add_argument("-o", "--output", help="Allows output to specified file.")
    parser.add_argument("-f", "--frames",  nargs='+', default=['+1','+2','+3',
                         '-1','-2','-3'], choices=['+1','+2','+3','-1','-2','-3'],
                        help="Allows a customized selection of reading frames for\
                        ORF search. Default is all reading frames.")
    parser.add_argument("-ml", "--minimum_length", type=int, default=300, help=
                        "Allows a customized minimum length threshold for ORF \
                        detection. Default is 300 nucleotides.")
    parser.add_argument("-a", "--ambiguities", help="Enables the detection of ORFs\
                         containing ambiguous nucleotides.", action="store_true")
    parser.add_argument("-fov", "--filter_overlaps", choices=['length', 'cai'],
                        help="Allows filtering overlapping ORFs based on the \
                        selected criteria ('length' or 'cai').")
    parser.add_argument("-t", "--type", default='peptide', choices=['peptide',
                        'cdna'], help="Allows a custom ORF output type. Default\
                         is in peptide sequences.")
    parser.add_argument("-fmt", "--format", type=int, default=70, help="Allows \
                        a customized line length for the ORF output. Default is \
                        70 characters per line.")
    args = parser.parse_args()

    # unique frames
    args.frames = list(set(args.frames))

    # check minimum length value
    if args.minimum_length <= 0:
        sys.exit("ERROR: Invalid minimum ORF length! Please, enter a positive integer value")

    # check output line length value
    if args.format <= 0:
        sys.exit("ERROR: Invalid line length for output! Please, enter a positive integer value")
    
    # check for already existing output file
    if args.output:
        if os.path.exists(args.output):
            sys.exit("ERROR: Please, select new output. File '{}' already exists!".format(args.output))

    #
    # read sequence from file
    #
    try:
        data = retrieve_data(args.input, args.batch)
    except FileNotFoundError as err:
        sys.exit("ERROR: File '{}' has not been found (check name or path to file)".format(err.filename))
    except IsADirectoryError as err:
        sys.exit("ERROR: Please, provide full path to file. To enable batch mode use '-b' option.")

    #
    # predict ORFs in sequence
    #
    all_orfs = []

    for name, seq in data.items():

        # get ORFs and total codon counts
        orfs, codon_counts  = find_ORFs_in_frames(args.frames, seq, args.minimum_length, args.ambiguities)

        # calculate relative codon frequencies within family
        codon_bias_table = calculate_codon_bias(codon_counts)

        # calculate Codon Adaptation Index (CAI)
        orfs = calculate_CAI(orfs, codon_bias_table)

        # filter overlapping ORFs
        if args.filter_overlaps:
            orfs = exclude_overlaps(orfs, args.filter_overlaps)

        # output formatting
        output = format_output(name, seq, orfs, args.type, args.format)
        all_orfs.extend(output)
    
    #
    # Output ORF predictions
    #   
    if args.output:
        with open (args.output, 'w') as F: 
            F.write('\n'.join(all_orfs))
    else:
        print('\n'.join(all_orfs))


def retrieve_data(input_files, batch_flag):

    """Reads in individual files found in list or directory"""

    input_data = {}

    if batch_flag:
        path = ' '.join(input_files)
        for file in glob.glob(os.path.join(path, '*.fasta')):
            seqs_to_add = read_file(file)
            for name in seqs_to_add.keys():
                if name in input_data.keys():
                    sys.exit("ERROR: At least 2 sequences found with the same header '>{}...'".format(name))
            input_data.update(seqs_to_add)
        if len(input_data) == 0:
            sys.exit("ERROR: Check path or file extension! No suitable files found in path '{}'.".format(path))

    else:
        for file in input_files:
            seqs_to_add = read_file(file)
            for name in seqs_to_add.keys():
                if name in input_data.keys():
                    sys.exit("ERROR: At least 2 sequences found with the same header '>{}...'".format(name))
            input_data.update(seqs_to_add)

    return(input_data)


def read_file(filename):

    """Reads in sequences from fasta or multifasta file formats"""

    seqs_in_file = {}
    strings = []
    name = ''

    with open (filename, 'r') as F:
        for line in F:
            if line.startswith('>'):
                if strings:
                    seqs_in_file[name] = ''.join(strings).upper()
                    strings = []
                name = line.rstrip('\n').split()[0][1:]
                if name in seqs_in_file.keys():
                    sys.exit("ERROR: At least 2 sequences found with the same header '>{}...'".format(name))
            elif name:
                strings.append(line.rstrip('\n'))
        
        if strings:
            seqs_in_file[name] = ''.join(strings).upper()
        else:
            sys.exit("ERROR: No sequences found in '{}'. Only files in FASTA format are supported.".format(filename))

    return(seqs_in_file)


def find_ORFs_in_frames(frames, seq, minimum_length, allow_ambiguity):

    """Finds potential ORFs in a given DNA sequence and registers codon counts"""

    start_codons = ('ATG')
    stop_codons = ('TAG', 'TGA', 'TAA')
    amb_nucl = ('N', 'B', 'D', 'H', 'V', 'K', 'Y', 'S', 'W', 'R', 'M')
    codon_count = {
    'AAA' : 0, 'AAC' : 0, 'AAG' : 0, 'AAT' : 0,
    'ACA' : 0, 'ACC' : 0, 'ACG' : 0, 'ACT' : 0,
    'AGA' : 0, 'AGC' : 0, 'AGG' : 0, 'AGT' : 0,
    'ATA' : 0, 'ATC' : 0, 'ATG' : 0, 'ATT' : 0,
    'CAA' : 0, 'CAC' : 0, 'CAG' : 0, 'CAT' : 0,
    'CCA' : 0, 'CCC' : 0, 'CCG' : 0, 'CCT' : 0,
    'CGA' : 0, 'CGC' : 0, 'CGG' : 0, 'CGT' : 0,
    'CTA' : 0, 'CTC' : 0, 'CTG' : 0, 'CTT' : 0,
    'GAA' : 0, 'GAC' : 0, 'GAG' : 0, 'GAT' : 0,
    'GCA' : 0, 'GCC' : 0, 'GCG' : 0, 'GCT' : 0,
    'GGA' : 0, 'GGC' : 0, 'GGG' : 0, 'GGT' : 0,
    'GTA' : 0, 'GTC' : 0, 'GTG' : 0, 'GTT' : 0,
    'TAA' : 0, 'TAC' : 0, 'TAG' : 0, 'TAT' : 0,
    'TCA' : 0, 'TCC' : 0, 'TCG' : 0, 'TCT' : 0,
    'TGA' : 0, 'TGC' : 0, 'TGG' : 0, 'TGT' : 0,
    'TTA' : 0, 'TTC' : 0, 'TTG' : 0, 'TTT' : 0, 
    }

    orf_count = 0
    orf_codon_list = []
    orfs_found = {}
    seq_length = len(seq)

    # search ORFs in different frames
    for frame in frames:
        dna = seq
        start = None
        amb_found = False

        # reverse complement sequence for negative frames
        if frame.startswith('-'):
            dna = reverse_complement(dna)

        for i in range(int(frame[1])-1, seq_length - 2, 3):
            codon = dna[i:i+3]

            # start codon
            if codon in start_codons and start == None:
                start = i
                orf_codon_list.append(codon)

            # keep track if ambiguities are present
            elif [n for n in codon if n not in 'ACGT'] and start != None:
                if [n for n in codon if n in amb_nucl]:
                    amb_found = True

            # store codons
            elif codon not in stop_codons and start != None:
                orf_codon_list.append(codon)

            # stop codon
            elif codon in stop_codons and start != None:
                if amb_found == True and not allow_ambiguity:
                    amb_found = False
                    start = None
                    orf_codon_list = []
                    continue
                
                # check ORF length
                stop = i + 2
                orf_length = stop - start + 1
                if orf_length > minimum_length:

                    # rescale positions for negative strand
                    if frame.startswith('-'):
                        start, stop = rescale(start, stop, seq_length)
                        
                    orf_count += 1
                    orf_codon_list.append(codon)
                    temp_codon_count = Counter(orf_codon_list)
                    orfs_found[orf_count] = {'frame' : frame, 'start' : start,
                                            'stop' : stop, 'length' : orf_length,
                                            'codon count' : temp_codon_count,
                                            }

                    # add to total codon count
                    for cod, count in temp_codon_count.items():
                        try:
                            codon_count[cod] += count
                        except KeyError as err:
                            sys.exit("ERROR: Sequence provided as input contains non-valid nucleotides '{}'".format(err.args[0]))
                
                # reset
                start = None
                orf_codon_list = []

    return(orfs_found, codon_count)


def reverse_complement(string):

    """Returns the reverse complement strand for a given DNA sequence"""
    
    return(string[::-1].translate(str.maketrans('ATGC','TACG')))


def rescale(start, stop, length_sequence):

    """Rescales negative strand start-stop positions to positive strand indexing"""

    rescaled_start = length_sequence - stop - 1
    rescaled_stop = length_sequence - start - 1

    return(rescaled_start, rescaled_stop)


def calculate_codon_bias(counts_per_codon):

    """Calculates relative codon frequencies for codons encoding the same amino acid"""

    codon_bias = {}
    codons_per_aa = {
    'K' : ['AAA', 'AAG'],
    'N' : ['AAC', 'AAT'],
    'T' : ['ACA', 'ACC', 'ACG', 'ACT'],
    'R' : ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
    'S' : ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
    'I' : ['ATA', 'ATC', 'ATT'],
    'M' : ['ATG'],
    'Q' : ['CAA', 'CAG'],
    'H' : ['CAC', 'CAT'],
    'P' : ['CCA', 'CCC', 'CCG', 'CCT'],
    'L' : ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
    'E' : ['GAA', 'GAG'],
    'D' : ['GAC', 'GAT'],
    'A' : ['GCA', 'GCC', 'GCG', 'GCT'],
    'G' : ['GGA', 'GGC', 'GGG', 'GGT'],
    'V' : ['GTA', 'GTC', 'GTG', 'GTT'],
    '*' : ['TAA', 'TAG', 'TGA'],
    'Y' : ['TAC', 'TAT'],
    'C' : ['TGC', 'TGT'],
    'W' : ['TGG', ],
    'F' : ['TTC', 'TTT'],
    }

    # codon counts per family
    for codons in codons_per_aa.values():
        list_codons_per_aa = []
        for codon in codons:
            list_codons_per_aa.append(counts_per_codon[codon])

        # find max
        max_codon_per_aa = max(list_codons_per_aa)
        
        # calculate relative frequencies
        for codon in codons:
            try:
                codon_bias[codon] = counts_per_codon[codon]/max_codon_per_aa
            except ZeroDivisionError as err:
                break

    return(codon_bias)


def calculate_CAI(orfs_found, codon_bias):

    """Calculates CAI for ORFs using codon count tables and codon relative frequencies"""

    for orf, orf_data in orfs_found.items():
        numerator = 0
        denominator = 0
        for codon, count in orf_data['codon count'].items():
            numerator += count * codon_bias[codon]
            denominator += count
        CAI = numpy.exp(numerator / denominator)
        orfs_found[orf]['cai'] = CAI

    return(orfs_found)


def exclude_overlaps(orfs_found, filter_type):

    """Filters ORF overlaps based on the criteria indicated"""
    
    blocked_orfs = []
    sorted_orfs = OrderedDict(sorted(orfs_found.items(), key = lambda kvls: kvls[1][filter_type], reverse=True))
    overlapping_orfs = find_overlaps(orfs_found)
    for orf in sorted_orfs.keys():
        if orf in blocked_orfs:
            continue
        blocked_orfs.extend(overlapping_orfs[orf])

    return({ orf : orfs_found[orf] for orf in orfs_found if orf not in blocked_orfs})


def find_overlaps(orfs_found):

    """Associates overlapping ORFs based on their start-stop intervals"""

    # initialize dictionary to store associations for each ORF
    overlaps_found = { i : [] for i in range(1, len(orfs_found)+1) }
    processed_orfs = []

    # sort ORFs by ascending start position
    sorted_orfs = OrderedDict(sorted(orfs_found.items(), key = lambda kvls: kvls[1]['start']))

    # Loop through all combinations and store associations
    for orf, data in sorted_orfs.items():
        processed_orfs.append(orf)
        for orf_to_compare, data_to_compare in sorted_orfs.items():
            if orf_to_compare in processed_orfs:
                continue
            if data['stop'] >= data_to_compare['start']:
                overlaps_found[orf].append(orf_to_compare)
                overlaps_found[orf_to_compare].append(orf)
            else:
                break

    return(overlaps_found)


def format_output(name, seq, orfs_found, output_type, char_per_line):

    """Extracts information from ORF dictionary and performs output formatting"""

    output_lines = []
    orf_num = 0
    sorted_orfs = OrderedDict(sorted(orfs_found.items(), key = lambda kvls: kvls[1]['frame']))

    for orf in sorted_orfs.values():
        orf_num += 1
        if orf['frame'].startswith('-'):
            start, stop, length = (str(orf['stop']), str(orf['start']), str(orf['length']))
            cdna = reverse_complement(seq[orf['start']:orf['stop']+1])
            peptide = translate(cdna)[0:-1]
        else:
            start, stop, length = (str(orf['start']), str(orf['stop']), str(orf['length']))
            cdna = seq[orf['start']:orf['stop']+1]
            peptide = translate(cdna)[0:-1]

        # append to output
        output_lines.append('>{}_F{}_{:05d}\t{}:{}\t{}'.format(name, orf['frame'], orf_num, start, stop, length))
        if output_type == 'cdna':
            lines_to_add = [cdna[i:i+char_per_line] for i in range(0, len(cdna), char_per_line)]
            output_lines.extend(lines_to_add)
        else:
            lines_to_add = [peptide[i:i+char_per_line] for i in range(0, len(peptide), char_per_line)]
            output_lines.extend(lines_to_add)

    return(output_lines)


def translate(seq):

    """Translates a given DNA string to an aminoacidic sequence
        including the codons containing ambiguous nucleotide bases"""

    translated_seq = []
    nonamb_codons = ('GCN', 'CGN', 'AGR', 'AAY', 'GAY', 'TGY', 'CAR',
                    'GAR', 'GGN', 'CAY', 'ATH', 'YTR', 'CTN', 'AAR',
                    'TTY', 'CCN', 'TCN', 'AGY', 'ACN', 'TAY', 'GTN')
    geneticode = {
    'AAA' : 'K', 'AAC' : 'N', 'AAG' : 'K', 'AAT' : 'N',
    'ACA' : 'T', 'ACC' : 'T', 'ACG' : 'T', 'ACT' : 'T',
    'AGA' : 'R', 'AGC' : 'S', 'AGG' : 'R', 'AGT' : 'S',
    'ATA' : 'I', 'ATC' : 'I', 'ATG' : 'M', 'ATT' : 'I',
    'CAA' : 'Q', 'CAC' : 'H', 'CAG' : 'Q', 'CAT' : 'H',
    'CCA' : 'P', 'CCC' : 'P', 'CCG' : 'P', 'CCT' : 'P',
    'CGA' : 'R', 'CGC' : 'R', 'CGG' : 'R', 'CGT' : 'R',
    'CTA' : 'L', 'CTC' : 'L', 'CTG' : 'L', 'CTT' : 'L',
    'GAA' : 'E', 'GAC' : 'D', 'GAG' : 'E', 'GAT' : 'D',
    'GCA' : 'A', 'GCC' : 'A', 'GCG' : 'A', 'GCT' : 'A',
    'GGA' : 'G', 'GGC' : 'G', 'GGG' : 'G', 'GGT' : 'G',
    'GTA' : 'V', 'GTC' : 'V', 'GTG' : 'V', 'GTT' : 'V',
    'TAA' : '*', 'TAC' : 'Y', 'TAG' : '*', 'TAT' : 'Y',
    'TCA' : 'S', 'TCC' : 'S', 'TCG' : 'S', 'TCT' : 'S',
    'TGA' : '*', 'TGC' : 'C', 'TGG' : 'W', 'TGT' : 'C',
    'TTA' : 'L', 'TTC' : 'F', 'TTG' : 'L', 'TTT' : 'F', }

    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        try:
            translated_seq.append(geneticode[codon])
        except KeyError:
            if codon in nonamb_codons:
                codon = codon.translate(str.maketrans('NRYH','AACA'))
                translated_seq.append(geneticode[codon])
            else:
                translated_seq.append('X')

    return(''.join(translated_seq))


if __name__ == '__main__':
    main()
