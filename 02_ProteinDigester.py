#!/usr/bin/env python3
# Solution based on Perl code provided by Prof. Simon Hubbard

import argparse, sys, os, glob, re
from collections import OrderedDict

def main():

    #
    # parser object and command line input arguments
    #
    parser = argparse.ArgumentParser(description="Digest protein sequences to peptides.")
    parser.add_argument("input", help="Enter fasta file with protein sequences to digest.",
                        nargs='+')
    parser.add_argument("-o", "--output", help="Allows output to specified file.")
    parser.add_argument("-e", "--enzyme", default='t', choices=['t', 'k', 'v', 'r', 'x', 'e'],
                        help="Allows selecting different enzymes for digestion: trypsin (t), \
                        lys-c (k), Staph-v8 (v), arg-c (r), trypsin-acetyl (x) or glu-c (e).")
    parser.add_argument("-m", "--missed", type=int, default=0, choices=[0, 1, 2, 3, 4, 5], 
                        help="Allows setting different numbers of missed cleavages (0-5 are allowed).")
    parser.add_argument("-t", "--type", default='all', choices=['N', 'I', 'C', 'all'],
                        help="Allows filtering for peptide origin: n-terminal (N), internal (I), \
                        c-terminal (C) or all.")
    parser.add_argument("-b", "--batch", help="Enables batch mode for processing\
                        all '.fasta' files in a directory. Path to directory must\
                        be specified as input.", action="store_true")
    args = parser.parse_args()

    # check for already existing output file
    if args.output:
        if os.path.exists(args.output):
            sys.exit("ERROR: Please, select new output. File '{}' already exists!".format(args.output))
    
    # define enzymes
    enzname = {
        't' : 'trypsin',
        'k' : 'lys-c',
        'v' : 'Staph-v8',
        'r' : 'arg-c',
        'x' : 'trypsin-acetyl',
        'e' : 'glu-c'
    }

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
    # get digested peptides
    #
    all_pepts = []

    for name, seq in data.items():

        if re.match(r'[BOUJZX]', seq):
            print('## WARNING: protein {} contains unusual amino acids\n'.format(name))

        pept_count = 0
        peptides = digest(seq, args.enzyme)

        #
        # first deal with fully cleaved peptides, then miscleaved
        #
        for i in range(0, len(peptides)):
            
            # set peptide type
            type = 'I'   # internal peptide
            if i == len(peptides)-1 : type = 'C'   # C-terminal peptide
            if i == 0 : type = 'N'   # N-terminal peptide

            pept_count += 1

            if args.type == 'all':
                all_pepts.append('>{} peptide={} enzyme={} missed=0 type={}'.format(name, pept_count, enzname[args.enzyme], type))
                all_pepts.append(peptides[i])
            elif args.type == type:
                all_pepts.append('>{} peptide={} enzyme={} missed=0 type={}'.format(name, pept_count, enzname[args.enzyme], type))
                all_pepts.append(peptides[i])
            #
            # add miscleaved peptides if required
            #
            miscleaved_seq = peptides[i]
            for j in range(1, args.missed+1):
                if i+j < len(peptides):
                    miscleaved_seq += peptides[i+j]
                    pept_count += 1
                    if i+j == len(peptides)-1 : type = 'C'
                    if args.type == 'all':
                        all_pepts.append('>{} peptide={} enzyme={} missed={} type={}'.format(name, pept_count, enzname[args.enzyme], j, type))
                        all_pepts.append(miscleaved_seq)
                    else:
                        if args.type != type : continue
                        all_pepts.append('>{} peptide={} enzyme={} missed={} type={}'.format(name, pept_count, enzname[args.enzyme], j, type))
                        all_pepts.append(miscleaved_seq)

    #
    # output digested peptides
    #   
    if args.output:
        with open (args.output, 'w') as F: 
            F.write('\n'.join(all_pepts))
    else:
        print('\n'.join(all_pepts))


def retrieve_data(input_files, batch_flag):

    """Reads in individual files found in list or directory"""

    input_data = OrderedDict()

    if batch_flag:
        path = ' '.join(input_files)
        for file in glob.glob(os.path.join(path, '*')):
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

    seqs_in_file = OrderedDict()
    strings = []
    name = ''

    with open (filename, 'r') as F:
        for line in F:
            if line.startswith('>'):
                if strings:
                    seqs_in_file[name] = ''.join(strings).upper()
                    strings = []
                name = line.rstrip('\n')[1:]
                if name in seqs_in_file.keys():
                    sys.exit("ERROR: At least 2 sequences found with the same header '>{}...'".format(name))
            elif name:
                strings.append(line.rstrip('\n'))
        
        if strings:
            seqs_in_file[name] = ''.join(strings).upper()
        else:
            sys.exit("ERROR: No sequences found in '{}'. Only files in FASTA format are supported.".format(filename))

    return(seqs_in_file)


def digest(sequence, enzyme):
    
    """Fragments sequence at all enzyme cutting sites"""

    specificity = {
        't' : r'([KR])([^P\n])',
        'k' : r'(K)([^P\n])',
        'r' : r'(R)([^P\n])',
        'v' : r'([DE])([^P\n])',
        'e' : r'(E)([^P\n])',
        'x' : r'(R)([^P\n])'
    }

    # run re.sub twice on sequence to effectively 'cut' overlapping matches
    sequence = re.sub(specificity[enzyme], r'\1\n\2', sequence)
    sequence = re.sub(specificity[enzyme], r'\1\n\2', sequence)
    peptides = sequence.split('\n')

    return(peptides)


if __name__ == '__main__':
    main()
