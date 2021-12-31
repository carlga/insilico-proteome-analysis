#!/usr/bin/env python3

import argparse, sys, os, glob, re
from collections import OrderedDict

def main():

    #
    # parser object and command line input arguments
    #
    parser = argparse.ArgumentParser(description="Calculate mass for peptide sequences.")
    parser.add_argument("input", help="Enter fasta file with peptide sequences.",
                        nargs='+')
    parser.add_argument("-o", "--output", help="Allows output to specified file.")
    parser.add_argument("-m", "--masses", type=int, default=0, choices=[0, 1, 2, 3, 4, 5, 6, 7], 
                        help="Allows to consider different modifications for calculating peptide mass:\
                        (0) Monoisotopic,\
                        (1) Average,\
                        (2) Monoisotopic + Phosphorylation (Ser, Tyr and Thr),\
                        (3) Average + Phosphorylation (Ser, Tyr and Thr),\
                        (4) Monoisotopic + Oxidation (Met),\
                        (5) Average + Oxidation (Met),\
                        (6) Monoisotopic + Phosphorylation (Ser, Tyr and Thr) + Oxidation (Met),\
                        (7) Average + Phosphorylation (Ser, Tyr and Thr) + Oxidation (Met)")
    parser.add_argument("-t", "--type", default='all', choices=['N', 'I', 'C', 'all'],
                        help="Allows filtering masses for peptide origin: n-terminal (N), internal (I), \
                        c-terminal (C) or all.")
    parser.add_argument("-b", "--batch", help="Enables batch mode for processing\
                        all '.fasta' files in a directory. Path to directory must\
                        be specified as input.", action="store_true")
    args = parser.parse_args()

    # check for already existing output file
    if args.output:
        if os.path.exists(args.output):
            sys.exit("ERROR: Please, select new output. File '{}' already exists!".format(args.output))

    #
    # read peptides from file
    #
    try:
        data = retrieve_data(args.input, args.batch)
    except FileNotFoundError as err:
        sys.exit("ERROR: File '{}' has not been found (check name or path to file)".format(err.filename))
    except IsADirectoryError as err:
        sys.exit("ERROR: Please, provide full path to file. To enable batch mode use '-b' option.")

    #
    # calculate peptide masses
    #
    pept_masses = ['id\tprotein\tpeptide\tmass-to-charge\tz\tmissed\tenzyme\ttype\tsequence']

    for header, seq in data.items():
        metadata = parseHeader(header)
        mass_to_charge = peptideMass(seq, mass_code=args.masses, ion_mass=1, ion_charge=1)
        
        if args.type == 'all':
            pept_masses.append('{}\t{}\t{}\t{:.3f}\t{}\t{}\t{}\t{}\t{}'.format(metadata['id'], metadata['protein'], metadata['peptide'], mass_to_charge, 1, metadata['missed'], metadata['enzyme'], metadata['type'], seq))
        elif args.type == metadata['type']:
            pept_masses.append('{}\t{}\t{}\t{:.3f}\t{}\t{}\t{}\t{}\t{}'.format(metadata['id'], metadata['protein'], metadata['peptide'], mass_to_charge, 1, metadata['missed'], metadata['enzyme'], metadata['type'], seq))
    
    #
    # output peptide masses
    #   
    if args.output:
        with open (args.output, 'w') as F: 
            F.write('\n'.join(pept_masses))
    else:
        print('\n'.join(pept_masses))


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


def parseHeader(header):

    """Extracts metadata from space-separated header"""

    header = header.split()
    id = header.pop(0)
    meta = dict(data.split('=') for data in header)
    meta['id'] = id

    return(meta)

def peptideMass(sequence, mass_code, ion_mass=1, ion_charge=1):

    """Calculates mass-to-charge ratio for a given peptide sequence"""

    # standard residue masses from http://www.matrixscience.com/help/aa_help.html
    aa_masses = {
        'A' : (71.0371, 71.08, 71.0371, 71.08, 71.0371, 71.08, 71.0371, 71.08),
        'C' : (103.0092, 103.14, 103.0092, 103.14, 103.0092, 103.14, 103.0092, 103.14),
        'D' : (115.0269, 115.09, 115.0269, 115.09, 115.0269, 115.09, 115.0269, 115.09),
        'E' : (129.0426, 129.12, 129.0426, 129.12, 129.0426, 129.12, 129.0426, 129.12),
        'F' : (147.0684, 147.18, 147.0684, 147.18, 147.0684, 147.18, 147.0684, 147.18),
        'G' : (57.0215, 57.05, 57.0215, 57.05, 57.0215, 57.05, 57.0215, 57.05),
        'H' : (137.0589, 137.14, 137.0589, 137.14, 137.0589, 137.14, 137.0589, 137.14),
        'I' : (113.0841, 113.16, 113.0841, 113.16, 113.0841, 113.16, 113.0841, 113.16),
        'K' : (128.0950, 128.17, 128.0950, 128.17, 128.0950, 128.17, 128.0950, 128.17),
        'L' : (113.0841, 113.16, 113.0841, 113.16, 113.0841, 113.16, 113.0841, 113.16),
        'M' : (131.0405, 131.19, 131.0405, 131.19, 147.0354, 147.19, 147.0354, 147.19),
        'N' : (114.0429, 114.10, 114.0429, 114.10, 114.0429, 114.10, 114.0429, 114.10),
        'P' : (97.0528, 97.12, 97.0528, 97.12, 97.0528, 97.12, 97.0528, 97.12),
        'Q' : (128.0586, 128.13, 128.0586, 128.13, 128.0586, 128.13, 128.0586, 128.13),
        'R' : (156.1011, 156.19, 156.1011, 156.19, 156.1011, 156.19, 156.1011, 156.19),
        'S' : (87.0320, 87.08, 166.9983, 167.06, 87.0320, 87.08, 166.9983, 167.06),
        'T' : (101.0477, 101.10, 181.014, 181.08, 101.0477, 101.10, 181.014, 181.08),
        'V' : (99.0684, 99.13, 99.0684, 99.13, 99.0684, 99.13, 99.0684, 99.13),
        'W' : (186.0793, 186.21, 186.0793, 186.21, 186.0793, 186.21, 186.0793, 186.21),
        'Y' : (163.0633, 163.18, 243.0296, 243.16, 163.0633, 163.18, 243.0296, 243.16),
        '\s' : (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        "*" : (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    }

    # iterate through residues to get total mass
    mass = 0
    for aa in sequence:
        mass += aa_masses[aa][mass_code]
    
    # add masses of terminating groups (H and OH)
    if mass_code % 2 == 0:
        mass += 18.0106
    else:
        mass += 18.0153
    
    # calculate m/z
    mass_to_charge = (mass + ion_mass) / ion_charge

    return(mass_to_charge)
    

if __name__ == '__main__':
    main()
