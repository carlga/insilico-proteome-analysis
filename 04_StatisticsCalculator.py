#!/usr/bin/env python3

import argparse, sys, os
import numpy as np
from collections import OrderedDict

def main():

    #
    # parser object and command line input arguments
    #
    parser = argparse.ArgumentParser(description="Calculate mass spectrometry statistics for ionized peptides.")
    parser.add_argument("input", help="Enter tab-delimited file with peptide information.")
    parser.add_argument("-o", "--output", help="Allows output to specified file.")
    parser.add_argument("-r", "--range", type=int, nargs=2, default=[1000,1500], 
                        help="Allows specifying a custom range for peptide detection. Default is 1000-1500 Da.")
    parser.add_argument("-s", "--step", type=int, default=5, 
                        help="Allows to select step size for sliding window analysis. Default is 5 Da.")
    parser.add_argument("-w", "--window", type=int, default=10, 
                        help="Allows to select window size for sliding window analysis. Default is 10 Da.")
    args = parser.parse_args()

    # check for already existing output file
    if args.output:
        if os.path.exists(args.output):
            sys.exit("ERROR: Please, select new output. File '{}' already exists!".format(args.output))
    
    # check detection range
    if args.range[0] < 0 or args.range[1] < 0:
        sys.exit('ERROR: Detection range limits must be positive!')
    if args.range[0] >= args.range[1]:
        sys.exit('ERROR: Invalid range limits! Upper limit must be of a higher value than lower one.')

    # check window and step size
    if args.step <= 0 or args.window <= 0:
        sys.exit('ERROR: Window and step size must be positive!')
    if args.step > args.window:
        sys.exit('ERROR: Step must be lower or equal to the window size!')

    #
    # read peptide data from file
    #
    try:
        data = read_table(args.input)
    except FileNotFoundError as err:
        sys.exit("ERROR: File '{}' has not been found (check name or path to file)".format(err.filename))
    except IsADirectoryError as err:
        sys.exit("ERROR: Please, provide full path to file.")
    
    # sort peptides by mass-to-charge ratio (ascending)
    sorted_data = OrderedDict(sorted(data.items(), key = lambda kvls: kvls[1]['mass-to-charge'], reverse=False))

    # bin peptides
    bins = bin_peptides(sorted_data, args.range, args.step, args.window)

    #
    # output bins
    #   
    if args.output:
        with open (args.output, 'w') as F: 
            F.write('\n'.join(bins))
    else:
        print('\n'.join(bins))


def read_table(filename):

    """Reads in table from tab-delimited file"""

    data = OrderedDict()

    with open(filename, 'r') as F:
        lines = F.readlines()
    
    header = lines.pop(0).split('\t')
    for line in lines:
        line = line.split('\t')
        pept_data = { header[i]:line[i] for i in range(len(header)) }
        data[pept_data['protein']+'.'+pept_data['peptide']] = pept_data

    return(data)


def bin_peptides(data, range, step, window):

    """Bins peptides using sliding window method"""

    bins = ['range (Da)\taverage\tcount\tpeptides']

    for mass in np.arange(range[0], range[1], step):
        peptides = []
        if mass+window > range[1]:
            break
        for peptide, info in data.items():
            if mass > float(info['mass-to-charge']):
                continue
            elif mass+window >= float(info['mass-to-charge']):
                peptides.append(peptide)
            else:
                break
        if peptides:
            bin_sum = sum([float(data[peptide]['mass-to-charge']) for peptide in peptides])
            bin = '{}:{}\t{:.3f}\t{}\t{}'.format(mass, mass+window, bin_sum/len(peptides), len(peptides), ';'.join(peptides))
        else:
            bin = '{}:{}\t{}\t{}\t{}'.format(mass, mass+window, 0, len(peptides), 'NA')
        bins.append(bin)
    
    return(bins)


if __name__ == '__main__':
    main()
