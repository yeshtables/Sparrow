#!/usr/bin/python3

import argparse
from Bio import SeqIO
import re
import os.path

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Read Splitter')
     
    parser.add_argument('inFileList', metavar='INPUT_FILE', nargs='+',
                        help='List of input files containing sequence reads to be split.')
     
    parser.add_argument('-f', '--format', dest='format', default='fasta',
                        help='Format of input files (see INPUT_FILE).')
    
    parser.add_argument('-p', '--pattern', dest='pattern',
                        help='Pattern to split.')
    
    parser.add_argument('-m', '--matchout', dest='matchOut', default='match.fasta',
                        help='File with matches (default = match.fasta).')
    
    parser.add_argument('-n', '--nomatchout', dest='noMatchOut', default='nomatch.fasta',
                        help='File with no matches (default = nomatch.fasta')
     
    args = parser.parse_args()
    
    # Set temporary variables
    matchList = []
    mismatchList = []
    
    # Open files
    matchFile = open(args.matchOut, 'w')
    mismatchFile = open(args.noMatchOut, 'w')
    
    # Process input files
    for inFile in args.inFileList:
        seq = SeqIO.parse(inFile, args.format)
        
        for record in seq:
            if re.search(args.pattern, record.description) != None:
                matchList.append(record)
            else:
                mismatchList.append(record)
    
    # Write sequences
    SeqIO.write(matchList, matchFile, 'fasta')
    SeqIO.write(mismatchList, mismatchFile, 'fasta')
    
    # Close output
    matchFile.close()
    mismatchFile.close()
    