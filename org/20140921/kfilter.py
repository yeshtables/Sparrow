#!/usr/bin/python3

# Imports
import argparse
import sys
import subprocess
import gzip

# Constants
ERR_NONE = 0
ERR_USAGE = 1
ERR_IO = 2

# Function: err
def err(msg, ret):
    print("{0}: {1}".format(sys.argv[0], msg))
    sys.exit(ret)

# Function: nextKmer
def nextKmer(inFile):
    
    while True:
        line = inFile.readline()
        
        if (isinstance(line, bytes)):
            line = line.decode()
    
        if (line == ''):
            return None
        
        tok = line.strip().split('\t')
        
        if (len(tok) >= 2):
            return tok[0]

# Main
if __name__ == '__main__':
    
    # Declarations
    filterSet = set()
    
    # Get command line arguments
    parser = argparse.ArgumentParser(description='K-mer filter')

    parser.add_argument('inFileList', metavar='INPUT_FILE', nargs='+',
                        help='List of input files.')
    #parser.add_argument('-i', '--in', dest='inFileName', required=True,
    #                    help='Sequence file to read.')
    
    parser.add_argument('-k', '--ksize', dest='kSize', default=21,
                        help='Size of output k-mers.')
    
    parser.add_argument('-s', '--subksize', dest='subKSize', default=16,
                        help='Size of sub-k-mers. Search k-mers and filter if any sub-kmers match the filtering set.')
    
    parser.add_argument('-f', '--format', dest='format', default='fasta',
                        help='Format of the input file (fasta, fastq).')

    parser.add_argument('-o', '--out', dest='outFileName', required=True,
                        help='Output file of filtered k-mers.')

    parser.add_argument('-r', '--filter', dest='filterFileName', required=True,
                        help='File of k-mers to be filtered. K-mers must be of the sub-kmer size.')
    
    parser.add_argument('-t', '--tempdir', dest='tempDirName', default='temp',
                        help='Location of temporary files.')
    
    parser.add_argument('-e', '--kanloc', dest='kanalyze', default='kanalyze',
                        help='Location of KAnalyze. "count" and "stream" should be found in this directory.')
    
    parser.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true',
                        help='Set verbose output')
    
    parser.add_argument('-V', '--noverbose', dest='verbose', action='store_false',
                        help='No verbose output (default)')
    
    parser.add_argument('-z', '--filtergz', dest='filterGz', default=True, action='store_true',
                        help='Filter file is gzipped (default).')
    
    parser.add_argument('-Z', '--nofiltergz', dest='filterGz', action='store_false',
                        help='Filter file is not gzipped.')

    args = parser.parse_args()
    
    verbose = args.verbose
    
    # Set locations
    targetSubKc = args.tempDirName + "/targetsub.kc"
    
    if (verbose):
        print('Target sub k-mer file: {0}'.format(targetSubKc))
        print('Opening output file {0}'.format(args.outFileName))
  
    # Open output file
    try:
        outFile = open(args.outFileName, 'w')

    except OSError as ex:
        err('Error opening output file "{0}": {1}'.format(args.outFileName, ex.strerror), ERR_IO)
        
        if (verbose):
            print('Opening filter file {0}'.format(args.filterFileName))

    # Open filter file
    try:
        if (args.filterGz):
            filterFile = gzip.open(args.filterFileName, 'r')
        else:
            filterFile = open(args.filterFileName, 'r')

    except OSError as ex:
        err('Error opening filter file "{0}": {1}'.format(args.filterFileName, ex.strerror), ERR_IO)
    
    if (verbose):
        print('Writing sub-kmers of target sequence')
    
    # Write sub-kmer of target sequence
    kanProc = subprocess.Popen([args.kanalyze + '/count', '-k', str(args.subKSize), '-f', args.format, '-o', targetSubKc] + args.inFileList)
    
    retcode = kanProc.wait()
    
    if (retcode != 0):
        err('Target sub-kmer process terminated with code {0}'.format(retcode), retcode)
    
    if (verbose):
        print('Opening full sub k-mer target file {0}'.format(targetSubKc))
    
    # Open sub-kmer file
    try:
        subKFile = open(targetSubKc, 'r')
    
    except OSError as ex:
        err('Error opening sub k-mer file "{0}": {1}'.format(args.filterFileName, ex.strerror), ERR_IO)
    
    if (verbose):
        print('Filtering sub k-mers and loading hash')
    
    # Traverse filter and sub k-mer file
    subKmer = nextKmer(subKFile)
    filterKmer = nextKmer(filterFile)
    
    while (subKmer != None and filterKmer != None):
        
        if (subKmer == filterKmer):
            filterSet.add(subKmer)
            
            subKmer = nextKmer(subKFile)
            filterKmer = nextKmer(filterFile)
        
        elif (subKmer < filterKmer):
            subKmer = nextKmer(subKFile)
        
        else:
            filterKmer = nextKmer(filterFile)
    
    if (verbose):
        print('Read {0} k-mers into hash memory'.format(len(filterSet)))
        print('Closing sub k-mer file')
            
    subKFile.close()
    
    if (verbose):
        print('Closing filter file')
    
    filterFile.close()
    
    if (verbose):
        print('Processing target full k-mers')
    
    # Filter target by sub k-mers
    kanProc = subprocess.Popen([args.kanalyze + '/count', '-k', str(args.kSize), '-f', args.format, '--stdout'] + args.inFileList, stdout=subprocess.PIPE)
    
    for line in kanProc.stdout:
        tok = line.strip().split()
        writeKmer = True
        
        if (len(tok) < 2):
            continue
        
        kanSubProc = subprocess.Popen([args.kanalyze + '/stream', '-k', str(args.subKSize), '--input', tok[0], '--stdout'], stdout=subprocess.PIPE)
        
        for subLine in kanSubProc.stdout:
            subLine = subLine.decode().strip()
                  
            if (subLine in filterSet):
                writeKmer = False
        
        outFile.write('{0},{1}\n'.format(tok[0].decode(), str(1 if writeKmer else 0)))
        
        kanSubProc.wait()
        
        if (kanSubProc.returncode != 0):
            err('KAnalyze per-k-mer sub-kmer process died with return code {0}'.format(kanSubProc.retcode), kanSubProc.retcode)
        
    kanProc.wait()
    
    if (kanProc.returncode != 0):
        err('KAnalyze input kmer process died with return code {0}'.format(kanSubProc.retcode), kanSubProc.retcode)
    
    if (verbose):
        print('Closing output file')

    outFile.close()
    