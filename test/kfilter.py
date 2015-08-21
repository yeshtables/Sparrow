#!/usr/bin/python3

# Imports
import argparse
import sys
import subprocess
import gzip
import os

# Constants
ERR_NONE = 0
ERR_USAGE = 1
ERR_IO = 2

kSize = 21
subKSize = 16

# Function: err
def err(msg, ret):
    errMsg(msg)
    sys.exit(ret)

# Function: errMsg
def errMsg(msg):
    print("{0}: {1}".format(sys.argv[0], msg))

# Function: nextKmer
def nextKmer(inFile):
    
    line = inFile.readline()
        
    if (line == ''):
        return None
    
    return line[0:kSize]

# Main
if __name__ == '__main__':
    
    # Declarations
    filterSet = set()
    
    # Get command line arguments
    parser = argparse.ArgumentParser(description='K-mer filter')

    parser.add_argument('inFileList', metavar='INPUT_FILE', nargs='+',
                        help='List of input files.')
    
    parser.add_argument('-k', '--ksize', dest='kSize', default=21,
                        help='Size of output k-mers.')
    
    parser.add_argument('-s', '--subksize', dest='subKSize', default=16,
                        help='Size of sub-k-mers. Search k-mers and filter if any sub-kmers match the filtering set.')
    
    parser.add_argument('-f', '--format', dest='format', default='fasta',
                        help='Format of the input file (default = fasta).')

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
    # Note: Append PID to tempDirName. Use it as the temp file location for KAnalyze
    tempDirName = args.tempDirName + "/run." + str(os.getpid())
    targetSubKc = tempDirName + "/targetsub.kc"
    
    kSize = args.kSize
    subKSize = args.subKSize

    
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
            filterFile = gzip.open(args.filterFileName, 'rb')
        else:
            filterFile = open(args.filterFileName, 'rb')

    except OSError as ex:
        err('Error opening filter file "{0}": {1}'.format(args.filterFileName, ex.strerror), ERR_IO)
    
    if (verbose):
        print('Writing sub-kmers of target sequence')
    
    # Write sub-kmer of target sequence
    kanProc = subprocess.Popen([args.kanalyze + '/count', '-k', str(subKSize), '-t', tempDirName, '-f', args.format, '-o', targetSubKc] + args.inFileList)
    
    retcode = kanProc.wait()
    
    if (retcode != 0):
        err('Target sub-kmer process terminated with code {0}'.format(retcode), retcode)
    
    if (verbose):
        print('Opening full sub k-mer target file {0}'.format(targetSubKc))
    
    # Open sub-kmer file
    try:
        subKFile = open(targetSubKc, 'rb')
    
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
    
    # Remove sub-k-mer file
    if (verbose):
        print('Removing temporary file: ' + targetSubKc)
    
    try:
        os.remove(targetSubKc)
        
    except:
        errMsg('Warning: Error closing temporary file: ' + targetSubKc + ': ' + sys.exc_info()[0])
    
    # Remove temporary directory (fails if it contains files, which it will if k-analyze or removing the temporary file fails)
    if (verbose):
        print('Removing temporary directory: ' + tempDirName)
    
    try:
        os.rmdir(tempDirName)
        
    except:
        errMsg('Warning: Error removing temporary directory: ' + tempDirName + ': ' + sys.exc_info()[0])
    
    if (verbose):
        print('Closing filter file')
    
    filterFile.close()
    
    if (verbose):
        print('Processing target full k-mers')
    
    # Filter target by sub k-mers
    kanProc = subprocess.Popen([args.kanalyze + '/stream', '-k', str(kSize), '-f', args.format, '--stdout', '--index'] + args.inFileList, stdout=subprocess.PIPE)
    
    for line in kanProc.stdout:
        line = line.decode().strip()
        tok = line.split('\t')
        
        if (len(tok) < 2):
            continue
        
        writeKmer = True
        
        for i in range(len(tok[0]) - subKSize + 1):
            if (tok[0][i:(i + 16)] in filterSet):
                writeKmer = False
        
        outFile.write('{0},{1},{2}\n'.format(tok[0], str(1 if writeKmer else 0), str(tok[1])))
    
    if (verbose):
        print('Closing output file')

    outFile.close()
    