#!/usr/bin/python3

# Imports
import argparse
import subprocess
import sys
import os
from Bio import SeqIO

# Function: err
def err(msg, ret):
    errMsg(msg)
    sys.exit(ret)

# Function: errMsg
def errMsg(msg):
    print('{0}: {1}'.format(sys.argv[0], msg))

# Main
if (__name__ == '__main__'):
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Find k-mer-wise conservation by reference.')
    
    parser.add_argument('inFileList', metavar='INPUT_FILE', nargs='+',
                        help='List of input files.')
    
    parser.add_argument('-e', '--kanloc', dest='kanalyze', default='kanalyze',
                        help='Location of KAnalyze. "count" and "stream" should be found in this directory.')
    
    parser.add_argument('-f', '--format', dest='format', default='fasta',
                        help='Format of the input file (default = fasta).')
    
    parser.add_argument('-k', '--ksize', dest='kSize', default=21,
                        help='Size of output k-mers.')
    
    parser.add_argument('-m', '--multi', dest='multiSeq', default=True, action='store_true',
                        help='Count input sequences from multi-sequence files (default).')
    
    parser.add_argument('-M', '--nomulti', dest='multiSeq', action='store_false',
                        help='Count input sequences from individual files.')
    
    parser.add_argument('-o', '--out', dest='outFileName', required=True,
                        help='Output file name.')
    
    parser.add_argument('-t', '--tempdir', dest='tempDirName', default='temp',
                        help='Location of temporary files.')
    
    parser.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true',
                        help='Set verbose output.')
    
    parser.add_argument('-V', '--noverbose', dest='verbose', action='store_false',
                        help='Unset verbose output (default).')
    
    args = parser.parse_args()
    verbose = args.verbose
    
    tempDirName = args.tempDirName + "/temp.cons." + str(os.getpid())
    
    # Read each input file and save number of occurrences
    recordCount = 0
    
    if verbose:
        print('Getting k-mer frequencies')
    
    kmerFrequency = {}
    
    for file in args.inFileList:
        
        if verbose:
            print('Opening sequence file: ' + file)
        
        if args.multiSeq:
            
            inFile = open(file)
            
            for (index, record) in enumerate(SeqIO.parse(inFile, args.format)):
                recordCount += 1
                
                if verbose:
                    print('Processing read: ' + str(recordCount))
                
                kanProc = subprocess.Popen([args.kanalyze + '/count', '-k', str(args.kSize), '-t', tempDirName, '--stdout', '--input', str(record.seq)], stdout=subprocess.PIPE)
                
                for line in kanProc.stdout:
                    line = line.decode().strip()
                    kmer = line.split('\t')[0]
                
                    if kmer in kmerFrequency:
                        kmerFrequency[kmer] += 1
                    else:
                        kmerFrequency[kmer] = 1
            
                kanProc.wait()
        
                if (kanProc.returncode != 0):
                    err('KAnalyze input kmer process died with return code {0}'.format(kanProc.returncode), kanProc.returncode)
            
            inFile.close()
            
        else:
            recordCount += 1
            
            kanProc = subprocess.Popen([args.kanalyze + '/count', '-k', str(args.kSize), '-t', tempDirName, '-f', args.format, '--stdout', file], stdout=subprocess.PIPE)
        
            for line in kanProc.stdout:
                line = line.decode().strip()
                kmer = line.split('\t')[0]
            
                if kmer in kmerFrequency:
                    kmerFrequency[kmer] += 1
                else:
                    kmerFrequency[kmer] = 1
        
            kanProc.wait()
    
            if (kanProc.returncode != 0):
                err('KAnalyze input kmer process died with return code {0}'.format(kanProc.returncode), kanProc.returncode)
    
    # Write reference conservation
    if verbose:
        print('Opening output file: ' + args.outFileName)
    
    outFile = open(args.outFileName, 'w')
    
    outFile.write('#kmer,score\n')
    
    for kmer in kmerFrequency:
        outFile.write('{0},{1}\n'.format(kmer, kmerFrequency[kmer] / recordCount))
    
    if verbose:
        print('Closing output file')
    
    outFile.close()
    
    # Remove temporary directory (fails if it contains files, which it will if k-analyze or removing the temporary file fails)
    if (verbose):
        print('Removing temporary directory: ' + tempDirName)
    
    try:
        os.rmdir(tempDirName)
        
    except:
        errMsg('Warning: Error removing temporary directory: ' + tempDirName + ': ' + sys.exc_info()[0])
    