#!/usr/bin/python3

# Imports
import argparse
import sys
import re

# Constants
ERR_NONE = 0
ERR_USAGE = 1
ERR_IO = 2

# Function: err
def err(msg, ret):
    print("{0}: {1}".format(sys.argv[0], msg))
    sys.exit(ret)

# Class: AnnoElement
class AnnoElement:
    
    def __init__(self, title, maxScore):
        self.title = title
        self.maxScore = maxScore
    
    def check(self, kmer):
        return True

# Class: AnnoStartAA (k-mer does not start with AA)
class AnnoStartAA(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'No Start AA', 0)
    
    def check(self, kmer):
        if (re.match('^AA.*', kmer) == None):
            return True
        
        return False

# Class: AnnoPos10AU (K-mer has A or U at position 10) 
class AnnoPos10AU(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'Pos 10 A/U', 0.05)
    
    def check(self, kmer):
        if (re.match('^.{9}[AU].*', kmer) != None):
            return 0.05
        
        return 0

# Class AnnoPos19 (k-mer has A or U at position 19) 
class AnnoPos19AU(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'Pos 19 A/U', 0.15)
    
    def check(self, kmer):
        if (re.match('^.{18}[AT].*', kmer) != None):
            return 0.15
        
        return 0

# Class AnnoPos15to20AUx3 (k-mer 3 A or U between positions 15 and 20 inclusive)
class AnnoPos15to20AUx3(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'Pos 15-20 3 A/U', 0.50)
    
    def check(self, kmer):
        if (len(re.findall('A|T', kmer[14:20])) >= 3):
            return 0.50
        
        return 0

# Class AnnoPos17to19AUCount (0.10 for each A/U between positions 17 and 19 inclusive)
class AnnoPos17to19AUCount(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'Pos 17-19 A/U Ct', 0.30)
    
    def check(self, kmer):
        return (len(re.findall('A|T', kmer[16:19])) * 0.1)

# Class AnnoGC30to52 (False for any k-mer without GC content between 30% and 52% inclusive)
class AnnoGC30to52(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'GC 30-52%', 0)
    
    def check(self, kmer):
        gcContent = len(re.findall('G|C', kmer)) / len(kmer)
        
        if (gcContent >= 0.30 and gcContent <= 0.52):
            return True
        
        return False

# Class AnnoRestrictionSites (False for any with restriction sites)
class AnnoRestrictionSites(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'No Restr Site', 0)
    
    def check(self, kmer):
        if (re.match('.*(GGTACC)|(GAATTC)|(CTCGAG)|(CATATG)|(ACTAGT)|(GGTAC)|(GAATT)|(GTACC)|(TACC)|(CTAGT).*', kmer) == None):
            return True
        
        return False

# Class Anno7GCRun (False for any with a run of 7 or more GC)
class Anno7GCRun(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'No 7 GC Run', 0)
    
    def check(self, kmer):
        if (re.match('.*[GC]{7}.*', kmer) == None):
            return True
        
        return False
    
# Main
if (__name__ == '__main__'):
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Annotate 21-mers for RNAi effectiveness')
    
    parser.add_argument('-i', '--in', dest='inFileName', required=True,
                        help='Input file name. A list of k-mers with no counts.')
    
    parser.add_argument('-o', '--out', dest='outFileName', required=True,
                        help='Output file name.')
    
    parser.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true',
                        help='Set verbose output.')
    
    parser.add_argument('-V', '--noverbose', dest='verbose', action='store_false',
                        help='Unset verbose output (default).')
    
    args = parser.parse_args()
    
    verbose = args.verbose
    
    # Create list of annotations
    annoList = []
    annoList.append(AnnoStartAA())
    annoList.append(AnnoPos10AU())
    annoList.append(AnnoPos19AU())
    annoList.append(AnnoPos15to20AUx3())
    annoList.append(AnnoPos17to19AUCount())
    annoList.append(AnnoGC30to52())
    annoList.append(AnnoRestrictionSites())
    annoList.append(Anno7GCRun())
    
    # Get score divisor (normalizes so the max score is 1)
    maxScore = 0
    for anno in annoList:
        maxScore += anno.maxScore
    
    # Open input file
    if (verbose):
        print('Max score: {0}'.format(maxScore))
        print('Opening input file: {0}'.format(args.inFileName))
    
    try:
        inFile = open(args.inFileName, 'r')

    except OSError as ex:
        err('Error opening input file "{0}": {1}'.format(args.outFileName, ex.strerror), ERR_IO)
    
    # Open output file
    if (verbose):
        print('Opening output file: {0}'.format(args.outFileName))
    
    try:
        outFile = open(args.outFileName, 'w')

    except OSError as ex:
        err('Error opening output file "{0}": {1}'.format(args.outFileName, ex.strerror), ERR_IO)
    
    # Write header line
    outFile.write('kmer,position,filter')
    
    for anno in annoList:
        outFile.write(',{0}'.format(anno.title))
    
    outFile.write(',score\n')

    # Read each line
    for line in inFile:
        tok = line.strip().split(',')
        score = 0
        disqualify = False
                
        if (len(tok) < 2):
            continue
        
        kmer = tok[0]
        filterPass = tok[1]
        position = tok[2]
        
        if (filterPass == '0'):
            disqualify=True
        
        outFile.write('{0},{1},{2}'.format(kmer, position, 'T' if filterPass == '1' else 'F'))
        
        # Write filter results
        for anno in annoList:
            annoScore = anno.check(kmer)
            
            if (isinstance(annoScore, bool)):
                if (annoScore == False):  # Disqualify k-mer if a filter returns False
                    disqualify = True
                    annoScore = 'F'
                
                else:  # For True, do not disqualify, but do not add to score
                    annoScore = 'T'
                
            else:  # All else, add score
                score += annoScore
            
            outFile.write(',{0}'.format(annoScore))
        
        if (disqualify):
            score = 0
        else:
            score /= maxScore
        
        outFile.write(',{0}\n'.format(score))
    
    # Close files
    if (verbose):
        print('Closing input file')
    
    inFile.close()
    
    if (verbose):
        print('Closing output file')
    
    outFile.close()
