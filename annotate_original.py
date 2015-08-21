#!/usr/bin/python3

# Imports
import argparse
import sys
import re

# Constants
ERR_NONE = 0
ERR_USAGE = 1
ERR_IO = 2

# Globals
revDict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# Function: Reverse complement
def revComplStr(s):
    sr = ''
    l = len(s)
    
    for i in range(l):
        sr = revDict[s[i]] + sr
    
    return sr

# Function: err
def err(msg, ret):
    print("{0}: {1}".format(sys.argv[0], msg))
    sys.exit(ret)

# Class: AnnoElement
class AnnoElement:
    
    def __init__(self, title, maxScore):
        self.title = title
        self.maxScore = maxScore
    
    def check(self, kmer, pos):
        return True

# Class: AnnoEndGC (k-mer does not end with GC)
class AnnoEndGC(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'No GC Clamp', 0)
    
    def check(self, kmer, pos):
        if (re.match('^.*GC$', kmer) == None):
            return True
        
        return False

# Class: AnnoDinucleotide (k-mer does not contain 3x Dinucleotide Repeat)
class AnnoDiNuc(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'No 3x Dinucleotide Repeats', 0)
    
    def check(self, kmer, pos):
        return (re.search('(..)\\1\\1', kmer) == None)

# Class: AnnoHomopolymer (k-mer does not contain 4x Homopolymer Repeat)
class AnnoHomoPol(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'No 4x Hompolymer Repeats', 0)
    
    def check(self, kmer, pos):
        return (re.search('(.)\\1{3}', kmer) == None)

# Class: AnnoSelfComplementarity (k-mer is not Self-Complementary)
class AnnoSelfComp(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'Not Self-Complementary', 0)
    
    def check(self, kmer, pos):
        return (re.search('(...).*-.*\\1', '{0}-{1}'.format(kmer, revComplStr(kmer))) == None)
	
# Class: AnnoReverseKmer (prints reverse kmer)
class AnnoReverseKmer(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'Not Self-Complementary', 0)
    count = 0
    def check(self, kmer, pos):
        if pos >= 300:
            return (kmer)
		
# Class: AnnoPos10AU (K-mer has A or U at position 10) 
class AnnoPos10AU(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'Pos 10 A/U', 0.05)
    
    def check(self, kmer, pos):
        if (re.match('^.{9}[AU].*', kmer) != None):
            return 0.05
        
        return 0

# Class AnnoPos19 (k-mer has A or U at position 19) 
class AnnoPos19AU(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'Pos 19 A/U', 0.15)
    
    def check(self, kmer, pos):
        if (re.match('^.{18}[AU].*', kmer) != None):
            return 0.15
        
        return 0

# Class AnnoPos15to20AUx3 (k-mer 3 A or U between positions 15 and 20 inclusive)
class AnnoPos15to20AUx3(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'Pos 15-20 3 A/U', 0.50)
    
    def check(self, kmer, pos):
        if (len(re.findall('A|U', kmer[14:20])) >= 3):
            return 0.50
        
        return 0

# Class AnnoPos17to19AUCount (0.10 for each A/U between positions 17 and 19 inclusive)
class AnnoPos17to19AUCount(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'Pos 17-19 A/U Ct', 0.30)
    
    def check(self, kmer, pos):
        return (len(re.findall('A|U', kmer[16:19])) * 0.1)

# Class AnnoGC40to60 (False for any k-mer without GC content between 40% and 60% inclusive)
class AnnoGC40to60(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'GC 40-60%', 0)
    
    def check(self, kmer, pos):
        gcContent = len(re.findall('G|C', kmer)) / len(kmer)
        
        if (gcContent >= 0.40 and gcContent <= 0.60):
            return True
        
        return False

# Class AnnoRestrictionSites (False for any with restriction sites)
class AnnoRestrictionSites(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'No Restr Site', 0)
    
    def check(self, kmer, pos):
        if (re.match('.*(GGUACC)|(GAAUUC)|(CUCGAG)|(CAUAUG)|(ACUAGU)|(GGUAC)|(GAAUU)|(GUACC)|(UACC)|(CUAGU).*', kmer) == None):
            return True
        
        return False

# Class Anno7GCRun (False for any with a run of 7 or more GC)
class Anno7GCRun(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'No 7 GC Run', 0)
    
    def check(self, kmer, pos):
        if (re.match('.*[GC]{7}.*', kmer) == None):
            return True
        
        return False

class AnnoCons(AnnoElement):
    def __init__(self, consFileName):
        AnnoElement.__init__(self, 'Conservation', 0)
        
        self.consDict = {}
        
        self.consFileName = consFileName
        self.__readConsFile()
    
    def check(self, kmer, pos):
        
        if kmer in self.consDict:
            return (self.consDict[kmer], 0)
        
        return (0, 0)
    
    def __readConsFile(self):
        lineCount = 0
        
        consFile = open(self.consFileName, 'r')
        
        for line in consFile:
            lineCount += 1
            
            line = line.strip()
            
            # Skip empty lines and comments
            if (len(line) == 0 or line[0] == '#'):
                continue
            
            tok = line.split(',')
            kmer = tok[0].upper().replace('T', 'U')
            conScore = float(tok[1])
            
            if conScore < 0 or conScore > 1:
                raise ValueError('Bad conservation record on line {0} ({1}): Entry must be between 0 and 1: '.format(lineCount, self.consFileName, tok[1]))
            
            self.consDict[kmer] = conScore

class AnnoConsFilter(AnnoElement):
    def __init__(self, consDict, threshold=0.9):
        AnnoElement.__init__(self, 'Cons Thresh', 0)
        
        if consDict is None:
            consDict = {}
        
        self.consDict = consDict
        self.threshold = threshold
    
    def check(self, kmer, pos):
        
        if kmer in self.consDict:
            if (self.consDict[kmer] > self.threshold):
                return True
        
        return False

class GC_Content(AnnoElement):
    def __init__(self):
        AnnoElement.__init__(self, 'kmer GC Content', 0)
    
    def check(self, kmer, pos):
        gcContent = len(re.findall('G|C', kmer)) / len(kmer)
        
       # if (gcContent >= 0.30 and gcContent <= 0.52):
       #     return True
        
       # return False
        return(gcContent, 0)

# Class AnnoGene (Determine if k-mer is genic)
class AnnoGene(AnnoElement):
    def __init__(self, kSize, geneFileName, filterNonGenic=True):
        AnnoElement.__init__(self, 'Gene', 0)
        
        self.kSize = kSize
        self.geneFileName = geneFileName
        
        self.__readGeneFile()
        self.__notFound = not bool(filterNonGenic)
    
    def check(self, kmer, pos):
        
        for gene in self.__geneList:
            if pos > gene[2]:
                return ('', self.__notFound)
            
            if pos >= gene[1]:
                if pos <= gene[2]:
                    return (gene[0], True)
        
        return ('', self.__notFound)
    
    def __readGeneFile(self):
        
        lineCount = 0
        
        geneList = None
        
        geneFile = open(self.geneFileName, 'r')
        
        for line in geneFile:
            lineCount += 1
            
            line = line.strip()
            
            # Skip empty lines and comments
            if (len(line) == 0 or line[0] == '#'):
                continue

            # Split and check            
            tok = line.split('\t')
            tok[1] = int(tok[1])
            tok[2] = int(tok[2])
            
            if (len(tok) != 3):
                raise ValueError('Bad gene record on line {0} ({1}): Must contain 3 tab-separated entries'.format(lineCount, self.geneFileName))
            
            if tok[1] <= 0 or tok[2] <= 0:
                raise ValueError('Bad gene record on line {0} ({1}): Gene positions must not be 0 or negative'.format(lineCount, self.geneFileName))
            
            if tok[1] >= tok[2]:
                raise ValueError('Bad gene record on line {0} ({1}): Gene must stop at a position greater than its start'.format(lineCount, self.geneFileName))
            
            # Check gene name
            tok[0] = tok[0].strip()
            
            if len(tok[0]) == 0:
                raise ValueError('Bad gene record on line {0} ({1}): Gene name must not be empty'.format(lineCount, self.geneFileName))
            
            # Correct k-mer (k-mers must fall within the gene)
            tok[2] = tok[2] - self.kSize + 1
            
            if tok[2] < 0:
                raise ValueError('Bad gene record on line {0} ({1}): Gene length is shorter than the k-mer size (no k-mers in gene)'.format(lineCount, self.geneFileName))
            
            # Create list if empty
            if geneList is None:
                geneList = [tok]
                continue

            # Add to gene list (sorted)        
            for i in range(len(geneList)):
                if tok[2] < geneList[i][1]:
                    pass # Move to next case
                
                elif tok[1] > geneList[i][2]:
                    if i > 0:
                        geneList.insert(i - 1, tok)
                    else:
                        geneList.insert(0, tok)
                
                else:
                    raise ValueError('Bad gene record on line {0} ({1}): Record overlaps with another entry ({0}: {1} - {2}, {3}: {4} - {5}'.format(lineCount, self.geneFileName, tok[0], tok[1], tok[2], geneList[i][0], geneList[i][1], geneList[i][2]))
        
        self.__geneList = geneList
        
        geneFile.close()
    
# Main
if (__name__ == '__main__'):
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Annotate 21-mers for RNAi effectiveness')
    
    parser.add_argument('-a', '--annotateonly', dest='annotateOnly', default=False, action='store_true',
                        help='When including a gene file or conservation scores, use it for annotation only (do not alter scores).')
    
    parser.add_argument('-A', '--noannotateonly', dest='annotateOnly', action='store_false',
                        help='When filtering by a gene file, annotate and alter scores (default).')
    
    parser.add_argument('-c', '--conservation', dest='consFileName', default=None,
                        help='File of k-mers and conservation scores.')
    
    parser.add_argument('-g', '--gene', dest='geneFileName', default=None,
                        help='Gene file name')
    
    parser.add_argument('-i', '--in', dest='inFileName', required=True,
                        help='Input file name. A list of k-mers with no counts.')
    
    parser.add_argument('-o', '--out', dest='outFileName', required=True,
                        help='Output file name.')
    
    parser.add_argument('-r', '--reverse', dest='revCompl', default=False, action='store_true',
                        help='Show k-mer and its reverse complement')
    
    parser.add_argument('-R', '--noreverse', dest='revCompl', action='store_false',
                        help='Do not show the reverse complement of k-mers (default)')
    
    parser.add_argument('-s', '--score', dest='filtScore', default=False, action='store_true',
                        help='Write score regardless of disqualify filters. This option adds another score column.')
    
    parser.add_argument('-S', '--noscore', dest='filtScore', action='store_false',
                        help='Do not write an additional column for the score regardless of disqualifying filters (default).')
    
    parser.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true',
                        help='Set verbose output.')
    
    parser.add_argument('-V', '--noverbose', dest='verbose', action='store_false',
                        help='Unset verbose output (default).')
    
    args = parser.parse_args()
    
    verbose = args.verbose
    revCompl = args.revCompl
    filtScore = args.filtScore
    
    # Create list of annotations
    annoList = []
    annoList.append(AnnoEndGC())
    annoList.append(AnnoDiNuc())
    annoList.append(AnnoHomoPol())
    annoList.append(AnnoSelfComp())
    #annoList.append(AnnoReverseKmer())
    #annoList.append(AnnoPos17to19AUCount())
    annoList.append(AnnoGC40to60())
    annoList.append(GC_Content())
    #annoList.append(AnnoRestrictionSites())
    #annoList.append(Anno7GCRun())
    
    if args.consFileName is not None:
        consFilter = AnnoCons(args.consFileName)
        annoList.append(consFilter)
        
        if not args.annotateOnly:
            annoList.append(AnnoConsFilter(consFilter.consDict))
    
    if args.geneFileName is not None:
        annoList.append(AnnoGene(21, args.geneFileName, not args.annotateOnly))
    
    # Get score divisor (normalizes so the max score is 1)
    maxScore = 0
    for anno in annoList:
        maxScore += anno.maxScore
    
    if maxScore == 0:
        maxScore = 1
 
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

    outFile.write('Kmer,Rev,Index,Filter')
    
    for anno in annoList:
        outFile.write(',Fwd {0}'.format(anno.title))
		
    if (filtScore):
        outFile.write(',Fwd Score,Fwd No Filt Score')
    else:
        outFile.write(',Fwd Score')
    
        for anno in annoList:
            outFile.write(',Rev {0}'.format(anno.title))
	
    if (filtScore):
        outFile.write(',Rev Score,Rev No Filt Score\n')
    else:
        outFile.write(',Rev Score\n')

    # Read each line
    for line in inFile:
        tok = line.strip().split(',')
        score = 0
        disqualify = False
                
        if (len(tok) < 2):
            continue
        
        kmer = tok[0].upper().replace('T', 'T')
        revKmer = revComplStr(kmer)
        filterPass = tok[1]
        position = int(tok[2])
        
        if (filterPass == '0'):
            disqualify=True
        

        outFile.write('{0},{1},{2},{3}'.format(kmer, revComplStr(kmer), position, 'T' if filterPass == '1' else 'F'))
       
        # Write filter results (fwd primer)
        for anno in annoList:
            annoRet = anno.check(kmer, position)
            
            if isinstance(annoRet, tuple):
                # Returned tuple
                if isinstance(annoRet[1], bool):
                    annoScore = (annoRet[0], 0, annoRet[1])
                else:
                    annoScore = (annoRet[0], annoRet[1], True)
            else:
                # Returned score or boolean
                if isinstance(annoRet, bool):
                    annoScore = ('T' if annoRet else 'F', 0, annoRet)
                else:
                    annoScore = (str(annoRet), annoRet, True)
            
            # Write and increment score
            outFile.write(',{0}'.format(annoScore[0]))
            score = 0
            
            # Check disqualify flag
            if not annoScore[2]:
                disqualify = True
                
        score /= maxScore
        
        # Correct score
        if (disqualify):
            filteredScore = 0
        else:
            filteredScore = score
        
        # Write score
        if (filtScore):
            outFile.write(',{0},{1}'.format(filteredScore, score))
        else:
            outFile.write(',{0}'.format(filteredScore))
    
		# Write filter results (rev primer)
        score = 0
        disqualify = False
        
        for anno in annoList:
            annoRet = anno.check(revKmer, position)
            
            if isinstance(annoRet, tuple):
                # Returned tuple
                if isinstance(annoRet[1], bool):
                    annoScore = (annoRet[0], 0, annoRet[1])
                else:
                    annoScore = (annoRet[0], annoRet[1], True)
            else:
                # Returned score or boolean
                if isinstance(annoRet, bool):
                    annoScore = ('T' if annoRet else 'F', 0, annoRet)
                else:
                    annoScore = (str(annoRet), annoRet, True)
            
            # Write and increment score
            outFile.write(',{0}'.format(annoScore[0]))
            score += annoScore[1]
            
            # Check disqualify flag
            if not annoScore[2]:
                disqualify = True
                
        score /= maxScore
        
        # Correct score
        if (disqualify):
            filteredScore = 0
        else:
            filteredScore = score
        
        # Write score
        outFile.write(',{0}\n'.format(filteredScore))
		
    # Close files
    if (verbose):
        print('Closing input file')
    
    inFile.close()
    
    if (verbose):
        print('Closing output file')
    
    outFile.close()
