#!/usr/bin/python3
#Imports
import argparse
import sys
import re
from collections import defaultdict

#Constants
ERR_NONE = 0
ERR_USAGE = 1
ERR_IO = 2

#Main 
if (__name__ == '__main__'):

	#Parse command line arguments
	parser = argparse.ArgumentParser(description='Pick Primer Pairs for Target Sequence Based on Product Length')
	
	parser.add_argument('-i', '--in', dest='inFileName', required=True,
						help='Input file name. The output file of the first script.')
	
	parser.add_argument('-q', '--query', dest='queryFileName', required=True,
						help='Output file name.')
	
	parser.add_argument('-l', '--productlength', dest='productlength', required=False, 
						help='Desired Product Length') 
	
	parser.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true',
						help='Set verbose output.')
	
	args = parser.parse_args()
	
	#Open input file
inFile = open(args.inFileName, 'rt')
	
	#Write header line
	#outFile.write('Kmer,Rev,Score')

#Create Dictionary                      
a = defaultdict(float)
                      
for line in inFile:
        tok = line.strip().split(',')
        if not line.find("#"):
            continue
        if len(tok)<2:
            continue
        a[tok[0]] = float(tok[1])
                      
print (a['TCTAACGCAACATAATAAACT'])

query = open(args.queryFileName)

for line in query:
	tok = line.strip().split(',')
	if not line.find("#"):
	continue
	print (line
