#!/usr/bin/python3

#Imports
import argparse
import sys
import re
import csv
import operator


#Constants
ERR_NONE = 0
ERR_USAGE = 1
ERR_IO = 2

#Open input file
#inFile = open(args.inFileName, 'r')



#Main 
if (__name__ == '__main__'):

	#Parse command line arguments
	parser = argparse.ArgumentParser(description='Pick Primer Pairs for Target Sequence Based on Product Length')
	
	parser.add_argument('-i', '--in', dest='inFileName', required=True,
						help='Input file name. The output file of the first script.')
	
	parser.add_argument('-o', '--out', dest='outFileName', required=True,
						help='Output file name.')
	
	parser.add_argument('-l', '--productlength', dest='productlength', required=True, 
						help='Desired Product Length') 
	
	parser.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true',
						help='Set verbose output.')
	
	args = parser.parse_args()
	
#Open input file
#if (verbose):
#	print('Opening input file: {0}'.format(args.inFileName)
	
#Write header line
	
#outFile.write('Kmer,Rev,Score')

x = []
row_index=0
f= open(args.outFileName, 'wt')

print("Forward Kmer, Reverse Kmer, Cumulative Score", file=f)
with open(args.inFileName, 'rt') as csvfile:
	samplereader = csv.reader(csvfile, delimiter =',', quotechar=' ')
	for row in samplereader:
              x.append(row)
             
for row in x:
        second_row=int(row_index+250)
        if second_row < len(x):
                score=0
                if (x[row_index][4] == "T"): score=score+1
                if (x[row_index][5] == "T"): score=score+1
                if (x[row_index][6] == "T"): score=score+1
                if (x[row_index][7] == "T"): score=score+1
                if (x[row_index][8] == "T"): score=score+1
                if (x[second_row][11] == "T"): score=score+1
                if (x[second_row][12] == "T"): score=score+1
                if (x[second_row][13] == "T"): score=score+1
                if (x[second_row][14] == "T"): score=score+1
                if (x[second_row][15] == "T"): score=score+1
                if (x[row_index][8] == "T" and x[second_row][15] == "T" and (float(x[row_index][9])-float(x[second_row][16])==0)): score = score+4
                if row_index !=0:
                        print("x[row_index][0],x[row_index][2],x[second_row][1],x[second_row][2],score", file=f)
                        #print(x[row_index][0],",",x[row_index][2],",",x[second_row][1],",",x[second_row][2],",",score, file=f)
                row_index=(row_index+1)

csvfile.close()
f.close()
#with open(args.outFileName, 'rt') as csvfile:
#	reader = csv.reader(csvfile, delimiter =',', quotechar=' ')
#	sortedlist = sorted(reader, key=operator.itemgetter(2), reverse=True)
#	csvfile.close()


#f= open(args.outFileName, 'w')
#print(sortedlist, file=f)
#f.close()

#for line in inFile:
#	tok = line.strip().split(',')
#	if not line.find("Kmer"):
#		continue
#	forward_kmer = tok[0]
#	if (tok[4] == 'T' and tok[5] == 'T'and tok[6] == 'T' and tok[7] == 'T' and tok[8] == 'T'):
#		print line
#		score =0
#		disqualify = False
	
#Create list of annotations
#print(x[1][1])
