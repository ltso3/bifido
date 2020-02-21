#!/usr/bin/env python
#----------------USAGE-----------------------
#by Dan Barshis
#This script parses a .sam alignment file created using BWA and returns a tab-delimited text file with the 
#number of multiply aligned reads and good quality (above a certain mapping quality threshold and read length) singly aligned 
#reads for every contig/gene in the reference sequence database. 
#The input text for the command line is:
#countxpression.py mapqualitythreshold, lengththreshold, outputstatsfilename, anynumberofinputfiles ...  

def countxpression(infilename, threshold, lengththresh, zeroforheader1forno, tosort, outfilename, countsoutputname):
	#Opens an infile specified by the user. Should be a .sam file 
	IN = open(infilename, 'r')
	print(infilename)

	threshold=float(threshold) #inputs a mapping quality threshold for counting reads as sucessfully mapped
	lenthresh=float(lengththresh)
	
	#Opens an output text file as specified by user 
	OUT = open(outfilename, 'w')
	
	countsout = open(countsoutputname, 'a')
	
	totreads=0
	notaligned=0
	multi=0
	single=0
	aligned=0
	contigswzeros=0
	goodmaps=0
	
	#Starts a for loop to read through the infile line by line
	multicontigs={}
	contigs={}
	for line in IN:
	
		line=line.rstrip()
		cols=line.split('\t')
		if cols[0]=="@SQ": #Generate dictionary of contig names
			contigname=cols[1].split(':')[2] # SN:NC_011593.1:c2595355-2593832
			contigs[contigname]=[0, 0, 0]
			
		if cols[0][0] != '@' : #to skip past header
			if 'RG:' in line:
				columntocount=14
			else:
				columntocount=13 #UPDATE: changed from 13 to 11
			unmapped_flags = [77, 141, 133, 69, 165, 101]
			if float(cols[1]) in unmapped_flags: #this is to find reads that are unaligned
				notaligned+=1
				totreads+=1	
			else: #this is to find reads that had >=1 alignment
				aligned+=1
				totreads+=1
				if float(cols[columntocount].split(':')[2]) > 1: # for reads with >1 optimal hit
					multi+=1
					contigs[cols[2].split(':')[1]][1]+=1     #add a count to that contig for the 'top' hit for reads that optimally align to multiple contigs (all the other good hits are often not printed with current bwa settings)
				if float(cols[columntocount].split(':')[2]) == 1: # for reads with == 1 optimal hit
					single+=1
					if float(cols[4]) >= threshold and len(cols[9]) >= lenthresh:
						goodmaps+=1
						contigs[cols[2].split(':')[1]][0]+=1
	for item in contigs:
		contigs[item][2]=contigs[item][0]+contigs[item][1]
		if contigs[item][0] == 0:
			contigswzeros+=1
	
	contigsmatched=len(contigs.keys())-contigswzeros
	propqualaligned=float(goodmaps)/float(totreads)
		
	if int(zeroforheader1forno) == 0:
		countsout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % ('Filename', 'Total#Reads', 'NumContigsMatched', 'NumUnaligned', 'NumAligned', 'NumMultiAligned', 'NumSingleAligned', 'NumQualSingles', 'PropQualAligned'))
	else:
		countsout.write('\n%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f' % (infilename, totreads, contigsmatched, notaligned, aligned, multi, single, goodmaps, propqualaligned))

	
	OUT.write('ContigName'+'\t'+'UniqueTotReads'+'\t'+'MultiTotReads'+'\t'+'totalreadsgoodmapped')
	l=[]

	#This is for sorting the contigs into numerical order and is written specifically to handle a coral and symbiodinium mixed reference
	#you will need to modify this part based on your unique "gene" or "contig" names used in your reference assembly
	#this is for general sorting based on text sorting rules
	if tosort=='text':
		for key,value in contigs.items():
			l.append((key,value))
		l.sort()
		for item in l:
			# print item
			OUT.write('\n%s\t%d\t%d\t%d' % (item[0], item[1][0], item[1][1], item[1][2])) #writes each line of the tuple as separate tab delimited text
		
		OUT.close()

if __name__=="__main__":
	import sys
	# Usage countxpression(infilename, threshold, lengththreshold, zeroforheader1forno, tosort, outfilename, countsoutputname)
	filelist=sys.argv[4:]
	numfiles=0
	for file in filelist:
		numfiles+=1
		if numfiles==1:
			countxpression(file, sys.argv[1], sys.argv[2], 0, 'text', file[:-4]+'_counts.txt', sys.argv[3])
		if numfiles>1:
			countxpression(file, sys.argv[1], sys.argv[2], 1, 'text', file[:-4]+'_counts.txt', sys.argv[3])
