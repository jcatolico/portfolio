#!/usr/bin/env python
#usage: python undupe_comp.py inputFilePath outputFilePath t/f
#example: python undupe_comp.py /home/jct/sampleCoverage.tsv /home/jct/sampleOutput.tsv t
#t = true if file contains header line, f = false if it doesn't.
#add input check

import commands
import re
from fnmatch import fnmatch
import sys,os

inputFile = sys.argv[1]
outputFile = sys.argv[2]
header = sys.argv[3]

#define the files
#Sample experimental design with 4 experimental groups with 4 biological replicates each.
#These filenames correspond to SAM/BAM files from STAR aligner output.
dirdict={0:"Group1_Rep1",1:"Group1_Rep2",2:"Group1_Rep3",3:"Group1_Rep4",4:"Group2_Rep1",5:"Group2_Rep2",6:"Group2_Rep3",7:"Group2_Rep4",8:"Group3_Rep1",9:"Group3_Rep2",10:"Group3_Rep3",11:"Group3_Rep4",12:"Group4_Rep1",13:"Group4_Rep2",14:"Group4_Rep3",15:"Group4_Rep4"}

#define threshold for filtering (using Sanger/Phred33 base quality)
threshold=25

#helper function to return a range for subsetting sequence/quality scores.
#input: cigar string
#output: (start,end position) and truncated cigar string for subsetting sequence/quality scores. 
#a range of (3,6) acknowledges positions 3, 4, 5, and 6 to be covered by this read.
#0 is considered the first position.
#Example: singleCigarToRange("3S44M") => ([3, 46], '44M') 
def singleCigarToRange(cigarStringInput):
	currentPosition=0 
	range=[]
	trimCigar=""
	openCoverage = False
	cigarElements = filter(None,re.compile("([0-9]+)([A-Z]+)").split(cigarStringInput))
	while len(cigarElements)>0:
		cigarNumber=int(cigarElements.pop(0))
		cigarSymbol=cigarElements.pop(0)
		if cigarSymbol in 'HS':
			if cigarSymbol is 'H':
				trimCigar+='H'
			if openCoverage is True:
				range.append(currentPosition-1) #close out the range if it hits a tailing S or H 
				return range,trimCigar
			elif cigarSymbol is 'S': #add to the start position if it hits a S at the beginning
				currentPosition+=cigarNumber
		elif cigarSymbol in 'MDI': 
			if openCoverage is False: #open range up at the first M
				openCoverage=True
				range.append(currentPosition)
			if cigarSymbol in 'MI':
				currentPosition+=cigarNumber #Adds M to the range. 
				#if it hits an I, in this case, just add to the range as well. This will be corrected for in the following step.
				#if it's a D, no range is added.
			trimCigar+=str(cigarNumber)+cigarSymbol 
	#close out the coverage when cigar string ends.
	if openCoverage is True:
		range.append(currentPosition-1)
	return range,trimCigar	

#Main program starts here.
#Open the coverage data file
#Coverage data was kept in a TSV file in the following format:
#1_3342593_A_T	1,4	3,3	2,5	3,5
#1_3342598_A_T	1,2	2,4	1,3	1,2
#The first column contains information about a detected SNP - chromosome 1, position 3342593, A to T
#The following columns contain reference/alternate base counts for each site in each sample
#1 A, 4 T
with open(inputFile,'r') as getfile:
	#The results will be output in this file:
	with open(outputFile,'w') as writefile:
		#optionally skip a header line... 
		if 't' is header:
			writefile.write(getfile.next())
		
		#For every candidate site (each row)
		for line in getfile:
			splitLine = line.strip().split('\t') #0 = variant id, 1-16 = sample groups 
			splitVariant = splitLine[0].split('_') #0 = chr, 1 = pos, 2 = ref, 3 = alt 
			saveVariant = splitLine[0]
			splitLine = splitLine[1:] #subset the line to remove the variant element.
			
			#for each sample group... ignore "0,0" as there should be no data there. 
			for na in [inum for inum, entry in enumerate(splitLine) if entry is not "0,0"]:
				#get the sam file data for coverage at chromosome:editsite-editsite 
				data = commands.getoutput('samtools view /home/jct/proj/alignment/STAR/'+dirdict[na]+'/split_reads.bam \"'+str(splitVariant[0])+':'+str(splitVariant[1])+'-'+str(splitVariant[1])+'\"')
				#need to parse each sam entry...
				entryDict={} #create dictionary to hold the information. Keys will be [start,length], data will be [sequence, quality scores]

				refAltCounter=[0,0]
				#parse the sam data and trim using coordinates from singleCigarToRange
				#builds entryDict using position, length, and cigar as dictionary keys, and the actual sequence data and quality scores as values.
				for samline in data.strip().split("\n"):
					if len(data)>0:
						splitsamline = samline.split("\t") 
						genomicStartPos = int(splitsamline[3]) #start
						cigarStr = splitsamline[5] #cigar
						samSeq = splitsamline[9] #sequence
						samQual = splitsamline[10] #quality score
						trimCoords,trimCigar = singleCigarToRange(cigarStr)
						trimStart,trimEnd = trimCoords
						trimLength=(int(trimEnd)-int(trimStart))+1 #length. 1 to 3 is 3 - 1 + 1 = length 3
						trimSeq=samSeq[int(trimStart):int(trimEnd)+1] #second number passed to range is exclusive, must add 1.
						trimQual=samQual[int(trimStart):int(trimEnd)+1] #second number passed to range is exclusive, must add 1.
					
						trimVariantPos = int(splitVariant[1])-genomicStartPos
						
						try: 
							if ord(trimQual[trimVariantPos])-33 >= threshold:
								#passes threshold, put in dictionary
								if trimSeq[trimVariantPos] in splitVariant[2:]:
									variant = trimSeq[trimVariantPos] is splitVariant[3] #1 if alt, 0 if ref
									key = (variant,genomicStartPos,trimCigar,trimLength)
									if key not in entryDict:
										entryDict[key] = []
									entryDict[key].append([trimSeq,trimQual]) 
						except IndexError:
							pass
							
				#All entries are now in the dictionary.
				for key in entryDict:
					if len(entryDict[key]) is 1: 
						#only 1 entry here, means we can count it and leave straight away
						refAltCounter[key[0]]+=1
					else: 
						#Multiple similar entries at the same position must now be compared basewise
						#in order to properly remove duplicates based on SNP content.
						#Here, "words" are created using dissimilar bases
						#that pass the base quality threshold, and are compared at the end.
						#For example, comparing the following:
						#AATACAC -> AA
						#AATGCTC -> GT
						#AATACAC -> AA
						#   ^ ^ dissimilar bases at position 3 and 5 produce two letter words assuming all bases pass quality threshold
						#Bases that don't pass the quality threshold can't be used to distinguish reads from one another,
						#however, the read may still contain information at other positions.
						wordarray=["" for w in xrange(len(entryDict[key]))]
						for i in range(0,key[3]): #for each position... (build wordarray in this loop)
							if not all(entry[0][i] == entryDict[key][0][0][i] for entry in entryDict[key]): #If there is a mismatch ... 
								qualcheck={}
								qualseq=[]
								recordMismatch=False
								for j in range(0,len(entryDict[key])): #iterate over the list of lists[seq,qual]
									if ord(entryDict[key][j][1][i])-33>=threshold and entryDict[key][j][0][i] is not "N":
										if entryDict[key][j][0][i] not in qualcheck:
											qualcheck[entryDict[key][j][0][i]]=0
										qualcheck[entryDict[key][j][0][i]]+=1
										qualseq.append(entryDict[key][j][0][i])
									else: #if it doesn't pass quality check, or it is an N...
										qualseq.append('?') #? is used later on as a wildcard when checking for duplicates
								#Must delete unsupported bases from qualcheck
								#Each dissimilar base must be supported by at least 2 different reads with sufficient quality.
								#The following code uses information from ALL samples to determine whether 
								#a base is supported.
								for unconfirmedBase in [qual for qual in qualcheck if qualcheck[qual] < 2]: 
									#must check all reads now...
									#build qualcheck but across all reads...
									evaluatePos = int(key[1])+i 
									data2 = commands.getoutput('samtools view /home/jct/proj/alignment/STAR/'+dirdict[na]+'/split_reads_2.bam \"'+str(splitVariant[0])+':'+str(evaluatePos)+'-'+str(evaluatePos)+'\"')
									#need to parse each sam entry...
									confirmed = False
									for samline2 in data2.split("\n"):
										if len(data2)>0:
											splitsamline2 = samline2.split("\t") 
											genomicStartPos2 = int(splitsamline2[3]) #start
											cigarStr2 = splitsamline2[5] #cigar
											samSeq2 = splitsamline2[9] #sequence
											samQual2 = splitsamline2[10] #quality score
											trimCoords2,trimCigar2 = singleCigarToRange(cigarStr2)
											trimStart2,trimEnd2 = trimCoords2
											targetTrimBase = evaluatePos-genomicStartPos2
											try: 
												if samSeq2[targetTrimBase] is qual and ord(samQual2[targetTrimBase])-33 >= threshold:
													#base is confirmed now, record this position. 
													recordMismatch=True
													confirmed = True
													break
											except IndexError:
												pass
		
									#following is run if we didn't break upon confirmation...
									if not confirmed:
										del(qualcheck[unconfirmedBase])	
										qualseq = [base if base is not unconfirmedBase else '?' for base in qualseq]	
								#done with this position, add information to wordarray if a difference remains 
								if recordMismatch:
									for idx,sequence in enumerate(qualseq):
										wordarray[idx]+=qualseq[idx]
						#now check word array for dupes
						if len(wordarray[0])>0:
							for iternum in range(len(wordarray)): 
								iternum2=iternum+1
								while iternum2 < len(wordarray):
									if fnmatch(wordarray[iternum],wordarray[iternum2]):
										del(wordarray[iternum2])
									else:
										iternum2+=1
							#for each unique word, refAltCounter[key[0]]+=1
							#The output of this program is similar to the input,
							#where Ref/Alt counts are presented in a #,# format.
							#The numbers represent the number of unduplicated reads that support that site with sufficient quality.
							for unduplicatedWord in wordarray:
								refAltCounter[key[0]]+=1
						else: #there were no mismatches detected 
							refAltCounter[key[0]]+=1
				splitLine[na] = str(refAltCounter[0])+","+str(refAltCounter[1])
			writefile.write(saveVariant+'\t'+'\t'.join(splitLine)+'\n')
	writefile.close()		
getfile.close()
