'''
Quick check!

Quick set of commands. Get the distribution read length of fragments at the denoted region
Requires Pysam.

This should be integrated into STRFinder as an initial QC evaluation.
'''

import pysam
region = ['11', 2148865, 2148950]
url = "ftp://ftp.sanger.ac.uk/pub/rd/humanSequences/CV.bam"

sam = pysam.Samfile(url, "rb")

list = []
for item in sam.fetch(region[0], region[1], region[2]):
	list.append(item.seq)
	
d = {}
for i in list:
	if len(i) in d:
		d[len(i)] += 1
	else:
		d[len(i)] = 1
		
print d		