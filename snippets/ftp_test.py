'''
Quick check of the Pysam install and configuration for FTP.
If no errors, you are all set!
'''

import pysam
region = ['1', 1, 10000]

url = 'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00116/alignment/HG00116.chrom1.ILLUMINA.bwa.GBR.low_coverage.20100517.bam'

sam = pysam.Samfile(url, "rb")

list = []
for item in sam.fetch(region[0], region[1], region[2]):
	list.append([item.qname, item.seq])
	
print list[0:10]
