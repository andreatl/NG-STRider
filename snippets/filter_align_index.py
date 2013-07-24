#ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp

'''
July 2010
Andrea Lai, Infinity Pharmaceuticals

Parse the alignment index found at the above link to provide a list of files for analysis using STRFinder.

== Usage ==
python filter_align_index.py -c [chrs] -p [populations] 

example options:
chrs: 1-4,X,Y
chrs: 1,3,5,8

populations:
CEU,CHB

*** Note, different chromosomes or populations should not be separated by any spaces!

== Output ==
Returns align_filtered.out

* This should be eventually merged directly into the 1KHG profiler.py workflow
'''

import csv, sys, re, optparse, string, os

def parse_file(file, coverage):
	"""Parses a tab-delimited file into a list."""
	sample = open(file, "r")
	readSample = csv.reader(sample, dialect="excel-tab")
	sampleList = []
	
	# filters by the coverage type (low, high, exon-targetted) selected
	for row in readSample:
		bam_loc = row[0]
		if re.search('('+coverage+')', bam_loc):
			sampleList.append(bam_loc)
	return sampleList

def create_csv(file, data):
	"""produce a tab delimited file from a list"""
	sample = open(file, 'w')
	writeSample = csv.writer(sample, lineterminator="\n", delimiter="\t")
	writeSample.writerows(data)
	sample.close()

def main():	
	p = optparse.OptionParser()
	p.add_option('--chr_filter', '-c', help="Filter by chromosomes -- no spaces! example: -c 1-4,X")
	p.add_option('--pop_filter', '-p', default = '', help="Filter by population! example: -p CEU,CHB")
	p.add_option('--output', '-o', default = 'align_filtered', help="Filter by population! example: -p CEU,CHB")
	options, args = p.parse_args()
	
	# select coverage type
	coverage = raw_input('What coverage: low, high, or exon-targetted? \n')
	
	# return all files of that coverage, prints to file
	try:
		list = parse_file('alignment.index', coverage)
	except IOError:
		# if alignment index not present, download from 1KHG FTP
		print 'Alignment index not found! Now downloading via FTP'
		os.system("wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/alignment.index")
		
	print str(len(list))+' '+coverage+ ' coverage files in the alignment index'
	
	new = []
	table = []
	
	# parse input population list to regex search string
	pop_filter = options.pop_filter.replace(',', '|')
	
	# filter by chromosome and population
	for item in list:
		if re.search(pop_filter, item):
			if options.chr_filter != None:
				if re.search("(chrom)["+options.chr_filter+"]", item):			
					split = re.split('[.]', item)
					table.append([item, split[4]])
			else:
				split = re.split('[.]', item)
				table.append([item, split[4]])
	
	# print filtered list to file
	create_csv(options.output+'.out', table)
	print str(len(table))+' files in the filtered set'
	
	# generate sample sizes of different populations		
	dict = {}
	for item in table:	
		if item[1] not in dict:
			dict[item[1]] = 1
		else:
			dict[item[1]] += 1

	print 'Distribution of populations:'
	print dict

if __name__ == '__main__':
	main()

