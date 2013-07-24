'''
Create custom config files for STRfinder.
Input is a text file with a single line, comma-separate list of probes to be included, or just a series of probes as a command line argument.

Not designating '-o' only prints results to console.

== usage ==
python sts_mapping.py -l listin.txt -o out
python sts_mapping.py D3S1358 TH01 D21S11

== example ==
listin.txt
D3S1358, TH01, D21S11, D18S51, Penta E, D5S818, D13S317, D7S820, D16S539

== Note ==
1. Need to download the STSMap files from UCSC Genome Browser
http://hgdownload.cse.ucsc.edu/downloads.html

2. After tab-delimited file is generated, you still need to enter the repeat sequences by hand to a 4th tab delimited column. Repeat sequences can be found here:
http://www.cstl.nist.gov/strbase/str_fact.htm

If it is a complex repeat, leave the repeat column blank, but make sure it has an entry in the 'complex' config file.
'''

#!/usr/bin/env python

import csv, sys
import optparse

#usage = "Usage: %prog [options]"

def parse_stsMap(file):
	"""Parses an stsMap (from GenomeBrowser) to return ___"""
	sample = open(file, "r")
	readSample = csv.reader(sample, dialect="excel-tab")

	sampleDict = {}
	for row in readSample:
		chromosome = row[0].replace("chr","")
		start = row[1]
		end = row[2]
		name = row[3]		
		sampleDict[name] = [chromosome, start, end]

	return sampleDict
	sample.close()

def parse_cytoMap(file):
	#would require regEx to get all of, say, 21q5.* or a range of bands - more thought needed. 
	"""Parses an stsMap (from GenomeBrowser) to return ___"""
	sample = open(file, "r")
	readSample = csv.reader(sample, dialect="excel-tab")

	sampleDict = {}
	for row in readSample:
		chromosome = row[0].replace("chr","")
		start = row[1]
		end = row[2]
		name = chromosome+row[3]		
		sampleDict[name] = [chromosome, start, end]

	return sampleDict
	sample.close()
	
def parse_alias(file):
	"""Parses an stsMap (from GenomeBrowser) to return ___"""
	sample = open(file, "r")
	readSample = csv.reader(sample, dialect="excel-tab")

	sampleDict = {}
	for row in readSample:
		alias = row[0]
		trueName = row[2]		
		sampleDict[alias] = trueName

	return sampleDict
	sample.close()

def convert_sts(arg_list, options):
	v_no = options.version
	mapping = parse_stsMap("mapFiles/stsMap_hg"+v_no+".txt")
	aliasDict = parse_alias("mapFiles/stsAlias.txt")

	config_list = {}
	exclude = []
	
	out_list = []
	for item in arg_list:
		if item in mapping:
			out_list.append([item] + mapping[item])
		elif item in aliasDict:
			trueName = aliasDict[item]
			out_list.append([item] + mapping[trueName])
			#print item, trueName, mapping[trueName]
		else:
			exclude.append(item)	

	return out_list, exclude

def convert_cyto(arg_list):
	mapping = parse_cytoMap("mapFiles/cytoBand_hg18.txt")
	config_list = {}
	exclude = []
	
	out_list = []
	for item in arg_list:
		if item in mapping:
			out_list.append([item] + mapping[item])
		elif item in aliasDict:
			trueName = aliasDict[item]
			out_list.append([item] + mapping[trueName])
			#print item, trueName, mapping[trueName]
		else:
			exclude.append(item)	

	return out_list, exclude
	
def create_csv(file, data):
	"""produce a tab delimited file from a list"""
	sample = open(file, 'w')
	writeSample = csv.writer(sample, lineterminator="\n", delimiter="\t")
	writeSample.writerows(data)
	sample.close()
	
def main():
	p = optparse.OptionParser()
	p.add_option('--version', '-v', default="18")
	p.add_option('--output', '-o', default = "custom")
	p.add_option('--list', '-l', action="store_true")
	options, args = p.parse_args()
	print "mapped to human genome version " + options.version
	
	sts_list = args
	
	if options.list == True:
		file = args[0]
		f = open(file, "r")
		sts_list = f.read().replace(' ','').split(',')
	
	out_list, exclude = convert_sts(sts_list, options)

	create_csv(options.output+"_v"+options.version+".config", out_list)
	if len(exclude) > 1:
		print "Items excluded from config: ", str(list(exclude)).replace("'", '').strip('[]')
		
if __name__ == '__main__':
	main()
