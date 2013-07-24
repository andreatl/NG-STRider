'''
July 2010
Andrea Lai, Infinity Pharamceuticals
Find common (CODIS) STR repeats in next-gen sequencing data for cell typing, etc.

strider_1khg.py
For use in parsing 1000 Genomes Project data.

==Usage==
Note that all input files must be sorted bam files.
Output is a tab-delimited file containing the input filename, query id of the probe, the probe id, repeat count, sequence (with 'n' nucleotides before and after the repeat sequence), rad length, and quality control flags.

python str_finder.py -o [output file] -n [# nt] list_file.txt

see README.txt for more details
'''

from operator import itemgetter
import csv, sys, pysam, re, optparse, string
from decimal import *

class Result():
	def __init__(self, list, n):
		""" infer genotypes by returning the top 'n' items in your sorted list"""
		if len(list)>0:
			self.counts = [[k for k,v in item] for item in list][0][0:n]
		else:
			self.counts = 0
	
	def summary(self, flag = 0):
		""" reformats the list for printout in summary file"""
		if flag == 0:
			results = [str(item) for item in self.counts]
			self.result = ", ".join(results)
		else:
			self.result = ''
		return self.result

class Match():
	def __init__(self, re_object, read, b):
		self.repeat = re_object.group(0)
		self.match_len = len(self.repeat)
		self.start = re_object.start()
		self.end = re_object.end()
		self.start_buffered = self.start - b
		if self.start_buffered < 0:
			self.start_buffered = 0
		if self.start_buffered < 0:
			self.start_buffered = 0

		self.read = read
		self.buffered_repeat = self.read.seq[self.start_buffered:(self.end+b)]

	def qc(self):
		""" assigns a quality control; if no flag is raised, indicates that the repeat is contained within the entire fragment within a 3-4 nt bound"""
		qc = ''
		if self.start - 3 < 0:
			qc += '^'
		if self.read.len < self.end + 4:
			qc += '*'
		return qc

	def get_count(self, probe):
		""" get the repeat counts -- note that certain if exceptions may need to be made for complex probes""" 
		if probe.id == "D21S11":
			self.count = (self.match_len-11)/4
		else:
			self.count = self.match_len/probe.len
			self.remainder = self.match_len%probe.len
			if self.remainder > 0:
				self.count = str(self.count)+"."+str(self.remainder)
				self.count = Decimal(self.count)
		return self.count

	def get_chr_pos(self, probe):
		""" get the chromosomal position of the repeat sequence, 1-indexed, not python 0-indexing"""
		self.read_start = self.read.pos + self.start + 1
		self.read_end = self.read_start + self.end
		self.region = "chr"+probe.chr+" : "+str(self.read_start)+" - "+str(self.read_end)
		return self.region

class Read():
	def __init__(self, list):
		self.qname = list[0]
		self.flag = list[1]
		self.seq = list[2]
		self.pos = list[3]
		self.len = len(self.seq)
	
	def PerfectRegex(self, probe, read, options):
		""" return a list of all matches in the read using perfect regular expressions"""
		list = [Match(m, read, options.buffer) for m in re.finditer(probe.regex(), self.seq)]
		return list
	
	def FuzzyRegex(self, probe):
		matches = re.finditer(probe.regex(), self.seq)
		return matches
		
class Probe():
	def __init__(self, list, options):
		"""parse probe config file."""
		self.id = list[0]
		self.chr = str(list[1])
		self.start = int(list[2])
		self.end = int(list[3])
		try:
			self.repeat = list[4]
			self.len = len(self.repeat)
			self.complseq = complement(self.repeat)
		except IndexError:
			# this addresses an issue with blank repeat sequences in probe list
			self.repeat = ''
			self.len = 4
		if self.len == 0:
			self.len = 4

		try:
			self.cutoff = int(list[5])
		except IndexError:
			# default cut-off
			self.cutoff = options.min_cutoff
		except ValueError:
			self.cutoff = options.min_cutoff
			
	def regex(self):
		""" get the associated regular expression for the probe's repeat sequence"""
		regexpr = get_regexpr(self.id, self.repeat)
		return regexpr

def get_options():
	# set up command line options
	usage = "python str_finder.py [options] -o strs -c hg19 [args] bam_list.txt"

	p = optparse.OptionParser(usage)
	p.add_option('--output', '-o', default="result", help="designate output file name")
	p.add_option('--buffer', '-b', default=6, type = "int", help="display b nucleotides before and after the repeat sequence")
	p.add_option('--top_n', '-n', default=2, type = "int" , help="return top n allele counts")
	p.add_option('--min_cutoff', '-m', default=4, type = "int", help="Minimum threshold for reported repeat counts if not given in the config file.")
	p.add_option('--probe', '-p', default="codis_hg19.config", help="Probe config file - provides chromosomal location data and repeat sequence.")
	options, args = p.parse_args()
	return options, args

def clear_files(data_file, summary_file):
	f = open(data_file, "w")
	f.write("FILENAME	QUERY ID	PROBE ID	ALLELE	READ LENGTH CHR LOCATION	REPEAT SEQ \n")
	f.close()
	g = open(summary_file, "w")
	g.write("FILENAME	PROBE ID	GENOTYPE	# READS	# MATCHED READS	COVERAGE \n")
	g.close	
		
def get_probe_list(options):	
	"""get probe list from config file"""
	try:
		probe_list = parse_file(options.probe)
	except IOError:
		print 'Config file not found! Try again'
		probe_list = get_file()
	return probe_list
	
def get_input_list(args):
	try:
		input_list = flatten(parse_file(args[0]))
	except IndexError:
		print 'No input list given. Input list file?'
		input_list = flatten(get_file())
	except IOError:
		print 'Input file not found'
		input_list = flatten(get_file())
	return input_list

def check_inputs(file):
	"""check that input files exist and identify appropriate naming convention for filtering reads by chromosome"""
	try:
		samfile = pysam.Samfile(file, "rb")
	except IOError:
		print 'file not found'
		pass
		
	try:
		pysam.view("-X", file, "chr11:1000-1100")
		head = "chr"
	except NameError:
		head = ""
	
	return samfile, head
	
def hg_inputs(file):
	"""for handling 1000 Genomes Project alignment indices as input"""

	samfile = pysam.Samfile('ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/'+file, 'rb')
	head = ""

	return samfile, head
	
def parse_file(file, flag = 0):
	"""Parses a tab-delimited file into a list. If flag is 1, parses file to dictionary"""
	sample = open(file, "rb")
	readSample = csv.reader(sample, dialect="excel-tab")
	sampleList = []
	sampleDict = {}
	
	if flag == 0:
		for row in readSample:
			sampleList.append(row)
		return sampleList
	if flag == 1:
		for row in readSample:
			sampleDict[row[0]] = row[1]
		return sampleDict
	sample.close()

def complement(seq):  
	"""returns the reversed complementary sequence"""
	complseq = seq.translate(string.maketrans("ATCG", "TAGC"))[::-1]
	return complseq
	
def get_regexpr(probe, repeat):
	"""Generate a regular expression to search for the repeat string. Complex repeats are stored in complex.config"""
	comp_seq = parse_file('complex.config', 1)
	try:
		regex = comp_seq[probe]
	except KeyError:	
		regex = "("+repeat+")+.?.?.?("+repeat+")+|("+complement(repeat)+")+.?.?.?("+complement(repeat)+")+"
	#print probe, regex
	return regex

def flatten(x):
    """http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks -- reduces levels of nesting"""
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def get_file():
	i = raw_input('')
	try:
		list = parse_file(i)
	except IOError:
		print 'File not found. Try again'
		list = get_file()
	return list

def main():			
	options, args = get_options()
	
	b = options.buffer
	n = options.top_n
	
	data_file = options.output+'_data.out'
	summary_file = options.output+'_summary.out'

	clear_files(data_file, summary_file)
	
	probe_list = get_probe_list(options)
	input_list = get_input_list(args)
	
	# LOOP STARTS HERE
	for STRlocus in probe_list:
		probe = Probe(STRlocus, options)
		new_list = []
		for file in input_list:
			expression = "(chrom)"+str(probe.chr)+"\."
			if re.search(expression, file):
				new_list.append(file)
		k = 1	
		for file in new_list:
			#print file
			samfile, head = hg_inputs(file)

			# FILTER READS
			# fetch all reads in designated region and add them to a list. An additional QC control can be implemented here by filtering the query quality measure (AlignedRead.qual)
			read_list = []		
			for read in samfile.fetch(head+probe.chr, probe.start, probe.end):
				read_list.append([read.qname, read.flag, read.seq, read.pos])

				# FIND REPEATS
				count_hist = {}
				frag_list = []
				for r in read_list:
					read = Read(r)
					# returns all the matches contained within one read 
					match_per_read = read.PerfectRegex(probe, read, options)
					
					# If a match is found, return the longest match.
					if len(match_per_read) > 0:
						list = [i.match_len for i in match_per_read]
						match = match_per_read[list.index(max(list))]
						# if no qc flags are raised, add read to list of all matches and get its count.
						if match.qc() == '':
							frag_list.append(match)
							count = match.get_count(probe)				
						else:
							count = 0
						# INFER GENOTYPE
						# produce a histogram of the counts
						if count >= probe.cutoff:
							if count in count_hist:
								count_hist[count] += 1
							else:
								count_hist[count] = 1
				# COVERAGE INFO
				total_reads = len(read_list)
				total_matches = len(frag_list)
				getcontext().prec = 2
				ratio = Decimal(len(frag_list))/Decimal(len(read_list))

				#potential filter mechanism -- only return items that have frequency_count > 1 (good coverage), unless ...
				"""
				del_list = []
				if total_reads < 125 and ratio > 0.15:
					count_hist = count_hist
				else:
					del_list = [item for item in count_hist if count_hist[item] > 1]
					for item in del_list:
						count_hist.pop(item)
				
				print len(read_list), len(frag_list), Decimal(len(frag_list))/Decimal(len(read_list))
				print count_hist
				"""
			
				if len(count_hist):
					flag = 0
					# Sort by frequency, then by number of repeats, and return the top n. A more sophisticated peak finder may be eventually necessary.
					list = [sorted(count_hist.items(), key=itemgetter(1, 0), reverse=True)]
					result = Result(list, n)
					counts = result.counts	
					i = 0
					
					f = open(data_file, "a")
					for match in frag_list:
						# print all the non-qc-flagged matches to data file, regardless of count
						data_row = [file, read.qname, probe.id, match.count, match.read.len, match.get_chr_pos(probe), match.buffered_repeat]
						writeSample = csv.writer(f, lineterminator="\n", delimiter="\t")
						writeSample.writerow(data_row)
					i += 1
					f.close()			
				
				else:
					# if no matches are found, return no counts in the summary file
					list = []
					result = Result(list, n)
					flag = 1
				
			summary = [file, probe.id, result.summary(flag), total_reads, total_matches, ratio]
			f = open(summary_file, "a")
			writeSample = csv.writer(f, lineterminator="\n", delimiter="\t")
			writeSample.writerow(summary)		
			f.close()	
			
			# print progress
			print probe.id+": "+str(k)+"/"+str(len(new_list))
			k += 1
			samfile.close()
		
if __name__ == '__main__':
	main()

# python strider_1khg.py -p powerplex16_hg19.config align_filtered.out
