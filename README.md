NG-STRIder
Last updated: August 4, 2010

=== QUICK START ===
Requires Python 2.5 or 2.6 and a Linux OS.
Dependencies: zlib, ncurses, Pyrex, Pysam.
For use with remote access files, see this Pysam update: 
http://code.google.com/p/pysam/source/detail?r=0ebb13e040

Usage details are on the project page:
http://code.google.com/p/strfinder

=== EXAMPLE CALLS ===
Get STR profiles for genomes on the Sanger FTP (listed in bam_list.txt)
   > python strider_main.py -o sanger bam_list.txt
   
	Returns sanger_data.out and sanger_summary.out
	
Get STR profiles for all Han Chinese in Beijing from 1000 Genomes Project Data
	First, get a filtered list of the BAM file locations.
	> python snippets/filter_align_index.py -p CHB -o 1khg_CHB
	
	Then, run the NG-STRIder (for 1KHG) on this alignment list
	> python strider_1khg -o profiles_chb 1khg_CHB.out 
	
	Returns profiles_chb_data.out and profiles_chb_summary.out
	
Get frequency distribution alleles for the CHB data
	> python calc_distribution.py -o CHB profiles_chb_summary.out
	
	Returns CHB_freqs.txt
