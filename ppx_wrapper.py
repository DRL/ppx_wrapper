#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File    	: ppx_wrapper.py
Author  	: 
Version 	: 0.1
Description : 
To do 		: 
"""

from __future__ import division
import sys
import time
import multiprocessing as mp 
import subprocess
import string

class ContigObject():
	def __init__(self, header, seq):
		self.header = header
		self.seq = seq

def parse_contigs_to_dict(contig_file):
	contigs = set()
	header, seq = '', ''
	with open(contig_file) as fh:
		for line in fh:
			if line.startswith(">"):
				if (seq):
					if len(seq) > 100000: 
						print "Parsed " + header
						contig = ContigObject(header, seq)
						contigs.add(contig) 
				seq = ''
				header = line.rstrip("\n").lstrip(">").replace(" ","_")	
			else:
 				#seq += line.translate(None,string.ascii_lowercase).rstrip("\n")
 				seq += line.rstrip("\n")
		if len(seq) > 100000: 
			print "Parsed " + header
			contig = ContigObject(header, seq)
			contigs.add(contig) 
	return contigs

def align_fasta(contig_file):
	pass

def make_msa_profile(contig_file):
	pass

def run_fastblocksearch(profile, header, seq):
	print "Searching " + profile + " in " + header

	temp_file = header + ".temp"
	temp = open(temp_file, 'w')
	temp.write(">" + header + "\n" + seq)
	temp.close()
	process = subprocess.Popen("/exports/software/augustus/augustus-3.0.3/bin/fastBlockSearch --cutoff=0.5 " + temp_file + " " + profile + " > " + header + ".result  &>log.txt", stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
	output, error = process.communicate()
	print output
def fastblocksearch(contigs, profile_file,):
	print str(len(contigs)) + " contigs"
	pool = mp.Pool(processes=40)
	results = [pool.apply(run_fastblocksearch, args=(profile_file, contig.header, contig.seq)) for contig in contigs]
	#output = [p.get() for p in results]
	#print(output)
	print "Done"

if __name__ == "__main__":
	try:
		contig_file = sys.argv[1]
		profile_file = sys.argv[2]
	except:
		sys.exit("Usage: ./ppx_wrapper.py [CONTIGFILE] [PROFILE]")
	
	contigs = parse_contigs_to_dict(contig_file)
	#print contigs
	fastblocksearch(contigs, profile_file)