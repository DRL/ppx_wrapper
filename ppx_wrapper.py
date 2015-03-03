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

class Contig():
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
					print "Parsed " + header
					contigs.add(Contig(header, seq)) 
					header = line.rstrip("\n").lstrip(">")	
					seq = ''
			else:
 				seq += line.translate(None,string.ascii_lowercase).rstrip("\n")
		print "Parsed " + header
		contigs.add(Contig(header, seq)) 
	return contigs

def align_fasta(contig_file):
	pass

def make_msa_profile(contig_file):
	pass

def run_fastblocksearch(profile, contig):
	print "Searching " + profile + " in " + contig.header
	temp_file = open(header + ".temp", 'w')
	temp_file.write(">" + contig.header + "\n" + contig.seq)
	temp_file.close()
	subprocess.Popen("/exports/software/augustus/augustus-3.0.3/bin/fastBlockSearch --cutoff=0.5 " + temp_file + " " + profile + " > " + contig.header + ".result"  , stdout=subprocess.PIPE, shell=True)

def fastblocksearch(profile_file, contigs):
	pool = mp.Pool(processes=4)
	results = [pool.apply(run_fastblocksearch, args=(profile_file, contig)) for contig in contigs]
	print(results)

if __name__ == "__main__":
	try:
		contig_file = sys.argv[1]
		profile_file = sys.argv[2]
	except:
		sys.exit("Usage: ./ppx_wrapper.py [CONTIGFILE] [PROFILE]")
	
	contigs = parse_contigs_to_dict(contig_file)
	print contigs
	fastblocksearch(contigs, profile_file)