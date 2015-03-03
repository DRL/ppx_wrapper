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

def parse_contigs_to_dict(contig_file):
	contigs = {}
	header, seq = '', ''
	with open(infile) as fh:
		for line in fh:
			if line.startswith(">"):
				if (seq):
					contigs[header] = seq
					seq = ''
				header = line.rstrip("\n").lstrip(">")	
			else:
				seq += line.rstrip("\n")
		contigs[header] = seq
	return contigs

def align_fasta(infile):
	pass

def make_msa_profile(aln_file):
	pass

def run_fastblocksearch(profile, contig):
	subprocess.Popen("/exports/software/augustus/augustus-3.0.3/bin/fastBlockSearch --cutoff=0.5 " + contig + " " + profile, stdout=subprocess.PIPE, shell=True)

def fastblocksearch(profile_file, contigs):
	pool = mp.Pool(processes=4)
	results = [pool.apply(run_fastblocksearch, args=(profile_file, contig)) for contig in contigs]
	print(results)

if __name__ == "__main__":
	try:
		contig_file = sys.argv[1]
		profile_file = sys.argv[2]
	except:
		sys.exit("Usage: ./fastats.py [FASTAFILE]")
	
	contigs = parse_contigs_to_dict(contig_file, profile_file)
	fastblocksearch(contigs)