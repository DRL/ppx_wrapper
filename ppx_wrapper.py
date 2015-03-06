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
import sys, os
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
	print "[STATUS] - Parsing contigs "
	with open(contig_file) as fh:
		for line in fh:
			if line.startswith(">"):
				if (seq):
					if len(seq) > 100000: 
						contig = ContigObject(header, seq)
						contigs.add(contig) 
				seq = ''
				header = line.rstrip("\n").lstrip(">").replace(" ","_")	
			else:
 				#seq += line.translate(None,string.ascii_lowercase).rstrip("\n")
 				seq += line.rstrip("\n")
		if len(seq) > 100000: 
			contig = ContigObject(header, seq)
			contigs.add(contig) 
	return contigs

def align_fasta(contig_file):
	pass

def make_msa_profile(contig_file):
	pass

def run_fastblocksearch(profile, header, seq):
	temp_file = ''
	temp_file = header + ".temp"
	temp = open(temp_file, 'w')
	temp.write(">" + header + "\n" + seq)
	temp.close()
	process = subprocess.Popen("/exports/software/augustus/augustus-3.0.3/bin/fastBlockSearch --cutoff=0.5 " + temp_file + " " + profile + " > fastblocksearch/" + header + ".result ", stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

def progress(counter, max_value):
	sys.stdout.write('\r')
	progress = int(counter)/int(max_value)
	print "\tProgress : " + format(float(progress),'.2%'),
	sys.stdout.flush()

def fastblocksearch(profile, contigs):
	print "[STATUS] - Running FastBlockSearch "
	pool = mp.Pool(processes = 10)
	counter, max_value = 0, len(contigs)
	for contig in contigs:
		counter += 1
		if counter % 5 == 0:
			progress(counter, max_value)
		pool.apply(run_fastblocksearch, args=(profile, contig.header, contig.seq,))
	sys.stdout.write('\r')
	print "\tProgress : 100.00%"


class Block():
	def __init__(self, contig, score, multi_score):
		self.contig = contig
		self.score = float(score)
		self.multi_score = float(multi_score)
		self.coordinates = []
		self.strand = ''

	def add_coordinate(self, coordinate):
		self.coordinates.append(coordinate)

	def get(self, position):
		if position == "start":
			return self.coordinates[0]
		elif position == "end":
			return self.coordinates[-1]
		else:
			sys.exit("postion... What?")

''' This has to be done for each contig in each species (by species)
	- create folder for each genome
	- infer best overall hits	
'''

def parseFastBlockSearchResult(results):
	list_of_blocks = []
	with open(results) as fh:
		contig, score, multi_score, coordinate, strand = '', 0.0, 0.0, 0, ''
		for line in fh:
			line = line.rstrip("\n")
			if line.startswith("Hits found in "):
				contig = line.lstrip("Hits found in ")
			elif line.startswith("Score:"):
				score = float(line.lstrip("Score:"))
			elif line.startswith("Mult. score:"):
				multi_score = float(line.lstrip("Mult. score:"))
				block = Block(contig, score, multi_score)
			elif line.startswith("--"):
				list_of_blocks.append(block)
			elif not line:
				pass
			else:
				coordinate, strand = line.split("\t")[0], line.split("\t")[2]
				block.add_coordinate(int(coordinate))
				block.strand = strand
	return list_of_blocks

def selectBestBlock(dict_of_blocks):
	for score in sorted(dict_of_blocks, reverse=True):
		block = dict_of_blocks[score]
		contig = block.contig
		start = block.get('start')
		end = block.get('end')
		strand = block.strand
		yield str(score) + " " + block.contig + " " + str(start) + " " + str(end) + " " + strand 
		#break


def runAugustusPPX():
	dict_of_blocks = {}
	for result in os.listdir("fastblocksearch/"):
		# For each FastBlockSearch result file ...
		result_file = "fastblocksearch/" + result
		if result_file.endswith(".result"):
			temp_file = result_file.replace(".result", ".temp")
			# Get results
			list_of_blocks = parseFastBlockSearchResult(result_file)
			if len(list_of_blocks) == 0:
				pass
			else:
				for block in list_of_blocks:
					dict_of_blocks[block.score] = block

	print selectBestBlock(dict_of_blocks)
	
if __name__ == "__main__":
	try:
		contig_file = sys.argv[1]
		profile = sys.argv[2]
	except:
		sys.exit("Usage: ./ppx_wrapper.py [CONTIGFILE] [PROFILE]")
	
	FASTBLOCKSEARCH_DIR = 'fastblocksearch/'
	AUGUSTUS_DIR = 'augustus/'

	contigs = parse_contigs_to_dict(contig_file)
	#print contigs
	fastblocksearch(profile, contigs)
	runAugustusPPX()