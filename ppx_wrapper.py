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

	def get_region(self, start, stop, strand):
		region = self.seq[start:stop]
		if strand == '-':
			complement = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
			region = "".join([complement.get(nt.upper(), '') for nt in region[::-1]])
		return region

def parse_contigs_to_dict(contig_file):
	contigs = set()
	header, seq = '', ''
	print "[STATUS] - Parsing contigs "
	with open(contig_file) as fh:
		for line in fh:
			if line.startswith(">"):
				if (seq):
					if len(seq) > 1000: 
						header = species + "." + header
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
	temp_file = TEMP_DIR + header + ".temp"
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

	def get(self, arg, buffer_range):
		if arg == "start":
			coordinates = self.coordinates[:]
			sort(coordinates)
			if coordinates[0] <= int(buffer_range):
				return 0
			else:
				return (self.coordinates[0] - int(buffer_range)) 
		elif arg == "end":
			coordinates = self.coordinates[:]
			sort(coordinates)
			return (self.coordinates[-1] + int(buffer_range))
		elif arg == 'strand':
			if self.strand == '+':
				return 'forward'
			elif self.strand == '-':
				return 'backward'
			else:
				sys.exit("[ERROR] - Strand : " + self.strand)
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
		print score, dict_of_blocks[score].__dict__
	for score in sorted(dict_of_blocks, reverse=True):
		block = dict_of_blocks[score]
		contig = block.contig

		start = block.get('start', 10000)
		end = block.get('end', 10000)
		strand = block.get('strand', 0)

		return block.contig, str(start), str(end), strand, str(score) 
		#break


def runAugustusPPX():
	dict_of_blocks = {}
	for result in os.listdir("fastblocksearch/"):
		# For each FastBlockSearch result file ...
		result_file = "fastblocksearch/" + result
		if result.startswith(species) and result_file.endswith(".result"):
			#temp_file = "temp/" + result_file.replace(".result", ".temp")
			# Get results
			list_of_blocks = parseFastBlockSearchResult(result_file)
			if len(list_of_blocks) == 0:
				pass
			else:
				for block in list_of_blocks:
					dict_of_blocks[block.score] = block

	contig, start, end, strand, score = selectBestBlock(dict_of_blocks)
	infile = TEMP_DIR + contig + ".temp"
	outfile = AUGUSTUS_DIR + contig + "." + query + ".gff3"
	print "[STATUS] - Calling protein \"" + query + "\" in contig \"" + contig + "\""
	process = subprocess.Popen("/exports/software/augustus/augustus-3.0.3/bin/augustus --species=caenorhabditis --gff3=on --proteinprofile=" + profile + " --predictionStart=" + start + " --predictionEnd=" + end + " --strand=" + strand + " " + infile + " > " + outfile , stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
	process.wait()
	print "[STATUS] - Writing proteins"
	proteins = parseProteinsFromGFF3(outfile)
	for protein_name, protein_seq in proteins.items():
		for motif in MOTIFS:
			if motif in protein_seq:
				print ">" + protein_name + "\n" + protein_seq
				break

	print "[STATUS] - Done."

def parseProteinsFromGFF3(gff3):
	#print ".".join(gff3.split(".")[0:-1])) + "aa.fa"
	#print gff3.split("/")[1].split(".")[1]
	#print gff3.split("/")[1].split(".")[0]
	contig, query, outfile = gff3.split("/")[1].split(".")[1], gff3.split("/")[1].split(".")[0], ".".join(gff3.split(".")[0:-1]) + "aa.fa"

	#fh = open(gff3 + ".proteins.fa", 'w')

	read_mode = 0
	temp = ''
	protein_name = ''
	protein_seq = ''
	dict_of_proteins = {}

	with open(gff3) as fh:
		for line in fh:
			if line.startswith("# start gene "):
				read_mode = 1
				protein_name = line.split(" ")[3].rstrip("\n")
			elif line.startswith("# end gene "):

				dict_of_proteins[protein_name] = protein_seq

				read_mode = 0
				protein_name = ''
				protein_seq = ''
			elif read_mode == 1 and line.startswith("#"):
				temp = line.lstrip("# protein sequence = [")
				temp = temp.lstrip("# ")
				temp = temp.rstrip("\n")
				temp = temp.rstrip("]")
				protein_seq += temp
			else:
				pass
	return dict_of_proteins

if __name__ == "__main__":
	try:
		contig_file = sys.argv[1]
		profile = sys.argv[2]
		species = sys.argv[3] # ID for contigs, etc
	except:
		sys.exit("Usage: ./ppx_wrapper.py [CONTIGFILE] [PROFILE]")
	
	query = profile.split("/")[-1].split(".")[0]
	
	GENOME_DIR = 'genome/'
	TEMP_DIR = 'temp/'
	FASTBLOCKSEARCH_DIR = 'fastblocksearch/'
	AUGUSTUS_DIR = 'augustus/'
	MOTIFS = ["ELEKEF", "WFQNRR"]
	contigs = parse_contigs_to_dict(contig_file)
	#print contigs
	fastblocksearch(profile, contigs)
	runAugustusPPX()