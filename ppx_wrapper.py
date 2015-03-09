#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File    	: ppx_wrapper.py
Author  	: 
Version 	: 0.1
Description : 
To do 		: Make global scoring for fastblocks for all profiles ...
"""

from __future__ import division
import sys, os
import time
import multiprocessing as mp 
import subprocess
import string
import collections

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

class Block():
	def __init__(self, contig, score, multi_score):
		self.contig = contig
		self.score = float(score)
		self.multi_score = float(multi_score)
		self.coordinates = []
		self.strand = ''
		self.profile = ''

	def get(self, arg, buffer_range):
		if arg == "start":
			coordinates = self.coordinates[:]
			coordinates = sorted(coordinates)
			if coordinates[0] <= int(buffer_range):
				return 0
			else:
				return (coordinates[0] - int(buffer_range)) 
		elif arg == "end":
			coordinates = self.coordinates[:]
			coordinates = sorted(coordinates)
			return (coordinates[-1] + int(buffer_range))
		elif arg == 'strand':
			if self.strand == '+':
				return 'forward'
			elif self.strand == '-':
				return 'backward'
			else:
				sys.exit("[ERROR] - Strand : " + self.strand)
		else:
			sys.exit("postion... What?")

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
						contig = ContigObject(header, seq.upper())
						contigs.add(contig) 
				seq = ''
				header = line.rstrip("\n").lstrip(">").replace(" ","_")	
			else:
 				#seq += line.translate(None,string.ascii_lowercase).rstrip("\n")
 				seq += line.rstrip("\n")
		if len(seq) > 1000: 
			header = species + "." + header
			contig = ContigObject(header, seq)
			contigs.add(contig) 
	return contigs

def align_fasta(contig_file):
	pass

def make_msa_profile(contig_file):
	pass

def progress(counter, max_value):
	sys.stdout.write('\r')
	progress = int(counter)/int(max_value)
	print "\tProgress : " + format(float(progress),'.2%'),
	sys.stdout.flush()

def fastblocksearch(profile, contigs):
	pool = mp.Pool(processes = 10)
	counter, max_value = 0, len(contigs)
	profile_name = profile.split("/")[-1].split(".")[0]
	print "[STATUS] - Running FastBlockSearch with profile : " + profile_name
	for contig in contigs:
		counter += 1
		progress(counter, max_value)
		pool.apply(run_fastblocksearch, args=(profile, contig,))
	sys.stdout.write('\r')
	print "\tProgress : 100.00%"

def run_fastblocksearch(profile, contig):
	out_file = FASTBLOCKSEARCH_DIR + contig.header + "." + profile.split("/")[-1].split(".")[0] + ".result"
	temp_file = TEMP_DIR + contig.header + ".temp"
	temp = open(temp_file, 'w')
	temp.write(">" + contig.header + "\n" + contig.seq)
	temp.close()
	process = subprocess.Popen("/exports/software/augustus/augustus-3.0.3/bin/fastBlockSearch --cutoff=0.5 " + temp_file + " " + profile + " > " + out_file + " ", stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)


''' This has to be done for each contig in each species (by species)
	- create folder for each genome
	- infer best overall hits	
'''

def parseFastBlockSearchResult(results):
	list_of_blocks = []
	contig, profile = results.split("/")[1].split(".")[1], results.split("/")[1].split(".")[2]
	raw = open(results).read()
	if len(raw.split("\n")) > 4:
		blocks = [filter(None, x.split('\n')) for x in raw.split("--") if len(x) > 2 ] 
		header = blocks[0].pop(0)
		#print blocks
		for hit in blocks:
			#print "New hit \n" + str(hit)
			score, multi_score, coordinate, strand = 0.0, 0.0, 0, ''
			score = float(hit[0].lstrip("Score:"))
			multi_score = float(hit[1].lstrip("Mult. score:"))
			block = Block(contig, score, multi_score)
			block.coordinates = sorted([int(x.split("\t")[0]) for x in hit[2:] if len(x.split("\t")) > 1])
			block.strand = hit[-2].split("\t")[2]
			block.profile = profile
		list_of_blocks.append(block)
	else:
		pass
	return list_of_blocks

def analyseBlocks(dict_of_blocks, dict_of_contigs):

	for profile in sorted(dict_of_blocks, reverse=True):
		for score in sorted(dict_of_blocks[profile], reverse=True):
			block = dict_of_blocks[profile][score]
			contig = block.contig
			start = block.get('start', 10000)
			end = block.get('end', 10000)
			strand = block.get('strand', 0)
			profile = block.profile

	for contig in dict_of_contigs:
		print contig 
		for profile in dict_of_contigs[contig]:
			print "\t" + profile,
			for score in dict_of_contigs[contig][profile]:
				print str(score) + str(dict_of_contigs[contig][profile][score].__dict__)

	#return block.contig, str(start), str(end), strand, str(score), profile 		#break

def runAugustusPPX():
	dict_of_blocks = AutoVivification()
	dict_of_contigs = AutoVivification()

	for result in os.listdir("fastblocksearch/"):
		# For each FastBlockSearch result file ...
		if result.startswith(species) and result.endswith(".result"):
			result_file = "fastblocksearch/" + result
			list_of_blocks = parseFastBlockSearchResult(result_file)
			if (list_of_blocks):
				for block in list_of_blocks:
					dict_of_blocks[block.profile][block.score] = block
					dict_of_contigs[block.contig][block.profile][block.score] = block

	analyseBlocks(dict_of_contigs)

	infile = TEMP_DIR + contig + ".temp"
	outfile = AUGUSTUS_DIR + contig + "." + profile + ".gff3"
	print "[STATUS] - Calling protein \"" + profile + "\" in contig \"" + contig + "\" from " + str(start) + " to " + str(end)  
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

def getProfiles(profile_dir):
	list_of_profiles = []
	for profile in os.listdir(profile_dir + "/"):
		if profile.endswith(".prfl"):
			list_of_profiles.append(profile_dir + profile)
	return list_of_profiles

if __name__ == "__main__":
	try:
		contig_file = sys.argv[1]
		#profile = sys.argv[2]
		profile_dir = sys.argv[2]
		species = sys.argv[3] # ID for contigs, etc
	except:
		sys.exit("Usage: ./ppx_wrapper.py [CONTIGFILE] [PROFILE_DIR] [SPECIES]")
	
	GENOME_DIR = 'genome/'
	TEMP_DIR = 'temp/'
	FASTBLOCKSEARCH_DIR = 'fastblocksearch/'
	AUGUSTUS_DIR = 'augustus/'
	MOTIFS = ["ELEKEF", "WFQNRR"]
	RESULTS_DIR = 'results/'

	list_of_profiles = getProfiles(profile_dir)
	contigs = parse_contigs_to_dict(contig_file)
	#print contigs
	for profile in list_of_profiles:
		#fastblocksearch(profile, contigs)
		runAugustusPPX()