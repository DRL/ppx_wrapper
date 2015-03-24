#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File    	: ppx_wrapper.py
Author  	: 
Version 	: 0.1
Description : 
To do 		:  ...
Bugs 		: Contig names in assemblies can't have "."
"""

from __future__ import division
import sys, os
import time
import multiprocessing as mp 
import subprocess
import string
import collections
from subprocess import Popen

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
	
	temp_dict = {}

	for contig in contigs:
		temp_file = TEMP_DIR + contig.header + ".temp"
		temp = open(temp_file, 'w')
		temp.write(">" + contig.header + "\n" + contig.seq)
		temp.close()
		temp_dict[contig.header]=temp_file
	print "[STATUS] - %s contigs found " %len(contigs)
	return temp_dict

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
	#pool = mp.Pool(processes = 10)
	counter, max_value = 0, len(contigs)
	profile_name = profile.split("/")[-1].split(".")[0]
	print "[STATUS] - Running FastBlockSearch with profile : " + profile_name
	processes = []

	for contig in contigs:
		time.sleep(0.1)
		temp_file = contigs[contig]
		out_file = FASTBLOCKSEARCH_DIR + contig + "." + profile.split("/")[-1].split(".")[0] + ".result"
		process = subprocess.Popen("/exports/software/augustus/augustus-3.0.3/bin/fastBlockSearch --cutoff=0.5 " + temp_file + " " + profile + " > " + out_file + " ", stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
		processes.append(process)
		progress(counter, max_value)
		counter += 1
		if counter % 100 == 0:
			for process in processes:
				process.wait()
				time.sleep(1)
	#for contig in contigs:
	#	counter += 1
	#	progress(counter, max_value)
	#	jobs = pool.apply(run_fastblocksearch, args=(profile, contig,))
	#	time.sleep(1)
	#pool.close()
	sys.stdout.write('\r')
	print "\tProgress : 100.00%"

#def run_fastblocksearch(profile, contig):
#	out_file = FASTBLOCKSEARCH_DIR + contig.header + "." + profile.split("/")[-1].split(".")[0] + ".result"
#	temp_file = TEMP_DIR + contig.header + ".temp"
#	temp = open(temp_file, 'w')
#	temp.write(">" + contig.header + "\n" + contig.seq)
#	temp.close()
#	process = subprocess.Popen("/exports/software/augustus/augustus-3.0.3/bin/fastBlockSearch --cutoff=0.5 " + temp_file + " " + profile + " > " + out_file + " ", stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
#	#process.wait()

''' This has to be done for each contig in each species (by species)
	- create folder for each genome
	- infer best overall hits	
'''

def parseFastBlockSearchResult(result_file):
	list_of_blocks = []
	contig, profile = result_file.split("/")[1].split(".")[1], result_file.split("/")[1].split(".")[2]
	raw = open(result_file).read()
	try:
		number_of_hits = len(raw.split("--"))
	except:
		number_of_hits = 0
	if number_of_hits >= 2:
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
			block.strand = hit[-1].split("\t")[2]
			block.profile = profile
		list_of_blocks.append(block)
	else:
		#print raw
		os.remove(result_file) 
	return list_of_blocks

def analyseBlocks(dict_of_blocks):


	# check those that overlap -> yes/no
	fastblockresults_dict = {}
	max_profile_count = 1
	profile_count = {}

	profile_hits = {}


	for score in sorted(dict_of_blocks, reverse=True):
		print str(score) + "\t" + str(dict_of_blocks[score].__dict__)
		block = dict_of_blocks[score]
		contig = block.contig
		profile = block.profile
		if not contig in fastblockresults_dict:
			profile_count[profile] = profile_count.get(profile, 0) + 1
			if profile_count[profile] <= max_profile_count:
				fastblockresults_dict[contig] = []
				fastblockresults_dict[contig].append(block)
				if not profile in profile_hits:
					profile_hits[profile]=[]
					profile_hits[profile].append(block)
				else: 
					profile_hits[profile].append(block)
		else: 
			for hit in fastblockresults_dict[contig]:
				hit_start, hit_end = int(hit.get('start', buffer_range)), int(hit.get('end', buffer_range))
				block_start, block_end = int(block.get('start', buffer_range)), int(block.get('end', buffer_range))
				coordinates = [hit_start, hit_end, block_start, block_end]
				sum_lengths = (hit_end - hit_start) + (block_end - block_start)
				if sum_lengths >= (max(coordinates) - min(coordinates)):
					# "Overlap"
					pass
				else:
					# "No Overlap"
					if (coordinates[0] >= coordinates[2] and coordinates[1] <= coordinates[3]):
						# hit in block
						pass
					elif (coordinates[2] >= coordinates[0] and coordinates[3] <= coordinates[1]):
						# block in hit
						pass
					else:
						profile_count[profile] = profile_count.get(profile, 0) + 1
						if profile_count[profile] <= max_profile_count:
							fastblockresults_dict[contig].append(block)
							if not profile in profile_hits:
								profile_hits[profile]=[]
								profile_hits[profile].append(block)
							else: 
								profile_hits[profile].append(block)

				#if not block.contig in fastblockresults_dict:
				#	fastblockresults_dict[contig] = block
				#else:
				#	pass
				#	#if fastblockresults_dict[contig].start > block.start and fastblockresults_dict[contig].end < block.start
				#block = dict_of_blocks[profile][score]
				#contig = block.contig
				#start = block.get('start', buffer_range)
				#end = block.get('end', buffer_range)
				#strand = block.get('strand', 0)
				#profile = block.profile

			#if not profile in fastblockresults_dict:
			#	fastblockresults_dict[profile] = block
			#else:
			#	if score > fastblockresults_dict[profile].score:
			#		fastblockresults_dict[profile] = block
#
				
			

	#for contig in dict_of_contigs:
	#	print contig 
	#	for profile in dict_of_contigs[contig]:
	#		print "\t" + profile,
	#		for score in dict_of_contigs[contig][profile]:
	#			print str(score) + str(dict_of_contigs[contig][profile][score].__dict__)
	for profile in profile_hits:
		print profile
		for hits in profile_hits[profile]:
			print hits.__dict__
	
	return profile_hits

def runAugustusPPX():
	dict_of_blocks = {}
	print "Buffer range : " + str(buffer_range) 
	for result in os.listdir("fastblocksearch/"):
		# For each FastBlockSearch result file ...
		if result.startswith(species) and result.endswith(".result"):
			result_file = "fastblocksearch/" + result
			print "[STATUS] - Parsing : " + result_file
			list_of_blocks = parseFastBlockSearchResult(result_file)
			if (list_of_blocks):
				for block in list_of_blocks:
					dict_of_blocks[block.score] = block

					#dict_of_contigs[block.contig][block.profile][block.score] = block

	profile_hits = analyseBlocks(dict_of_blocks)

	ppx_results = open(RESULTS_DIR + species + ".fa", "w")
	for profile in profile_hits:
		for hit in profile_hits[profile]:
			contig = hit.contig
			start = str(hit.get('start', buffer_range))
			end = str(hit.get('end', buffer_range))
			strand = hit.get('strand', 0)
			
			infile = TEMP_DIR + species + "." + contig + ".temp"
			outfile = AUGUSTUS_DIR + species + "." + contig + "." + profile + ".gff3"
			profile_file = dict_of_profiles[profile]
			print "[STATUS] - Calling protein \"" + profile_file + "\" in contig \"" + contig + "\" from " + str(start) + " to " + str(end)  
			#process = subprocess.Popen("/exports/software/augustus/augustus-3.0.3/bin/augustus --species=caenorhabditis --gff3=on --proteinprofile=" + profile + " --predictionStart=" + start + " --predictionEnd=" + end + " --strand=" + strand + " " + infile + " > " + outfile , stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
			process = subprocess.Popen("/exports/software/augustus/augustus-3.0.3/bin/augustus --species=caenorhabditis --gff3=on --proteinprofile=" + profile_file + " --predictionStart=" + start + " --predictionEnd=" + end + " --strand=" + strand + " " + infile + " > " + outfile , stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
			process.wait()
			print "[STATUS] - Writing proteins"
			proteins = parseProteinsFromGFF3(outfile)
			for protein_name, protein_seq in proteins.items():
				for motif in MOTIFS:
					if motif in protein_seq:
						print ">" + species + "." + profile + "." + contig + "." + protein_name + "\n" + protein_seq
						ppx_results.write(">" + species + "." + profile + "." + contig + "." + protein_name + "\n" + protein_seq + "\n")
						break
			break
	ppx_results.close()

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
	dict_of_profiles = {}
	for profile in os.listdir(profile_dir + "/"):
		if profile.endswith(".prfl"):
			profile_name = profile.split(".")[0]
			dict_of_profiles[profile_name]= profile_dir + profile
	return dict_of_profiles

if __name__ == "__main__":
	try:
		contig_file = sys.argv[1]
		#profile = sys.argv[2]
		profile_dir = sys.argv[2]
		species = sys.argv[3] # ID for contigs, etc
		modus = sys.argv[4] 
		buffer_range = sys.argv[5] 
	except:
		sys.exit("Usage: ./ppx_wrapper.py [CONTIGFILE] [PROFILE_DIR] [SPECIES] [SEARCH|NOSEARCH] [BUFFERRANGE]")
	
	GENOME_DIR = 'genome/'
	TEMP_DIR = 'temp/'
	FASTBLOCKSEARCH_DIR = 'fastblocksearch/'
	AUGUSTUS_DIR = 'augustus/'
	MOTIFS = ["ELEKEF", "WFQNRR"]
	RESULTS_DIR = 'results/'

	dict_of_profiles = getProfiles(profile_dir)
	contigs = parse_contigs_to_dict(contig_file)
	#print contigs
	for profile in dict_of_profiles:
		if modus == "SEARCH":
			fastblocksearch(dict_of_profiles[profile], contigs)
		else: 
			pass
	runAugustusPPX()