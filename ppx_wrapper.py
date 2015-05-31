#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File    	: ppx_wrapper.py
Author  	: 
Version 	: 0.2
Description : 
To do 		:  ...
Bugs 		: Contig names in assemblies can't have "."
"""

from __future__ import division
import sys
import os
import time
import multiprocessing as mp 
import subprocess
import string
import collections
from subprocess import Popen

#######################################

class Files():
	def __init__(self, contig):
		self.contig = contig
		self.contig_file = TEMP_DIR + contig + ".fa"
		self.profile = None
		self.profile_file = None
		self.fastblockfile = None
		self.augustus = None
		self.gff3 = None
		self.protein = None

	def update(self, profile):
		self.profile = profile
		self.profile_file = PROFILE_DIR + profile + ".prfl"
		self.fastblockfile = FASTBLOCKSEARCH_DIR + self.contig + "." + profile + ".result"
		self.augustus = AUGUSTUS_DIR + self.contig + "." + profile + ".gff3"
		self.gff3 = RESULTS_DIR + self.contig + "." + profile + ".gff3"
		self.protein = RESULTS_DIR + self.contig + "." + profile + ".faa"

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

#######################################

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

#######################################

def getProfiles(profile_dir):
	dict_of_profiles = {}
	for profile in os.listdir(profile_dir + "/"):
		if profile.endswith(".prfl"):
			profile_name = profile.rstrip(".prfl")
			dict_of_profiles[profile_name]= profile_dir + profile
	return dict_of_profiles

def parse_contigs_to_dict(contig_file):
	contigs = set()
	header, seq = '', ''
	print "[STATUS] - Parsing contigs "
	with open(contig_file) as fh:
		for line in fh:
			if line.startswith(">"):
				if (seq):
					if len(seq) > 5000: 
						header = species + "." + header
						contig = ContigObject(header, seq.upper())
						contigs.add(contig) 
				seq = ''
				header = line.rstrip("\n").lstrip(">").replace(" ","_")	
			else:
 				seq += line.rstrip("\n")
		if len(seq) > 5000: 
			header = species + "." + header
			contig = ContigObject(header, seq)
			contigs.add(contig)
	
	file_dict = {}

	for contig in contigs:
		filename_obj = Files(contig.header)
		contig_file = filename_obj.contig_file
		temp = open(contig_file, 'w')
		temp.write(">" + contig.header + "\n" + contig.seq)
		temp.close()
		file_dict[contig.header]=filename_obj
	print "[STATUS] - %s contigs found " %len(contigs)
	return file_dict

#################################

def fastblocksearch(profile, contigs):
	
	profile_name = os.path.basename(profile).rstrip(".prfl")
	print "[STATUS] - Running FastBlockSearch with profile : " + profile_name
	processes = []
	
	jobs = []
	for contig in contigs:
		in_file = contigs[contig].contig_file
		out_file = contigs[contig].fastblockfile
		cmd = "/exports/software/augustus/augustus-3.0.3/bin/fastBlockSearch --cutoff=0.5 " + in_file + " " + profile + " > " + out_file + " "
		#cmd = "fastBlockSearch --cutoff=0.5 " + in_file + " " + profile + " > " + out_file + " "
		jobs.append(cmd)

	results = run_jobs(jobs, 24, pause = 2, verbose = False)

def run_jobs(jobs, threads, pause=2, verbose=False):
	"""Takes list of cmd strings, returns dict with error levels."""
	counter, max_value = 0, len(jobs)
	pending = jobs[:]
	running = []
	results = {}

	while pending or running:
	#See if any have finished
		for (cmd, process) in running:
			return_code = process.poll() #non-blocking
			if return_code is not None:
				results[cmd] = return_code
		running = [(cmd, process) for (cmd, process) in running if cmd not in results]
	
		if verbose:
			print "%i jobs pending, %i running, %i completed" % (len(pending), len(running), len(results))
		#See if we can start any new threads
	
		while pending and len(running) < threads:
			progress(counter, max_value)
			counter += 1
			cmd = pending.pop(0)
			if verbose:
				print cmd
			process = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
			running.append((cmd, process))
		#Loop...
		time.sleep(pause)

	if verbose:
		print "%i jobs completed" % len(results)
	assert set(jobs) == set(results)
	sys.stdout.write('\r')
	print "\tProgress : 100.00%"
	return results

''' This has to be done for each contig in each species (by species)
	- create folder for each genome
	- infer best overall hits	
'''

def progress(counter, max_value):
	sys.stdout.write('\r')
	progress = int(counter)/int(max_value)
	print "\tProgress : " + format(float(progress),'.2%'),
	sys.stdout.flush()

def parseFastBlockSearchResult(result_file, profile_name, contig_name):
	list_of_blocks = []
	contig = contig_name
	profile = profile_name
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
	#else:
	#	os.remove(result_file) 

	return list_of_blocks

def analyseBlocks(dict_of_blocks):
	# dict_of_blocks contains all blocks found on each sequence ... can have multiple block per profile per sequence

	fastblockresults_dict = {} # where the "good" blocks (hits) are stored, archived by contig (dict of lists)
	profile_hits = {} # where the "good" blocks (hits) are stored, archived by profile (dict of lists)
	
	max_profile_count = 1 # number of how many "good" blocks are stored per profile

	profile_count = {}	# number of blocks stored per profile 

	for score in sorted(dict_of_blocks, reverse=True):
		# for every score (in decreasing numerical order)
		
		print str(score) + "\t" + str(dict_of_blocks[score].__dict__) # for debugging
		
		block = dict_of_blocks[score] # get the block

		contig = block.contig # get the name of the contig
		profile = block.profile # get the name of the profile

		# THIS IS NEW ...
		overlap_threshold = 1000 # buffer range for getting coordinates (was previously the same as buffer range for predicting)
		# IDEA : change this to 1000 to allow for less distance between blocks competing for the same region, tested it and looked okay on C. elegans ...
		if not contig in fastblockresults_dict:
			# if we haven't seen this contig before  
			profile_count[profile] = profile_count.get(profile, 0) + 1 # increase count of blocks for this profile
			if profile_count[profile] <= max_profile_count: 
				# as long as we have not exceeded max_profile_count
				
				fastblockresults_dict[contig] = [] # populate fastblockresults_dict : key = contig, value = an empty list  
				fastblockresults_dict[contig].append(block) # add current block to the list in fastblockresults_dict  
				if not profile in profile_hits:
					# if we haven't seen a hit for this profile before 
					profile_hits[profile]=[] # populate profile_hits : key = profile, value = an empty list
				profile_hits[profile].append(block) # add current block to the list in profile_hits   
		else: 
			# if we have seen this contig before 
			for hit in fastblockresults_dict[contig]:
				# for each hit ("sane" block) that has already been put into fastblockresults_dict (they all have better score than the current one)
				hit_start, hit_end = int(hit.get('start', overlap_threshold)), int(hit.get('end', overlap_threshold)) # get coordinates of hit
				block_start, block_end = int(block.get('start', overlap_threshold)), int(block.get('end', overlap_threshold)) # get coordinates of current block
				coordinates = [hit_start, hit_end, block_start, block_end] # make a list with the coordinates 
				sum_lengths = (hit_end - hit_start) + (block_end - block_start) # sum up the lengths of bot regions (existing hit on contig and new block to be added)
				if sum_lengths >= (max(coordinates) - min(coordinates)):
					# "overlap between the two" if the sum of the lengths is greater or equal to the difference between maximal and minimal coordinate 
					pass # do nothing (there is already one hit with a higher score in that region)
				else:
					# There is either complete overlap or none at all
					if (hit_start >= block_start and hit_end <= block_end):
						# Hit is contained within block
						pass # do nothing (although the region of the new block is bigger than the hit) 
						# IDEA: one could consider making the hit longer using the coordinates of the block
					elif (block_start >= hit_start and block_end <= hit_end):
						# Block is contained within hit
						pass # do nothing (the hit is longer than the block)
					else:
						# There is no overlap
						profile_count[profile] = profile_count.get(profile, 0) + 1  # increase count of blocks for this profile
						if profile_count[profile] <= max_profile_count:
							# as long as we have not exceeded max_profile_count
							fastblockresults_dict[contig].append(block) # add current block to the list in fastblockresults_dict  
							if not profile in profile_hits:
								# if we haven't seen a hit for this profile before 
								profile_hits[profile]=[] # populate profile_hits : key = profile, value = an empty list
							profile_hits[profile].append(block) # add current block to the list in profile_hits   

	# Debugging
	ppx_log = open(RESULTS_DIR + species + ".log", "w") 
	for profile in profile_hits:
		ppx_log.write(profile)
		for hits in profile_hits[profile]:
			ppx_log.write(str(hits.__dict__) + "\n")
	
	return profile_hits # return dict of lists with the hits archived by profile

def runAugustusPPX(files):
	dict_of_blocks = {}
	print "Buffer range : " + str(buffer_range) 
	for result in os.listdir(FASTBLOCKSEARCH_DIR):
		# For each FastBlockSearch result file ...
		if result.startswith(species) and result.endswith(".result"):
			result_file = FASTBLOCKSEARCH_DIR + result
			print "[STATUS] - Parsing : " + result_file
			#print result.split(".")
			profile_name = result.split(".")[-2]
			contig_name = ".".join(result.split(".")[0:-2])
			#print profile_name, contig_name
			list_of_blocks = parseFastBlockSearchResult(result_file, profile_name, contig_name)
			if (list_of_blocks):
				for block in list_of_blocks:
					#dict_of_blocks[block.score] = block
					#dict_of_blocks[block.multi_score] = block
					dict_of_blocks[block.contig][block.profile][block.score] = block

	profile_hits = analyseBlocks(dict_of_blocks)

	ppx_results = open(RESULTS_DIR + species + ".faa", "w") # file to which to write the resulting proteins for all profiles

	for profile in profile_hits:
		print '*'
		# for each profile 
		for hit in profile_hits[profile]:
			print hit.__dict__
			# for each hit (block object)
			contig = hit.contig # get contig name
			start = str(hit.get('start', buffer_range)) # get start with buffer range
			end = str(hit.get('end', buffer_range)) # get stop with buffer range
			strand = hit.get('strand', 0) # get strand
			
			infile = TEMP_DIR + contig + ".fa" # file containing sequence of contig
			outfile = AUGUSTUS_DIR + contig + "." + profile + ".gff3" # file to which the output gff3 is written
			gff_of_gene_file = RESULTS_DIR + contig + "." + profile + ".gff3" # other gff3 file to which only the good (motif-containing) gene models are written
			profile_file = dict_of_profiles[profile] # get filename of profile
		
			print "[STATUS] - Calling protein \"" + profile_file + "\" in contig \"" + contig + "\" from " + str(start) + " to " + str(end) + " : " + outfile
			# Start running the processes for gene finding by prociding start, stop, strand, contig sequence, profile, and output file
			process = subprocess.Popen("/exports/software/augustus/augustus-3.0.3/bin/augustus --species=onchocerca_gutturosa --gff3=on --proteinprofile=" + profile_file + " --predictionStart=" + start + " --predictionEnd=" + end + " --strand=" + strand + " " + infile + " > " + outfile , stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
			#process = subprocess.Popen("augustus --species=onchocerca_gutturosa --gff3=on --proteinprofile=" + profile_file + " --predictionStart=" + start + " --predictionEnd=" + end + " --strand=" + strand + " " + infile + " > " + outfile , stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
			process.wait()
			

			proteins = parseProteinsFromGFF3(outfile) # parse proteins from the GFF3 output file

			# screening for motifs in proteins
			for protein_name, protein_seq in proteins.items():
				# for (protein name, protein seq) pair in proteins 
				for motif in MOTIFS:
					# for each MOTIF
					if motif in protein_seq:
						print "[STATUS] - Writing proteins"
						# if motif appeas in protein seq
						parseGFFFromGFF3(outfile, gff_of_gene_file) # parse relevant part about this gene model from the gff3 output file
						# printing sequence to screen
						print ">" + species + "." + profile + "." + contig + "." + str(start) + "-" + str(end) + "." + strand + "." + protein_name + "\n" + protein_seq
						ppx_results.write(">" + species + "." + profile + "." + contig + "." + str(start) + "-" + str(end) + "." + strand + "." + protein_name + "\n" + protein_seq + "\n") # file to which to write the resulting proteins for all profiles
						break # if you found a protein with one motif stop there, do not do twice because both motifs are present
			break # only do one (best) hit
	ppx_results.close()

	print "[STATUS] - Done."

def parseProteinsFromGFF3(augustus):

	read_mode = 0 # reading flag
	temp = '' # 
	protein_name = '' 
	protein_seq = ''
	dict_of_proteins = {} # dict : key = protein name, vale = protein seq

	with open(augustus) as fh:
		# while reading the gff3 file
		for line in fh: # for each line
			if line.startswith("# start gene "): 
				# if you see "# start gene " at the beginning of the line
				read_mode = 1 # read mode on
				protein_name = line.split(" ")[3].rstrip("\n") # get protein name
			elif line.startswith("# end gene "):
				# if you see "# end gene " at the beginning of the line
				dict_of_proteins[protein_name] = protein_seq # add key = protein name, value = protein seq to dict_of_proteins
				read_mode = 0 # read mode off
				protein_name = '' # empty protein name
				protein_seq = '' # empty protein seq
			elif read_mode == 1 and line.startswith("#"):
				# if read mode is on and the line starts with "#" (that where the AA's are)
				if "]" in line:
					read_mode = 0
				temp = line.lstrip("# protein sequence = [") # remove clutter from string
				temp = temp.lstrip("# ") # remove clutter from string
				temp = temp.rstrip("\n") # remove clutter from string
				temp = temp.rstrip("]") # remove clutter from string
				protein_seq += temp # add AA to protein seq
			else:
				pass # do nothing
	return dict_of_proteins # return dict : key = protein name, vale = protein seq

def parseGFFFromGFF3(outfile, gff_of_gene_file): # (augustus gff3, gff3 output file containing the ones that passed motif filtering )

	read_mode = 0 # reading flag
	fh_out = open(gff_of_gene_file, 'w') 
	fh_out.write("##gff-version 3\n") # write header
	with open(outfile) as fh: 
		for line in fh: 
			if line.startswith("# start gene "):
				read_mode = 1
				fh_out.write(line)
			elif line.startswith("# end gene "):
				read_mode = 0
				fh_out.write(line)
			elif read_mode == 1:
				fh_out.write(line)
			else:
				pass
	fh_out.close()

if __name__ == "__main__":
	try:
		contig_file = sys.argv[1]
		profile_dir = sys.argv[2]
		species = sys.argv[3] # ID for contigs, etc
		modus = sys.argv[4] 
		buffer_range = sys.argv[5] 
	except:
		sys.exit("Usage: ./ppx_wrapper.py [CONTIGFILE] [PROFILE_DIR] [SPECIES] [SEARCH|NOSEARCH] [BUFFERRANGE]")
	
	PROFILE_DIR = profile_dir 
	
	TEMP_DIR = '_temp/'
	FASTBLOCKSEARCH_DIR = '_fastblocksearch/'
	AUGUSTUS_DIR = '_augustus/'
	RESULTS_DIR = '_results/'
	DIRS = [TEMP_DIR, FASTBLOCKSEARCH_DIR, AUGUSTUS_DIR, RESULTS_DIR]

	for DIR in DIRS:
		if not os.path.exists(DIR):
			os.makedirs(DIR)

	MOTIFS = ["ELEK", "WFQNRR"]

	# 1. Get profiles
	dict_of_profiles = getProfiles(profile_dir)
	# 2. parse assemblies
	files = parse_contigs_to_dict(contig_file)
	
	# 3. Search for blocks	
	for profile in dict_of_profiles:
		for contig in files:
			files[contig].update(profile)
		if modus == "SEARCH":
			fastblocksearch(dict_of_profiles[profile], files)
		else: 
			pass

	runAugustusPPX(files)