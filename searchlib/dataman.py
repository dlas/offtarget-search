
import gzip
import logging
from dataclasses import dataclass

import vcfpy
from Bio import SeqIO
@dataclass

class Variant:
	alternative:str
	rsid: str
	freq:float

class GenomeManager:
	""" This is a simple class that holds the concept of an augmented reference genome. 
	It can be loaded with a genome frmo a fasta file and a set of variants from a vcf file
	These are loaded into memory and we offer convenient functions to retrieve a subsequence
	or to retrive the variants at a specific locale"""	
	def __init__(self, fasta_path, vcf_path, contig):
		self.fasta_path = fasta_path
		self.vcf_path = vcf_path
		self.variants_for_position = {}
		self.contig = contig
		self.fasta_data = None
			
	def load(self):
		""" Actually load the data"""
		self.loadfa()
		self.loadvcf()

	def loadfa(self):
		""" Load the fasta file."""

		# Use the fasta parser to load the whole file into memory and stash the contig that
		# we care about.
		with gzip.open(self.fasta_path, "rt") as handle:
			all_fasta_data = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
			self.fasta_data = all_fasta_data[self.contig]
		
	def loadvcf(self):
		""" Load the VCF file"""

		#  Open the file
		reader = vcfpy.Reader.from_path(self.vcf_path)
		
		# Copy the contents of the VCF file into a map of Variant objects
		for record in reader:
			if record.is_snv():
				pos = record.POS
				for a in record.ALT:
					v_record = Variant(a.value, record.ID, 0)
					if (pos in self.variants_for_position):
						self.variants_for_position[pos].append(v_record)
					else:
						self.variants_for_position[pos]=[v_record]
	
	def subsequence(self, start, length):
		return self.fasta_data[start:start+length].seq

	def find_after(self, start, s):
		return self.fasta_data.seq.find(s, start)
	
	def vcfalternatives(self, position):
		if position in self.variants_for_position:
			return self.variants_for_position[position]
		else:
			return []

	def length(self):
		return len(self.fasta_data)

class MockGenomeManager:
	""" This is a mock GenomeManager class with a small amount of data that we can use
	for unit tests.
	"""
	def __init__(self):
		self.variants_for_position = {
			74: [Variant("A", "", 0)],
			144: [Variant("A", "", 0)],
			145: [Variant("A", "", 0)],
			146: [Variant("A", "", 0)],
			215: [Variant("G", "", 0)]

		}

		self.seq = (""
			"AAAAAAAAAAAAAAAAAAAAAGGG" #0
			"AAAAAAAAAAAAAAAAAAAAACGG" #24
			"AGAAAAAAAAAAAAAAAAAAACGG" #48
			"AGGAAAAAAAAAAAAAAAAAAGGG" #72
			"AAAAAAAAAAAAAAAAAAAAATTT" #96
			"TGGAAAAAAAAAAAAAAAAAAGGG" #120
			"GGGGAAAAAAAAAAAAAAAAAGGG" #144
			"CCCCCCCCCCCCCCCCCCCAAAAA" #168
			"AAAAAAAAAAAAAAAAAGAAAGGA") #192
		print(f"-{self.seq}-")

	def subsequence(self, start, length):
		return self.seq[start:start+length]

	def find_after(self, start, s):
		return self.seq.find(s, start)

	def vcfalternatives(self, position):
		if position in self.variants_for_position:
			return self.variants_for_position[position]
		else:
			return []

	def length(self):
		return len(self.seq)
