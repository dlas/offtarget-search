
from Bio import SeqIO
import gzip
import vcfpy
from dataclasses import dataclass

@dataclass

class Variant:
	alternative:str
	rsid: str
	freq:float

class GenomeManager:
	
	def __init__(self, fasta_path, vcf_path, contig):
		self.fasta_path = fasta_path
		self.vcf_path = vcf_path
		self.variants_for_position = {}
		self.contig = contig
		self.fasta_data = None
			
	def load(self):
		self.loadfa()
		self.loadvcf()

	def loadfa(self):
		with gzip.open(self.fasta_path, "rt") as handle:
			all_fasta_data = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
			self.fasta_data = all_fasta_data[self.contig]
		
		print("fa loaded")
	def loadvcf(self):

		reader = vcfpy.Reader.from_path(self.vcf_path)
		

		for record in reader:
			if record.is_snv():
				pos = record.POS
				for a in record.ALT:
					v_record = Variant(a.value, record.ID, 0)
					if (pos in self.variants_for_position):
						self.variants_for_position[pos].append(v_record)
					else:
						self.variants_for_position[pos]=[v_record]
	
		print("vcf loaded")


	def subsequence(self, start, length):
		#return self.fasta_data[contig][start:start+length]
		return self.fasta_data[start:start+length].seq

	
	def vcfalternatives(self, position):
		if position in self.variants_for_position:
			return self.variants_for_position[position]
		else:
			return []

	def length(self):
		return len(self.fasta_data)

class MockGenomeManager:
	def __init__(self):
		self.variants_for_position = {
			74: [Variant("A", "", 0)],
			144: [Variant("A", "", 0)],
			145: [Variant("A", "", 0)],
			146: [Variant("A", "", 0)]

		}

		self.seq = (""
			"AAAAAAAAAAAAAAAAAAAAAGGG" #0
			"AAAAAAAAAAAAAAAAAAAAACGG" #24
			"AGAAAAAAAAAAAAAAAAAAACGG" #48
			"AGGAAAAAAAAAAAAAAAAAAGGG" #72
			"AAAAAAAAAAAAAAAAAAAAATTT" #96
			"TGGAAAAAAAAAAAAAAAAAAGGG" #120
			"GGGGAAAAAAAAAAAAAAAAAGGG") # 144

	def subsequence(self, start, length):
		return self.seq[start:start+length]

	def vcfalternatives(self, position):
		if position in self.variants_for_position:
			return self.variants_for_position[position]
		else:
			return []

	def length(self):
		return len(self.seq)
