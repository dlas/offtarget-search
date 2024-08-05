#!/usr/bin/env python3

import argparse
import logging

import matplotlib.pyplot as plt
import numpy as np

from searchlib.search import Search, Match, EnumerateAmbiguous
from searchlib.dataman import GenomeManager, MockGenomeManager

def generate_report(results, spacer, pam, output):
	""" Generate a nice PDF report of matching results """


	# Count the results on a per PAM sequence basis.
	hit_counts = {}
	result_keys = []
	for r in results:
		if r.pam in hit_counts:
			hit_counts[r.pam]+=1
		else:
			hit_counts[r.pam]= 0
			result_keys.append(r.pam)
	
	# Render this as a simple bar graph
	fig, ax = plt.subplots()
	ax.barh(np.arange(len(result_keys)), [hit_counts[x] for x in result_keys])
	ax.set_yticks(np.arange(len(result_keys)), result_keys)
	ax.invert_yaxis()
	ax.set_title(f"spacer: {spacer} PAM: {pam}") 

	# Save to a file
	fig.savefig(output)


def do_main():
	logging.basicConfig(level=logging.INFO)
	parser = argparse.ArgumentParser()
	parser.add_argument("--fasta_reference_path")
	parser.add_argument("--vcf_augmentation_path")
	parser.add_argument("--contig")
	parser.add_argument("--pam")
	parser.add_argument("--spacer")
	parser.add_argument("--pdf_output")

	args = parser.parse_args()
	
	G = GenomeManager(args.fasta_reference_path, args.vcf_augmentation_path, args.contig)

	logging.info("Loading...")

	G.load()

	logging.info(f"Searching {args.fasta_reference_path} (contig: {args.contig} augmented by {args.vcf_augmentation_path} for PAM sequence {args.pam} after spacer {args.spacer}")

	reports = Search(G, args.pam, args.spacer)

	logging.info(f"Search completed with {len(reports)} matches")

	if (len(reports) == 0):
		print("No matches found.")
	else:
		print("locus\tpam")
		for r in reports:
			print(f"{args.contig}:{r.locus}\t{r.pam}")
			

		logging.info("Formatting PDF file {args.pdf_output}")
		generate_report(reports, args.spacer, args.pam, args.pdf_output)


if __name__ == "__main__":
	do_main()
