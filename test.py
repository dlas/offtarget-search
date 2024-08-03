
from searchlib.search import Search, Match, EnumerateAmbiguous

from searchlib.dataman import GenomeManager, MockGenomeManager
import matplotlib.pyplot as plt
import numpy as np

def generate_report(results, spacer, pam):


	hit_counts = {}
	result_keys = []
	for r in results:
		if r.pam in hit_counts:
			hit_counts[r.pam]+=1
		else:
			hit_counts[r.pam]= 0
			result_keys.append(r.pam)
	
	fig, ax = plt.subplots()

	ax.barh(np.arange(len(result_keys)), [hit_counts[x] for x in result_keys])
	ax.set_yticks(np.arange(len(result_keys)), result_keys)
	ax.invert_yaxis()
	ax.set_title(f"spacer: {spacer} PAM: {pam}") 
	fig.savefig("x1.pdf")



def test_everything():


	M = MockGenomeManager()

	r = Search(M, "GGG", "AAAAAAAAAAAAAAAAAAAAA")
	generate_report(r, "GGG", "A*")
	Search(M, "SGG", "AAAAAAAAAAAAATAAAAAA")

	x = []
	EnumerateAmbiguous(list("AYTW"), 0, x)
	print(x)
	G = GenomeManager("/home/dstaff/code/gene/data/tiny_chr21.fa.gz", "/home/dstaff/code/gene/data/gnomad_af01_snps_chr21.vcf.bgz", "chr21")

	G.load()

	Search(G, "AT", "TTGTGACTGAAGGGC")
	Search(G, "AT", "TTGTGACTGTAGGGC")
	Search(G, "AT", "TTGTGACTTAAGGGC")


def more_real_test():
	
	G = GenomeManager("/home/dstaff/code/gene/data/chr21.fa.gz",  "/home/dstaff/code/gene/data/gnomad_af01_snps_chr21.vcf.bgz", "chr21")

	G.load()

	spacer = "ACACGTGTA"
	pam = "CTAN"
	Search(G, pam, spacer)



more_real_test()
#test_everything()
