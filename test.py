
from py.search import Search, Match, EnumerateAmbiguous

from py.dataman import GenomeManager, MockGenomeManager


def test_everything():


	M = MockGenomeManager()

	Search(M, "GGG", "AAAAAAAAAAAAAAAAAAAA")
	Search(M, "SGG", "AAAAAAAAAAAAATAAAAAA")

	x = []
	EnumerateAmbiguous(list("AYTW"), 0, x)
	print(x)
	G = GenomeManager("/home/dstaff/code/gene/data/tiny_chr21.fa.gz", "/home/dstaff/code/gene/data/tiny_vcf.vcf.bgz", "chr21")

	G.load()

	Search(G, "AT", "TTGTGACTGAAGGGC")
	Search(G, "AT", "TTGTGACTGTAGGGC")
	Search(G, "AT", "TTGTGACTTAAGGGC")

test_everything()
