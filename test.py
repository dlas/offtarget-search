
from searchlib.search import Search, Match, EnumerateAmbiguous

from searchlib.dataman import GenomeManager, MockGenomeManager



def test_everything():


	M = MockGenomeManager()

	Search(M, "GGG", "AAAAAAAAAAAAAAAAAAAAA")
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
	
	G = GenomeManager("/home/dstaff/code/gene/data/chr21.fa.gz", "/home/dstaff/code/gene/data/tiny_vcf.vcf.bgz", "chr21")

	G.load()

	Search(G, "CTAN", "ACACGTGTA")



more_real_test()
#test_everything()
