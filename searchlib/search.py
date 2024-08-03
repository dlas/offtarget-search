
from .dataman import GenomeManager, MockGenomeManager
from copy import deepcopy
from dataclasses import dataclass

@dataclass
class OffTargetMatch:
	locus: str
	pam: str

def IsConcrete(s: str):
	for c in s:
		if c not in "AGCT":
			return False
	return True


def test_IsConcrete():
	assert(IsConcrete("AAAS") == False)
	assert(IsConcrete("ACCCTG") == True)



def EnumerateAmbiguous(pam: list[str], split, output):
	ALLOWED = "ACTGRYSWKMBDHVN"
	MAP = {"R": ["A", "G"],
		"Y": ["C", "T"],
		"S": ["G", "C"],
		"W": ["A", "T"],
		"K": ["G", "T"],
		"M": ["A", "C"],
		"B": ["C", "G", "T"],
		"D": ["A", "G", "C"],
		"H": ["A", "C", "T"],
		"V": ["A", "C", "G"], 
		"N": ["A", "C", "T", "G"]
	}

	if (IsConcrete(pam)):
		output.append("".join(pam))

	for i in range(split, len(pam)):
		c = pam[i]
		if (c not in ALLOWED):
			raise Exception("Unexpected base character in PAM string")
		if (c in MAP):
			newpam = deepcopy(pam)
			for x in MAP[c]:
				newpam[i] = x
				EnumerateAmbiguous(newpam, i, output)


def EnumerateAmbiguousWrapper(pam: str):
	output = []
	EnumerateAmbiguous(list(pam), 0, output)
	return sorted(output)

def test_EnumerateAmbiguous():
	assert(EnumerateAmbiguousWrapper("NRG") == ["AAG", "AGG", "CAG", "CGG", "GAG", "GGG", "TAG", "TGG"])

def Search(G: GenomeManager, pam: str, spacer:str):
	expanded_pams = []
	EnumerateAmbiguous(list(pam), 0, expanded_pams)

	matches = []

	for p in expanded_pams:
		matches += SearchExactPAM(G, p, spacer)

	return matches

def SearchExactPAM(G: GenomeManager, pam: str, spacer: str):
	#for idx in range(0, G.length() - (len(pam) + len(spacer))):

	matches = []

	#idx = G.find_after(0, pam)
	idx = len(spacer)
	while (idx >= 0 and (idx + len(pam) + len(spacer) < G.length())):
		matched = Match(G, idx-len(spacer), pam, spacer)
		if matched:
			print(f"{idx} {pam}")
			matches.append(OffTargetMatch(idx, pam))
		#idx = G.find_after(idx+1, pam)
		idx += 1

	return matches

def does_alternative_match(base, alts):
	for a in alts:
		if (base == a.alternative):
			return True
	return False

def CountMismatches(G, start, seq, test_seq, maximum):
	mismatches = 0
	for idx in range(0, len(test_seq)):
		if seq[idx] != test_seq[idx]:
			if not does_alternative_match(test_seq[idx], G.vcfalternatives(start+idx)):
				mismatches +=1
		if mismatches > maximum:
			return False
	return True

def Match(G: GenomeManager, start, pam, spacer):
	seq = G.subsequence(start, len(pam) + len(spacer))

	if not CountMismatches(G, start, seq[-(len(pam)):], pam, 0):
		return False

	#if seq[-((len(pam))):] != pam:
	#	return False

	mismatches = 0

	for idx in range(0, len(spacer)):
		if seq[idx] != spacer[idx]:
			if not does_alternative_match(spacer[idx], G.vcfalternatives(start+idx)):
				mismatches+=1

		if (mismatches > 1):
			return False

	return True

def test_Match():
	M = MockGenomeManager()

	assert(Match(M, 0, "GGG",        "AAAAAAAAAAAAAAAAAAAAA"))
	assert(Match(M, 0, "GGG", 	 "AAATAAAAAAAAAAAAAAAAA"))
	assert(not Match(M, 0, "TTT", 	 "AAAAAAAAAAAAAAAAAAAAA"))
	assert(not Match(M, 24, "GGG", 	 "AAAAAAAAAAAAAAAAAAAAA"))
	assert(Match(M, 96, "TTT", 	 "AAAAAAAAAAAAAAAAAAAAA"))
	assert(Match(M, 144, "GGG", 	 "AAAAAAAAAAAAAAAAAAAAA"))
	assert(Match(M, 144, "GGG", 	 "AGGAAAAAAAAAAAAAAAAAA"))
	assert(not Match(M, 144, "GGG",  "AAAAAAAAAAAAAAAGAAAAA"))
	assert(Match(M, 192, "GGG",  "AAAAAAAAAAAAAAAAAAAAA"))

