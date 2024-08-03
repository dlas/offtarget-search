
from .dataman import GenomeManager
from copy import deepcopy

print("hi")

def IsConcrete(s: str):
	for c in s:
		if c not in "AGCT":
			return False
	return True



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


def Search(G: GenomeManager, pam: str, spacer:str):
	expanded_pams = []
	EnumerateAmbiguous(list(pam), 0, expanded_pams)

	for p in expanded_pams:
		SearchExactPAM(G, p, spacer)

def SearchExactPAM(G: GenomeManager, pam: str, spacer: str):
	for idx in range(0, G.length() - (len(pam) + len(spacer))):


		matched = Match(G, idx, pam, spacer)
		if matched:
			print(f"{idx} {pam}")

	

def does_alternative_match(base, alts):
	for a in alts:
		if (base == a.alternative):
			return True
	return False

def Match(G: GenomeManager, start, pam, spacer):
	seq = G.subsequence(start, len(pam) + len(spacer))
	#print(seq)
	#print(seq[0:2] == "NN")

	if seq[-(len(pam)):] != pam:
		return False

	mismatches = 0

	for idx in range(0, len(spacer)):
		if seq[idx] != spacer[idx]:
			if not does_alternative_match(spacer[idx], G.vcfalternatives(start+idx)):
				mismatches+=1

		if (mismatches > 1):
			return False

	return True
