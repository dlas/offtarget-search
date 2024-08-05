from .dataman import GenomeManager, MockGenomeManager
from copy import deepcopy
from dataclasses import dataclass

import logging


@dataclass
class OffTargetMatch:
    locus: str
    pam: str


def IsConcrete(s: str) -> bool:
    """Is a string an unambiguous (ACTG) DNA sequence?"""
    for c in s:
        if c not in "AGCT":
            return False
    return True


def test_IsConcrete():
    assert IsConcrete("AAAS") == False
    assert IsConcrete("ACCCTG") == True


def EnumerateAmbiguous(pam: list[str], split, output):

    # What characters are allowed in a sequence?
    ALLOWED = "ACTGRYSWKMBDHVN"

    # What are the possible translations of ambiguous bases?
    MAP = {
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "C"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "C", "T", "G"],
    }

    # Is pam already concrete? then we're done. Add it to the output
    if IsConcrete(pam):
        output.append("".join(pam))

    # Iterate over every
    for i in range(split, len(pam)):
        c = pam[i]
        if c not in ALLOWED:
            raise Exception("Unexpected base character in PAM string")
        if c in MAP:
            newpam = deepcopy(pam)
            for x in MAP[c]:
                newpam[i] = x
                EnumerateAmbiguous(newpam, i, output)


def EnumerateAmbiguousWrapper(pam: str):
    """Take a string with an ambiguous sequence and generate
    an array of all possible realizations of that sequence.
    e.g. "AN" becomes ["AC", "AC", "AG", "AT"].
    """
    output = []
    EnumerateAmbiguous(list(pam), 0, output)
    return sorted(output)


def test_EnumerateAmbiguous():
    assert EnumerateAmbiguousWrapper("NRG") == ["AAG", "AGG", "CAG", "CGG", "GAG", "GGG", "TAG", "TGG"]
    assert EnumerateAmbiguousWrapper("ATTG") == ["ATTG"]


def Search(G: GenomeManager, pam: str, spacer: str):
    """Search for a given PAM + spacer combination in a given genome"""

    # Generate a list of all possible permutations of the PAM sequence to search for
    # This list is usually small (e.g. < 100) which would enable efficient seeding
    # for searches
    expanded_pams = []
    EnumerateAmbiguous(list(pam), 0, expanded_pams)

    matches = []

    # For each possible PAM sequence, search the entire genome for it
    # and the spacer and build up all the results in a nice array
    for p in expanded_pams:
        matches += SearchExactPAM(G, p, spacer)

    return matches


def SearchExactPAM(G: GenomeManager, pam: str, spacer: str):
    """Search for an exact PAM match and spacer match in an augmented genome"""

    matches = []

    # Simple linear search. idx is the place where the spacer ends
    # and the PAM sequence begins
    idx = len(spacer)
    while idx >= 0 and (idx + len(pam) + len(spacer) < G.length()):
        # Check if both match!
        matched = Match(G, idx - len(spacer), pam, spacer)
        if matched:
            logging.info(f"found: {idx} {pam}")
            matches.append(OffTargetMatch(idx, pam))
        # idx = G.find_after(idx+1, pam)
        idx += 1

    return matches


def does_alternative_match(base, alts):
    for a in alts:
        if base == a.alternative:
            return True
    return False


def CountMismatches(G, start, seq, test_seq, maximum):
    """Does 'test_seq' have more than maximum mismatches to
    seq?  (seq is augmented by variant information)
    """

    mismatches = 0
    for idx in range(0, len(test_seq)):
        ### Does test_seq match seq directly?
        if seq[idx] != test_seq[idx]:
            # It doesn't but check if any variants match
            if not does_alternative_match(test_seq[idx], G.vcfalternatives(start + idx)):
                mismatches += 1
        # Get out of here as soon as we hit the maximum
        if mismatches > maximum:
            return False
    return True


def Match(G: GenomeManager, start, pam, spacer):
    """do we consider [spacer][pam] as matching the genome at this locus?
    We consider a match if
            (1) PAM matches the augmented reference genome exactly
            (2) and, spacer matches the reference genome with at most one mismatch
    """

    seq = G.subsequence(start, len(pam) + len(spacer))

    if not CountMismatches(G, start + len(spacer), seq[-(len(pam)) :], pam, 0):
        return False

    return CountMismatches(G, start, seq, spacer, 1)


def test_Match():
    M = MockGenomeManager()

    assert Match(M, 0, "GGG", "AAAAAAAAAAAAAAAAAAAAA")
    assert Match(M, 0, "GGG", "AAATAAAAAAAAAAAAAAAAA")
    assert not Match(M, 0, "TTT", "AAAAAAAAAAAAAAAAAAAAA")
    assert not Match(M, 24, "GGG", "AAAAAAAAAAAAAAAAAAAAA")
    assert Match(M, 96, "TTT", "AAAAAAAAAAAAAAAAAAAAA")
    assert Match(M, 144, "GGG", "AAAAAAAAAAAAAAAAAAAAA")
    assert Match(M, 144, "GGG", "AGGAAAAAAAAAAAAAAAAAA")
    assert not Match(M, 144, "GGG", "AAAAAAAAAAAAAAAGAAAAA")
    assert Match(M, 192, "GGG", "AAAAAAAAAAAAAAAAAAAAA")
