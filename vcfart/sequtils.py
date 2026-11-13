"""Small sequence helpers: reverse complement, genetic code, codon choices."""

from itertools import product

BASES = "ACGT"

_COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def complement(seq):
    return seq.translate(_COMP)


# Standard genetic code (NCBI table 1)
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_ORDER = "TCAG"
CODON_TABLE = {}
_i = 0
for _a in _ORDER:
    for _b in _ORDER:
        for _c in _ORDER:
            CODON_TABLE[_a + _b + _c] = _AA[_i]
            _i += 1

AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "*": "Ter", "X": "Xaa",
}

# Most frequently used human codon per residue (Kazusa codon usage, Homo sapiens).
PREFERRED_CODON = {
    "A": "GCC", "R": "AGA", "N": "AAC", "D": "GAC", "C": "TGC", "Q": "CAG",
    "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC", "L": "CTG", "K": "AAG",
    "M": "ATG", "F": "TTC", "P": "CCC", "S": "AGC", "T": "ACC", "W": "TGG",
    "Y": "TAC", "V": "GTG", "*": "TGA",
}

# Human codon usage per thousand (Kazusa), used for tie-breaking.
CODON_USAGE = {
    "TTT": 17.6, "TTC": 20.3, "TTA": 7.7, "TTG": 12.9, "CTT": 13.2, "CTC": 19.6,
    "CTA": 7.2, "CTG": 39.6, "ATT": 16.0, "ATC": 20.8, "ATA": 7.5, "ATG": 22.0,
    "GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.1, "TCT": 15.2, "TCC": 17.7,
    "TCA": 12.2, "TCG": 4.4, "CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9,
    "ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1, "GCT": 18.4, "GCC": 27.7,
    "GCA": 15.8, "GCG": 7.4, "TAT": 12.2, "TAC": 15.3, "TAA": 1.0, "TAG": 0.8,
    "CAT": 10.9, "CAC": 15.1, "CAA": 12.3, "CAG": 34.2, "AAT": 17.0, "AAC": 19.1,
    "AAA": 24.4, "AAG": 31.9, "GAT": 21.8, "GAC": 25.1, "GAA": 29.0, "GAG": 39.6,
    "TGT": 10.6, "TGC": 12.6, "TGA": 1.6, "TGG": 13.2, "CGT": 4.5, "CGC": 10.4,
    "CGA": 6.2, "CGG": 11.4, "AGT": 12.1, "AGC": 19.5, "AGA": 12.2, "AGG": 12.0,
    "GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
}

ALL_CODONS = ["".join(c) for c in product(BASES, repeat=3)]
NNK_CODONS = [a + b + c for a in BASES for b in BASES for c in "GT"]

SYNONYMS = {}
for _cod, _aa in CODON_TABLE.items():
    SYNONYMS.setdefault(_aa, []).append(_cod)


def translate_codon(codon):
    return CODON_TABLE.get(codon.upper(), "X")


def translate(seq, stop_at_stop=False):
    out = []
    for i in range(0, len(seq) - 2, 3):
        aa = translate_codon(seq[i:i + 3])
        out.append(aa)
        if stop_at_stop and aa == "*":
            break
    return "".join(out)


def aa3(seq):
    return "".join(AA3.get(a, "Xaa") for a in seq)


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def pick_codon(aa, wt_codon, rule="usage"):
    """Choose one codon encoding ``aa`` (never the WT codon itself).

    rule = "usage"    -> most frequently used human codon (next best if it equals WT)
    rule = "min-edit" -> fewest nucleotide changes from WT, ties broken by usage
    Returns None if no alternative codon exists (e.g. Met/Trp synonymous).
    """
    cands = [c for c in SYNONYMS[aa] if c != wt_codon]
    if not cands:
        return None
    if rule == "min-edit":
        cands.sort(key=lambda c: (hamming(c, wt_codon), -CODON_USAGE[c]))
    else:
        pref = PREFERRED_CODON[aa]
        cands.sort(key=lambda c: (c != pref, -CODON_USAGE[c]))
    return cands[0]
