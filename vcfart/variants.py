"""Variant representation, normalization and the mutation generators."""

from dataclasses import dataclass, field
from itertools import product

from .sequtils import ALL_CODONS, BASES, NNK_CODONS, PREFERRED_CODON, pick_codon, revcomp


@dataclass(frozen=True)
class Edit:
    """Replace genomic [start, end] (1-based, inclusive) with ``alt`` (plus strand).

    Insertions are empty intervals: end == start - 1 means "insert between
    start-1 and start".
    """
    chrom: str
    start: int
    end: int
    alt: str

    @property
    def ref_len(self):
        return self.end - self.start + 1


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    edit: Edit                      # left-normalized edit (matches the VCF record)
    raw: Edit                       # edit as first generated
    origins: list = field(default_factory=list)
    raws: dict = field(default_factory=dict)   # region label -> edit as generated there
    sets: list = field(default_factory=list)
    mut_type: str = ""
    hgvs_c: str = ""
    hgvs_p: str = ""
    consequence: str = ""
    region: str = ""

    @property
    def key(self):
        return (self.chrom, self.pos, self.ref, self.alt)


def ref_of(edit, ref):
    return ref.fetch(edit.chrom, edit.start, edit.end)


def trim(edit, ref):
    r = ref_of(edit, ref)
    a = edit.alt
    s = edit.start
    while r and a and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    while r and a and r[0] == a[0]:
        r, a, s = r[1:], a[1:], s + 1
    return Edit(edit.chrom, s, s + len(r) - 1, a)


def shift(edit, ref, direction):
    """Shift a pure insertion/deletion left (-1) or right (+1) through repeats."""
    e = trim(edit, ref)
    r_len = e.ref_len
    chrom, s, t, a = e.chrom, e.start, e.end, e.alt
    if r_len > 0 and a:
        return e
    clen = ref.length(chrom)
    if r_len > 0:  # deletion
        if direction < 0:
            while s > 1 and ref.base(chrom, s - 1) == ref.base(chrom, t):
                s, t = s - 1, t - 1
        else:
            while t < clen and ref.base(chrom, t + 1) == ref.base(chrom, s):
                s, t = s + 1, t + 1
        return Edit(chrom, s, t, "")
    if not a:
        return e
    if direction < 0:  # insertion between s-1 and s
        while s > 1 and ref.base(chrom, s - 1) == a[-1]:
            a = a[-1] + a[:-1]
            s -= 1
    else:
        while s <= clen and ref.base(chrom, s) == a[0]:
            a = a[1:] + a[0]
            s += 1
    return Edit(chrom, s, s - 1, a)


def to_vcf(edit, ref):
    """Left-normalized edit -> (POS, REF, ALT) with an anchor base for indels."""
    r = ref_of(edit, ref)
    if r and edit.alt:
        return edit.start, r, edit.alt
    p = edit.start - 1
    b = ref.base(edit.chrom, p)
    return p, b + r, b + edit.alt


def classify(edit, ref):
    r_len, a_len = edit.ref_len, len(edit.alt)
    if r_len == 1 and a_len == 1:
        return "SNV"
    if r_len == a_len:
        return "MNV"
    if a_len == 0:
        return "DEL"
    if r_len == 0:
        return "INS"
    return "DELINS"


def apply_edit(seq, seq_start, edit):
    """Apply an edit to ``seq`` whose first base is genomic ``seq_start``."""
    i = edit.start - seq_start
    j = edit.end - seq_start + 1
    if i < 0 or j > len(seq) or i > len(seq):
        return None
    return seq[:i] + edit.alt + seq[j:]


# ================================================================= generators
def gen_snv(region, ref):
    for g in range(region.start, region.end + 1):
        b = ref.base(region.chrom, g)
        if b not in BASES:
            continue
        for alt in BASES:
            if alt != b:
                yield Edit(region.chrom, g, g, alt), "snv"


def gen_del(region, ref, sizes):
    for k in sizes:
        for g in range(region.start, region.end - k + 2):
            if "N" in ref.fetch(region.chrom, g, g + k - 1):
                continue
            yield Edit(region.chrom, g, g + k - 1, ""), f"del{k}"


def gen_ins(region, ref, sizes):
    for k in sizes:
        seqs = ["".join(p) for p in product(BASES, repeat=k)]
        for g in range(region.start, region.end):
            for s in seqs:
                yield Edit(region.chrom, g + 1, g, s), f"ins{k}"


def gen_codon(region, ref, tx, codon_set="aa", pick="usage", stats=None):
    """Codon-level saturation (deep mutational scanning style).

    codon_set: 'aa'  -> one codon per amino acid (20 aa + Ter) incl. a synonymous WT codon
               'all' -> all 63 alternative codons
               'nnk' -> NNK codons
    """
    if not tx.coding:
        return
    cds = tx.cds_seq(ref)
    n_codons = tx.cds_len // 3
    for ci in range(n_codons):
        gs = [tx.n_to_g(tx.cds_start_n + ci * 3 + j) for j in range(3)]
        if not all(region.start <= g <= region.end for g in gs):
            continue
        step = 1 if tx.plus else -1
        if gs[1] - gs[0] != step or gs[2] - gs[1] != step:
            if stats is not None:
                stats["split_codons"] = stats.get("split_codons", 0) + 1
            continue
        wt = cds[ci * 3:ci * 3 + 3]
        if "N" in wt:
            continue
        if codon_set == "all":
            targets = [c for c in ALL_CODONS if c != wt]
        elif codon_set == "nnk":
            targets = [c for c in NNK_CODONS if c != wt]
        else:
            targets = []
            for aa in PREFERRED_CODON:
                c = pick_codon(aa, wt, pick)
                if c is not None:
                    targets.append(c)
        lo, hi = min(gs), max(gs)
        for c in targets:
            alt = c if tx.plus else revcomp(c)
            yield Edit(region.chrom, lo, hi, alt), "codon"
