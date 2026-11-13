import os
import random
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vcfart.annotation import Transcript  # noqa: E402
from vcfart.reference import FastaReference  # noqa: E402

# ---------------------------------------------------------------------------
# Prototype-regression target: a minus-strand gene with a 132-bp intron /
# 162-bp exon / 150-bp intron block, laid out exactly like the authors' original
# MAVE scripts expected it (plus-strand bases listed in DECREASING genomic order).
# The sequence is synthetic (seeded) with canonical AG / GT splice motifs.
#   upstream intron   : chrT:10294 -> 10163  (c.301-132 .. c.301-1)
#   target exon       : chrT:10162 -> 10001  (c.301 .. c.462)
#   downstream intron : chrT:10000 -> 9851   (c.462+1 .. c.462+150)
# ---------------------------------------------------------------------------
PROTO_TOP = 10294
PROTO_BOTTOM = 9851
PROTO_EXONS_TX_ORDER = [
    (15001, 15350),   # exon 1: 50 nt 5'UTR + 300 coding  (c.-50 .. c.300)
    (10001, 10162),   # exon 2: c.301 .. c.462  (target)
    (5001, 5200),     # exon 3: 138 coding + 62 nt 3'UTR
]
PROTO_CONTIG_LEN = 20000


def _proto_blocks():
    """Plus-strand bases in decreasing genomic order: (intron_up, exon, intron_down)."""
    rng = random.Random(1234)

    def rand(n):
        return "".join(rng.choice("ACGT") for _ in range(n))
    # transcript strand: ...AG | exon | GT...
    up_t = rand(130) + "AG"
    ex_t = rand(162)
    dn_t = "GT" + rand(148)
    from vcfart.sequtils import complement
    # plus-strand base at each position, listed from high to low coordinate
    return complement(up_t), complement(ex_t), complement(dn_t)


PROTO_INTRON_UP, PROTO_EXON, PROTO_INTRON_DOWN = _proto_blocks()


def proto_genepred_line(bin_col=True):
    exons = sorted(PROTO_EXONS_TX_ORDER)
    cds_end = 15350 - 50                # 1-based highest coding base (c.1)
    cds_start = 5200 - 138 + 1          # 1-based lowest coding base (end of stop codon)
    f = ["NM_PROTO.1", "chrT", "-", str(exons[0][0] - 1), str(exons[-1][1]),
         str(cds_start - 1), str(cds_end), str(len(exons)),
         ",".join(str(s - 1) for s, _ in exons) + ",",
         ",".join(str(e) for _, e in exons) + ",",
         "0", "PROTO1", "cmpl", "cmpl", ",".join("0" for _ in exons) + ","]
    return "\t".join((["585"] if bin_col else []) + f)


def write_fasta(path, chroms, width=60):
    with open(path, "w") as fh:
        for name, seq in chroms:
            fh.write(f">{name} test\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


@pytest.fixture(scope="session")
def proto_env(tmp_path_factory):
    d = tmp_path_factory.mktemp("proto")
    rng = random.Random(7)
    seq = list("".join(rng.choice("ACGT") for _ in range(PROTO_CONTIG_LEN)))
    block = (PROTO_INTRON_UP + PROTO_EXON + PROTO_INTRON_DOWN)[::-1]   # increasing order
    assert len(block) == PROTO_TOP - PROTO_BOTTOM + 1
    seq[PROTO_BOTTOM - 1:PROTO_TOP] = list(block)
    fa = d / "mini.fa"
    write_fasta(fa, [("chrT", "".join(seq))])
    ann = d / "proto_refGene.txt"
    ann.write_text(proto_genepred_line() + "\n")
    return {"fasta": str(fa), "annotation": str(ann), "dir": d}


# ---------------------------------------------------------------------------
# synthetic genome with a real ORF split over exons, on both strands
# ---------------------------------------------------------------------------
def _random_orf(rng, n_codons):
    from vcfart.sequtils import CODON_TABLE
    sense = [c for c, a in CODON_TABLE.items() if a != "*"]
    return "ATG" + "".join(rng.choice(sense) for _ in range(n_codons)) + "TAA"


@pytest.fixture(scope="session")
def synth_env(tmp_path_factory):
    rng = random.Random(7)
    d = tmp_path_factory.mktemp("synth")
    from vcfart.sequtils import revcomp

    def rand(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    utr5, utr3 = rand(40), rand(60)
    cds = _random_orf(rng, 98)          # 300 nt CDS
    tx_seq = utr5 + cds + utr3          # 400 nt
    exon_lens = [70, 120, 90, 120]      # 4 exons
    exon_seqs, i = [], 0
    for L in exon_lens:
        exon_seqs.append(tx_seq[i:i + L])
        i += L

    def intron(n):
        # canonical GT...AG with a poly-A stretch for indel-normalization tests
        return "GT" + rand(10) + "AAAAAA" + rand(n - 20) + "AG"

    # plus-strand gene on chr1 starting at 1001
    pieces, exons, pos = [rand(1000)], [], 1001
    for k, es in enumerate(exon_seqs):
        pieces.append(es)
        exons.append((pos, pos + len(es) - 1))
        pos += len(es)
        if k < len(exon_seqs) - 1:
            it = intron(200)
            pieces.append(it)
            pos += len(it)
    pieces.append(rand(1000))
    chr1 = "".join(pieces)
    # compute CDS genomic bounds via a temporary transcript
    t_plus = Transcript("NM_PLUS.1", "chr1", "+", exons, None, None, "PLUSG")
    cds_lo = t_plus.n_to_g(41)
    cds_hi = t_plus.n_to_g(340)

    # minus-strand gene on chr2: build plus-strand layout of the reverse complement
    rc_pieces = []
    rc_exons_tx = []
    for k, es in enumerate(exon_seqs):
        rc_pieces.append(es)
        if k < len(exon_seqs) - 1:
            rc_pieces.append(intron(220))
    tx_region = "".join(rc_pieces)               # transcript-strand gene body
    lead = rand(1000)
    chr2 = lead + revcomp(tx_region) + rand(1000)
    # transcript-order exon coordinates on chr2
    top = len(lead) + len(tx_region)             # genomic pos of first tx base
    p = top
    for k, piece in enumerate(rc_pieces):
        if k % 2 == 0:
            rc_exons_tx.append((p - len(piece) + 1, p))
        p -= len(piece)
    t_minus = Transcript("NM_MINUS.2", "chr2", "-", rc_exons_tx, None, None, "MINUSG")
    m_hi = t_minus.n_to_g(41)
    m_lo = t_minus.n_to_g(340)

    fa = d / "synth.fa"
    write_fasta(fa, [("chr1", chr1), ("chr2", chr2)], width=50)
    ann = d / "synth_refGene.txt"

    def gp(name, chrom, strand, ex, lo, hi, gene):
        ex = sorted(ex)
        return "\t".join([name, chrom, strand, str(ex[0][0] - 1), str(ex[-1][1]), str(lo - 1),
                          str(hi), str(len(ex)), ",".join(str(s - 1) for s, _ in ex) + ",",
                          ",".join(str(e) for _, e in ex) + ",", "0", gene, "cmpl", "cmpl",
                          ",".join("0" for _ in ex) + ","])

    ann.write_text(gp("NM_PLUS.1", "chr1", "+", exons, cds_lo, cds_hi, "PLUSG") + "\n"
                   + gp("NM_MINUS.2", "chr2", "-", rc_exons_tx, m_lo, m_hi, "MINUSG") + "\n")
    return {"fasta": str(fa), "annotation": str(ann), "tx_seq": tx_seq, "cds": cds,
            "dir": d}


@pytest.fixture(scope="session")
def synth_ref(synth_env):
    return FastaReference(synth_env["fasta"])
