import random

import pytest
from conftest import write_fasta

from vcfart.annotation import choose_transcript, load_genepred, load_gtf, load_transcripts
from vcfart.cli import main
from vcfart.hgvs import Annotator
from vcfart.reference import FastaReference
from vcfart.regions import parse_region
from vcfart.sequtils import AA3, complement, revcomp, translate
from vcfart.variants import Edit, Variant, apply_edit, classify, shift, to_vcf


def get_tx(env, name, ref=None):
    tx, _ = choose_transcript(load_genepred(env["annotation"], name), name, ref=ref)
    return tx


@pytest.mark.parametrize("name", ["NM_PLUS.1", "NM_MINUS.2"])
def test_roundtrip_and_sequence(synth_env, synth_ref, name):
    tx = get_tx(synth_env, name, synth_ref)
    assert tx.transcript_seq(synth_ref) == synth_env["tx_seq"]
    assert tx.cds_seq(synth_ref) == synth_env["cds"]
    for g in range(tx.tx_start - 60, tx.tx_end + 61):
        n, off = tx.g_to_n(g)
        assert tx.n_to_g(n, off) == g
        assert tx.parse_pos(tx.prefix + tx.pos(g)) == g
        if off == 0 and 1 <= n <= tx.length:
            b = synth_ref.base(tx.chrom, g)
            assert synth_env["tx_seq"][n - 1] == (b if tx.plus else complement(b))


@pytest.mark.parametrize("name", ["NM_PLUS.1", "NM_MINUS.2"])
def test_cds_numbering(synth_env, synth_ref, name):
    tx = get_tx(synth_env, name, synth_ref)
    assert tx.pos(tx.n_to_g(41)) == "1"            # 40 nt 5'UTR
    assert tx.pos(tx.n_to_g(40)) == "-1"
    assert tx.pos(tx.n_to_g(340)) == "300"
    assert tx.pos(tx.n_to_g(341)) == "*1"
    # intron 1 (between exon1 [n 1-70] and exon2): first/last bases
    assert tx.pos(tx.n_to_g(70, 1)) == "30+1"
    assert tx.pos(tx.n_to_g(71, -1)) == "31-1"
    # odd/even midpoint rule
    kind, k = tx.region(tx.n_to_g(70, 1))
    assert (kind, k) == ("intron", 1)


@pytest.mark.parametrize("name", ["NM_PLUS.1", "NM_MINUS.2"])
def test_snv_annotation_matches_independent_translation(synth_env, synth_ref, name):
    tx = get_tx(synth_env, name, synth_ref)
    ann = Annotator(tx, synth_ref)
    cds, tseq = synth_env["cds"], synth_env["tx_seq"]
    prot = translate(cds)
    for c in range(1, 301):
        g = tx.n_to_g(tx.cds_start_n + c - 1)
        b = synth_ref.base(tx.chrom, g)
        for alt in "ACGT":
            if alt == b:
                continue
            v = Variant(tx.chrom, g, b, alt, Edit(tx.chrom, g, g, alt), None)
            ann.annotate(v)
            tb, ta = (b, alt) if tx.plus else (complement(b), complement(alt))
            assert v.hgvs_c == f"c.{c}{tb}>{ta}"
            mut = cds[:c - 1] + ta + cds[c:]
            k = (c - 1) // 3
            w, m = prot[k], translate(mut[k * 3:k * 3 + 3])
            if k == 0:
                exp = "p.Met1?"
            elif w == m:
                exp = f"p.{AA3[w]}{k + 1}="
            elif w == "*":
                tail = translate((mut + tseq[340:])[k * 3:], stop_at_stop=True)
                j = tail.find("*")
                exp = f"p.Ter{k + 1}{AA3[tail[0]]}extTer{j + 1 if j >= 0 else '?'}"
            elif m == "*":
                exp = f"p.{AA3[w]}{k + 1}Ter"
            else:
                exp = f"p.{AA3[w]}{k + 1}{AA3[m]}"
            assert v.hgvs_p == exp, (c, alt, v.hgvs_p, exp)


@pytest.mark.parametrize("name", ["NM_PLUS.1", "NM_MINUS.2"])
def test_indel_normalization_equivalence(synth_env, synth_ref, name):
    """VCF-left-normalized, raw and HGVS-3'-shifted edits all give the same genome."""
    tx = get_tx(synth_env, name, synth_ref)
    ann = Annotator(tx, synth_ref)
    rng = random.Random(1)
    lo, hi = tx.tx_start - 20, tx.tx_end + 20
    window = synth_ref.fetch(tx.chrom, lo - 100, hi + 100)
    base0 = lo - 100
    for _ in range(600):
        g = rng.randint(lo, hi)
        if rng.random() < 0.5:
            k = rng.randint(1, 4)
            e = Edit(tx.chrom, g, g + k - 1, "")
        else:
            k = rng.randint(1, 3)
            # bias toward repeat-creating insertions
            seq = synth_ref.fetch(tx.chrom, g - k, g - 1) if rng.random() < 0.5 else \
                "".join(rng.choice("ACGT") for _ in range(k))
            e = Edit(tx.chrom, g, g - 1, seq)
        left = shift(e, synth_ref, -1)
        right = shift(e, synth_ref, +1)
        e3 = ann.hgvs_edit(left)
        seqs = {apply_edit(window, base0, x) for x in (e, left, right, e3)}
        assert len(seqs) == 1, (e, left, right, e3)
        # 3' rule direction
        assert e3 == (right if tx.plus else left)
        pos, r, a = to_vcf(left, synth_ref)
        assert r[0] == a[0] == synth_ref.base(tx.chrom, pos)


def _mini(tmp_path, tx_seq, exon_lens, cds_n, strand, intron="GTAAGTCCCCCCCCCCCTTTTTTTTTTCAG"):
    """Build a tiny genome+genePred for a transcript sequence (transcript strand)."""
    lead, tail = "ACGTTGCA" * 25, "TTGACCGA" * 25
    body, i, parts = [], 0, []
    for k, L in enumerate(exon_lens):
        body.append(tx_seq[i:i + L])
        parts.append(("e", L))
        i += L
        if k < len(exon_lens) - 1:
            body.append(intron)
            parts.append(("i", len(intron)))
    gene = "".join(body)
    if strand == "+":
        chrom = lead + gene + tail
        p, exons = len(lead) + 1, []
        for kind, L in parts:
            if kind == "e":
                exons.append((p, p + L - 1))
            p += L
    else:
        chrom = lead + revcomp(gene) + tail
        p, exons = len(lead) + len(gene), []
        for kind, L in parts:
            if kind == "e":
                exons.append((p - L + 1, p))
            p -= L
    fa = tmp_path / f"mini{strand}.fa"
    write_fasta(fa, [("chr7", chrom)])
    from vcfart.annotation import Transcript
    t = Transcript("NM_MINI.1", "chr7", strand, exons, None, None, "MINI")
    a, b = t.n_to_g(cds_n[0]), t.n_to_g(cds_n[1])
    lo, hi = min(a, b), max(a, b)
    ex = sorted(exons)
    gp = "\t".join(["NM_MINI.1", "chr7", strand, str(ex[0][0] - 1), str(ex[-1][1]), str(lo - 1),
                    str(hi), str(len(ex)), ",".join(str(s - 1) for s, _ in ex) + ",",
                    ",".join(str(e) for _, e in ex) + ",", "0", "MINI"])
    ann = tmp_path / f"mini{strand}.txt"
    ann.write_text(gp + "\n")
    ref = FastaReference(str(fa))
    tx, _ = choose_transcript(load_genepred(str(ann), "NM_MINI.1"), "NM_MINI.1", ref=ref)
    return tx, ref, str(fa), str(ann)


# 5'UTR(10) | ATG AAA GAA GAA GAA TGG CCC CCC TAA | 3'UTR with a later stop
MINI_UTR5 = "CCCAGCCACC"
MINI_CDS = "ATGAAAGAAGAAGAATGGCCCCCCTAA"
MINI_UTR3 = "CAGCTGGTCTAAGCAGCAGCAGCAGCAG"
MINI_TX = MINI_UTR5 + MINI_CDS + MINI_UTR3


def _edit_at_c(tx, ref, c_from, c_to, alt_t):
    """Build a genomic Edit from transcript-strand c. coordinates (c_to < c_from: insertion)."""
    g1 = tx.parse_pos(f"c.{c_from}")
    if c_to < c_from:   # insertion before c_from
        g0 = tx.parse_pos(f"c.{c_to}")
        s = max(g0, g1)
        return Edit(tx.chrom, s, s - 1, alt_t if tx.plus else revcomp(alt_t))
    g2 = tx.parse_pos(f"c.{c_to}")
    lo, hi = min(g1, g2), max(g1, g2)
    return Edit(tx.chrom, lo, hi, alt_t if tx.plus else revcomp(alt_t))


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("c_from,c_to,alt,hgvs_c,hgvs_p,csq", [
    (10, 12, "", "c.13_15del", "p.Glu5del", "inframe_deletion"),
    (16, 15, "GAA", "c.13_15dup", "p.Glu5dup", "inframe_insertion"),
    (5, 5, "", "c.6del", "p.Glu3LysfsTer?", "frameshift_variant"),
    (17, 17, "A", "c.17G>A", "p.Trp6Ter", "stop_gained"),
    (16, 16, "C", "c.16T>C", "p.Trp6Arg", "missense_variant"),
    (25, 25, "C", "c.25T>C", "p.Ter9GlnextTer5", "stop_lost"),
    (2, 2, "C", "c.2T>C", "p.Met1?", "start_lost"),
    (19, 21, "TGG", "c.19_21delinsTGG", "p.Pro7Trp", "missense_variant"),
    (19, 21, "TAA", "c.19_21delinsTAA", "p.Pro7Ter", "stop_gained"),
    (21, 20, "TAA", "c.20_21insTAA", "p.Pro7_Pro8insAsn", "inframe_insertion"),
    (22, 21, "TAA", "c.21_22insTAA", "p.Pro8Ter", "stop_gained"),
    (22, 21, "CCC", "c.22_24dup", "p.Pro8dup", "inframe_insertion"),
    (4, 3, "G", "c.3dup", "p.Lys2GlufsTer?", "frameshift_variant"),
    (25, 27, "", "c.25_27del", "p.Ter9GlnextTer4", "stop_lost"),
    (25, 24, "G", "c.24_25insG", "p.Ter9ValextTer?", "frameshift_variant"),
])
def test_protein_cases(tmp_path, strand, c_from, c_to, alt, hgvs_c, hgvs_p, csq):
    tx, ref, _, _ = _mini(tmp_path, MINI_TX, [8, 40, len(MINI_TX) - 48], (11, 37), strand)
    assert tx.cds_seq(ref) == MINI_CDS
    e = _edit_at_c(tx, ref, c_from, c_to, alt)
    left = shift(e, ref, -1)
    pos, r, a = to_vcf(left, ref)
    v = Variant(tx.chrom, pos, r, a, left, e, mut_type=classify(left, ref))
    Annotator(tx, ref).annotate(v)
    assert v.hgvs_c == hgvs_c
    assert v.hgvs_p == hgvs_p
    assert v.consequence.split("&")[0] == csq or csq in v.consequence.split("&")


def test_codon_mode_counts(tmp_path):
    tx, ref, fa, ann = _mini(tmp_path, MINI_TX, [8, 40, len(MINI_TX) - 48], (11, 37), "-")
    for cset, per_codon in (("all", 63), ("aa", None), ("nnk", None)):
        out = tmp_path / f"cod_{cset}"
        main(["generate", "-q", "-f", fa, "-a", ann, "-t", "NM_MINI.1", "-r", "c.4_c.24",
              "-m", "codon", "--codon-set", cset, "-o", str(out)])
        rows = [l.split("\t") for l in open(str(out) + ".tsv") if not l.startswith("set")]
        if cset == "all":
            assert len(rows) == 7 * 63
        if cset == "aa":
            # Lys, Glu x3, Trp(no synonym), Pro x2 -> 21 each except Trp (20)
            assert len(rows) == 6 * 21 + 20
            assert any(r[10] == "p.Trp6Ter" for r in rows)
            assert any(r[10] == "p.Glu3=" for r in rows)
        if cset == "nnk":
            assert len(rows) == sum(32 - (c[2] in "GT") for c in
                                    ["AAA", "GAA", "GAA", "GAA", "TGG", "CCC", "CCC"])


def test_dedupe_homopolymer_deletions(tmp_path):
    tx, ref, fa, ann = _mini(tmp_path, MINI_TX, [8, 40, len(MINI_TX) - 48], (11, 37), "+")
    out = tmp_path / "d"
    # the CCCCCC run of codons 7-8 (c.19_c.24)
    main(["generate", "-q", "-f", fa, "-a", ann, "-t", "NM_MINI.1", "-r", "c.19_c.24",
          "-m", "del", "-o", str(out)])
    rows = [l.split("\t") for l in open(str(out) + ".vcf") if not l.startswith("#")]
    assert len(rows) == 1                 # six 1-bp deletions of C -> one variant
    assert rows[0][2] == "c.24del"        # 3'-most on the transcript


def test_regions_and_gtf(tmp_path, synth_env, synth_ref):
    tx = get_tx(synth_env, "NM_MINUS.2", synth_ref)
    r = parse_region("x=exon2:5:7", tx, synth_ref, 1)
    s, e = tx.exons[1]
    assert (r.start, r.end) == (s - 7, e + 5)        # minus strand: 5' flank is to the right
    r2 = parse_region("intron1", tx, synth_ref, 2)
    assert r2.hgvs.startswith("c.30+1_")
    r3 = parse_region(f"{tx.chrom}:{s}-{e}", tx, synth_ref, 3)
    assert (r3.start, r3.end) == (s, e)
    # GTF built from the genePred must give the identical model
    gtf = tmp_path / "x.gtf"
    with open(gtf, "w") as out:
        for st, en in tx.exons_g:
            out.write(f'{tx.chrom}\tt\texon\t{st}\t{en}\t.\t-\t.\tgene_id "MINUSG"; '
                      f'transcript_id "NM_MINUS.2"; gene_name "MINUSG";\n')
        out.write(f'{tx.chrom}\tt\tCDS\t{tx.cds_start + 3}\t{tx.cds_end}\t.\t-\t0\t'
                  f'gene_id "MINUSG"; transcript_id "NM_MINUS.2";\n')
        out.write(f'{tx.chrom}\tt\tstop_codon\t{tx.cds_start}\t{tx.cds_start + 2}\t.\t-\t0\t'
                  f'gene_id "MINUSG"; transcript_id "NM_MINUS.2";\n')
    t2 = load_gtf(str(gtf), "NM_MINUS.2")[0]
    assert t2.exons_g == tx.exons_g and (t2.cds_start, t2.cds_end) == (tx.cds_start, tx.cds_end)
    assert len(load_transcripts(str(gtf), "NM_MINUS")) == 1


def test_version_and_gene_lookup(synth_env, synth_ref):
    hits = load_genepred(synth_env["annotation"], "NM_PLUS.9")
    tx, notes = choose_transcript(hits, "NM_PLUS.9", ref=synth_ref)
    assert tx.name == "NM_PLUS.1" and notes
    hits = load_genepred(synth_env["annotation"], "PLUSG")
    tx, notes = choose_transcript(hits, "PLUSG", ref=synth_ref)
    assert tx.name == "NM_PLUS.1"


def test_chrom_alias_and_fai(tmp_path):
    fa = tmp_path / "g.fa"
    write_fasta(fa, [("16", "ACGT" * 100), ("MT", "A" * 50)], width=37)
    ref = FastaReference(str(fa))
    assert ref.resolve("chr16") == "16" and ref.resolve("chrM") == "MT"
    assert ref.fetch("16", 1, 8) == "ACGTACGT"
    assert ref.fetch("16", 36, 40) == ("ACGT" * 100)[35:40]
    assert (tmp_path / "g.fa.fai").exists()
    ref2 = FastaReference(str(fa))   # reuse index
    assert ref2.fetch("16", 395, 400) == ("ACGT" * 100)[394:400]


def test_exonic_only(tmp_path, synth_env):
    out_a, out_b = tmp_path / "a", tmp_path / "b"
    base = ["generate", "-q", "-f", synth_env["fasta"], "-a", synth_env["annotation"],
            "-t", "NM_MINUS.2", "-r", "c.1_c.300"]
    main(base + ["-o", str(out_a)])
    main(base + ["--exonic-only", "-o", str(out_b)])
    n_all = sum(1 for l in open(str(out_a) + ".vcf") if not l.startswith("#"))
    n_exo = sum(1 for l in open(str(out_b) + ".vcf") if not l.startswith("#"))
    assert n_exo == 300 * 3
    assert n_all > n_exo


def test_cds_keyword_and_intron_note(tmp_path, synth_env, capsys):
    base = ["generate", "-f", synth_env["fasta"], "-a", synth_env["annotation"], "-t", "NM_PLUS.1"]
    main(base + ["-r", "cds", "--exonic-only", "-o", str(tmp_path / "c")])
    n = sum(1 for l in open(str(tmp_path / "c.vcf")) if not l.startswith("#"))
    assert n == 300 * 3
    capsys.readouterr()
    main(base + ["-r", "c.30+5_c.30+40", "-o", str(tmp_path / "i")])   # inside intron 1 only
    assert "spans intron" not in capsys.readouterr().err
    main(base + ["-r", "c.20_c.40", "-o", str(tmp_path / "x")])        # exon1 .. exon2
    assert "spans intron(s) 1" in capsys.readouterr().err
