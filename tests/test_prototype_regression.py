"""VCFART must reproduce the loop logic of the authors' original MAVE scripts.

The original scripts walked a minus-strand intron / exon / intron block base by base,
wrote every SNV to a VCF and labelled it with a complemented c. description. Here the
same loops run on a synthetic block (see conftest.py) and must give exactly the
records and labels VCFART produces.
"""

from conftest import PROTO_EXON, PROTO_INTRON_DOWN, PROTO_INTRON_UP

from vcfart.cli import main
from vcfart.sequtils import complement


def prototype_records():
    """Re-implementation of the original loops (with complemented c. bases everywhere)."""
    exp = {}
    for i, b in enumerate(PROTO_INTRON_UP):
        for alt in "ATCG":
            if alt != b:
                exp[(10294 - i, b, alt)] = f"c.301-{132 - i}{complement(b)}>{complement(alt)}"
    for i, b in enumerate(PROTO_EXON):
        for alt in "ATCG":
            if alt != b:
                exp[(10162 - i, b, alt)] = f"c.{301 + i}{complement(b)}>{complement(alt)}"
    for i, b in enumerate(PROTO_INTRON_DOWN):
        for alt in "ATCG":
            if alt != b:
                exp[(10000 - i, b, alt)] = f"c.462+{i + 1}{complement(b)}>{complement(alt)}"
    return exp


def read_vcf(path):
    recs = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in f[7].split(";"))
            recs[(int(f[1]), f[3], f[4])] = (f[2], info, f)
    return recs


def test_reproduces_prototype(proto_env, tmp_path):
    out = tmp_path / "proto"
    main(["generate", "-q", "-f", proto_env["fasta"], "-a", proto_env["annotation"],
          "-b", "hg38", "-t", "NM_PROTO.1",
          "-r", "up=c.301-132_c.301-1", "-r", "ex=c.301_c.462",
          "-r", "down=c.462+1_c.462+150", "-o", str(out)])
    got = read_vcf(str(out) + ".vcf")
    exp = prototype_records()
    assert len(exp) == (132 + 162 + 150) * 3
    assert set(got) == set(exp)
    for key, label in exp.items():
        vid, info, _ = got[key]
        assert vid == label, (key, vid, label)
        assert info["VA_HGVSc"] == label


def test_exon_shortcut_equals_explicit_ranges(proto_env, tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    common = ["generate", "-q", "-f", proto_env["fasta"], "-a", proto_env["annotation"],
              "-t", "NM_PROTO"]
    main(common + ["-r", "exon2:132:150", "-o", str(a)])
    main(common + ["-r", "c.301-132_c.462+150", "-o", str(b)])
    ra, rb = read_vcf(str(a) + ".vcf"), read_vcf(str(b) + ".vcf")
    assert set(ra) == set(rb)
    assert len(ra) == 444 * 3


def test_splice_site_annotations(proto_env, tmp_path):
    out = tmp_path / "p"
    main(["generate", "-q", "-f", proto_env["fasta"], "-a", proto_env["annotation"],
          "-t", "NM_PROTO.1", "-r", "exon2:2:2", "-o", str(out)])
    by_c = {}
    for vid, info, _ in read_vcf(str(out) + ".vcf").values():
        by_c.setdefault(vid[:-3], []).append(info)       # strip 'X>Y'
    assert all(i["VA_CSQ"].startswith("splice_acceptor_variant") for i in by_c["c.301-1"])
    assert all(i["VA_REGION"] == "intron1" for i in by_c["c.301-2"])
    assert all(i["VA_CSQ"].startswith("splice_donor_variant") for i in by_c["c.462+1"])
    assert all(i["VA_REGION"] == "intron2" for i in by_c["c.462+2"])
    assert all(i["VA_REGION"] == "exon2" for i in by_c["c.301"])
    assert all("splice_region_variant" in i["VA_CSQ"] for i in by_c["c.301"])


def test_oligo_matches_prototype_full_sequence(proto_env, tmp_path):
    """--oligo reproduces adapter + mutated intron + adapter in transcript sense."""
    out = tmp_path / "o"
    main(["generate", "-q", "-f", proto_env["fasta"], "-a", proto_env["annotation"],
          "-t", "NM_PROTO.1", "-r", "up=c.301-132_c.301-1", "--oligo",
          "--adapter5", "AAAA", "--adapter3", "CCCC", "-o", str(out)])
    seqs = {}
    with open(str(out) + ".oligos.fa") as fh:
        name = None
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                seqs[name] = line.strip()
    assert seqs["up|WT"] == "AAAA" + complement(PROTO_INTRON_UP) + "CCCC"
    b0 = PROTO_INTRON_UP[0]                      # plus-strand base at c.301-132
    alt = "A" if b0 != "A" else "C"
    mut = alt + PROTO_INTRON_UP[1:]
    label = f"c.301-132{complement(b0)}>{complement(alt)}"
    assert seqs[f"up|{label}"] == "AAAA" + complement(mut) + "CCCC"
    assert len(seqs) == 1 + 132 * 3
