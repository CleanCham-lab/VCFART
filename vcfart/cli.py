"""Command-line interface.

  vcfart generate  -f hg38.fa -a hg38_refGene.txt -t NM_007294.4 -r exon17:100 -o BRCA1_ex17
  vcfart info      -a hg38_refGene.txt -t NM_007294.4 [-f hg38.fa]
  vcfart fetch-annotation -b hg38
"""

import argparse
import os
import shlex
import signal
import sys
from collections import Counter

from . import __version__, utils
from .annotation import choose_transcript, load_transcripts
from .fetch import TRACKS, cached_annotation, fetch
from .hgvs import Annotator
from .reference import open_reference
from .regions import parse_region
from .sequtils import complement, revcomp, translate
from .utils import (AMBER, BOLD, DIM, GREEN, GREY, LILAC, PURPLE, SEVERITY_COLOR, die, log,
                    paint, warn)
from .variants import (Variant, apply_edit, classify, gen_codon, gen_del, gen_ins, gen_snv,
                       shift, to_vcf)
from .writers import write_fasta, write_tsv, write_vcf

BANNER_ART = [
    " ██╗   ██╗ ██████╗███████╗ █████╗ ██████╗ ████████╗",
    " ██║   ██║██╔════╝██╔════╝██╔══██╗██╔══██╗╚══██╔══╝",
    " ██║   ██║██║     █████╗  ███████║██████╔╝   ██║",
    " ╚██╗ ██╔╝██║     ██╔══╝  ██╔══██║██╔══██╗   ██║",
    "  ╚████╔╝ ╚██████╗██║     ██║  ██║██║  ██║   ██║",
    "   ╚═══╝   ╚═════╝╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝   ╚═╝",
]
BANNER_SHADES = ["38;5;141", "38;5;141", "38;5;135", "38;5;98", "38;5;97", "38;5;61"]


def banner():
    lines = [""] + [paint(t, c) for t, c in zip(BANNER_ART, BANNER_SHADES)]
    lines.append("   " + paint("💨 뿌뿡", BOLD, LILAC) + paint("  artificial VCFs for MAVE / SGE / DMS", GREY)
                 + "   " + paint(f"v{__version__}", PURPLE))
    return "\n".join(lines) + "\n"


def _bar_chart(counter, tag):
    """Consequence counts as a coloured horizontal bar chart (most common first)."""
    items = counter.most_common()
    top = items[0][1] if items else 1
    width = max(len(k) for k, _ in items) if items else 0
    for i, (k, n) in enumerate(items):
        col = SEVERITY_COLOR.get(k, GREY)
        bar = "█" * max(1, round(22 * n / top))
        log(f"{'📊 ' if i == 0 else '   '}{k.ljust(width)} {paint(bar, col)} {n:,}", tag if i == 0 else "")

MUT_TYPES = ("snv", "del", "ins", "codon")


def _sizes(text, name):
    try:
        vals = sorted({int(x) for x in text.split(",") if x.strip()})
    except ValueError:
        die(f"--{name} expects comma-separated integers, got '{text}'")
    if not vals or min(vals) < 1:
        die(f"--{name} values must be >= 1")
    return vals


def _resolve_annotation(args):
    if args.annotation:
        if not os.path.exists(args.annotation):
            die(f"annotation not found: {args.annotation}")
        return args.annotation
    if args.build:
        p = cached_annotation(args.build)
        if p:
            log(f"📚 using cached {p}", "annotation")
            return p
        die(f"no annotation given and none cached for {args.build}. "
            f"Run `vcfart fetch-annotation -b {args.build}` or pass -a "
            "(UCSC refGene/ncbiRefSeq table, ANNOVAR humandb/*_refGene.txt, or GTF).")
    die("give -a/--annotation (or -b/--build with a cached annotation)")


def _load_tx(args, ref=None):
    ann = _resolve_annotation(args)
    log(f"📚 {ann}", "annotation")
    hits = load_transcripts(ann, args.transcript)
    try:
        tx, notes = choose_transcript(hits, args.transcript, ref=ref, chrom=args.chrom)
    except LookupError as e:
        die(str(e))
    for n in notes:
        if n.startswith("gene symbol"):
            log(n, "transcript")
        else:
            warn(n)
    if ref is not None:
        tx.chrom = ref.resolve(tx.chrom)
    return tx


def _sanity(tx, ref):
    if not tx.coding:
        return
    cds = tx.cds_seq(ref)
    if "N" in cds:
        warn("CDS contains N bases in the reference")
    if len(cds) % 3:
        warn(f"CDS length {len(cds)} is not a multiple of 3 (incomplete CDS?) -- p. may be off")
    if cds[:3] != "ATG":
        warn(f"CDS does not start with ATG ({cds[:3]}) -- check genome build vs annotation!")
    if translate(cds[-3:]) != "*":
        warn(f"CDS does not end with a stop codon ({cds[-3:]}) -- check build/annotation")
    motifs = []
    for k in range(tx.n_exons - 1):
        last, first = tx.exon_n[k][1], tx.exon_n[k + 1][0]
        donor = "".join(ref.base(tx.chrom, tx.n_to_g(last, o)) for o in (1, 2))
        acc = "".join(ref.base(tx.chrom, tx.n_to_g(first, o)) for o in (-2, -1))
        if not tx.plus:
            donor, acc = complement(donor), complement(acc)
        motifs.append(donor in ("GT", "GC") and acc == "AG")
    if motifs and sum(motifs) / len(motifs) < 0.5:
        warn(f"only {sum(motifs)}/{len(motifs)} introns have GT-AG/GC-AG splice motifs -- "
             "FASTA and annotation builds probably do not match (hg19 vs hg38?)")


def _intron_inside(tx, k, r):
    """True if the region reaches exonic bases on BOTH sides of intron k (c. range across it)."""
    def overlaps(e):
        return r.start <= e[1] and e[0] <= r.end
    return overlaps(tx.exons[k - 1]) and overlaps(tx.exons[k])


def _exonic_edit(tx, edit):
    """True if the edit touches exonic bases of a single exon only (cDNA-style saturation)."""
    gs = [edit.start - 1, edit.start] if edit.ref_len == 0 else [edit.start, edit.end]
    regs = {tx.region(g) for g in gs}
    return len(regs) == 1 and next(iter(regs))[0] == "exon"


# ===================================================================== generate
def cmd_generate(args):
    ref = open_reference(args.fasta)
    log(f"🧬 {ref.path}", "reference")
    tx = _load_tx(args, ref)
    chrom = tx.chrom
    log("📍 " + paint(f"{tx.name} {tx.gene}", BOLD) + f"  {tx.chrom}:{tx.tx_start:,}-{tx.tx_end:,} ({tx.strand})"
        f"  {tx.n_exons} exons" + (f"  CDS {tx.cds_len:,} bp" if tx.coding else "  non-coding"),
        "transcript")
    _sanity(tx, ref)

    types = [t.strip().lower() for t in args.mut.split(",") if t.strip()]
    bad = [t for t in types if t not in MUT_TYPES]
    if bad:
        die(f"unknown mutation type(s) {bad}; choose from {', '.join(MUT_TYPES)}")
    if "codon" in types and not tx.coding:
        die("codon mode needs a coding transcript")
    del_sizes = _sizes(args.del_sizes, "del-sizes")
    ins_sizes = _sizes(args.ins_sizes, "ins-sizes")

    regions = []
    for i, spec in enumerate(args.region, 1):
        try:
            regions.append(parse_region(spec, tx, ref, i))
        except ValueError as e:
            die(str(e))
    labels = [r.label for r in regions]
    if len(set(labels)) != len(labels):
        die(f"region labels must be unique: {labels}")

    lw = max(len(r.label) for r in regions)
    hw = max(len(r.hgvs) for r in regions)
    for r in regions:
        a, b = (r.start, r.end) if tx.plus else (r.end, r.start)
        locs = [tx.region_label(a), tx.region_label(b)]
        span = locs[0] if locs[0] == locs[1] else f"{locs[0]} .. {locs[1]}"
        log("🎯 " + paint(r.label.ljust(lw), BOLD, PURPLE) + f"  {r.hgvs.ljust(hw)}  {r.chrom}:{r.start:,}-{r.end:,}  "
            f"{r.length:,} bp  " + paint(span.replace(" .. ", " → "), GREY), "region")
        if r.length > 100000:
            warn(f"region {r.label} is {r.length:,} bp -- this will be a big VCF")
        inner = [k for k in range(1, tx.n_exons) if _intron_inside(tx, k, r)]
        if inner and not args.exonic_only:
            log(paint(f"   ↳ spans intron(s) {', '.join(map(str, inner))}, included "
                      "(--exonic-only for cDNA positions only)", GREY), "")

    variants = {}
    per_set = {r.label: [] for r in regions}
    stats = {}
    raw_count = 0
    for r in regions:
        gens = []
        if "snv" in types:
            gens.append(gen_snv(r, ref))
        if "del" in types:
            gens.append(gen_del(r, ref, del_sizes))
        if "ins" in types:
            gens.append(gen_ins(r, ref, ins_sizes))
        if "codon" in types:
            gens.append(gen_codon(r, ref, tx, args.codon_set, args.codon_pick, stats))
        seen_in_set = set()
        for g in gens:
            for edit, origin in g:
                if args.exonic_only and not _exonic_edit(tx, edit):
                    continue
                raw_count += 1
                norm = shift(edit, ref, -1)
                if norm.ref_len == 0 and not norm.alt:
                    continue
                pos, vref, valt = to_vcf(norm, ref)
                key = (chrom, pos, vref, valt)
                v = variants.get(key)
                if v is None:
                    v = Variant(chrom, pos, vref, valt, norm, edit, mut_type=classify(norm, ref))
                    variants[key] = v
                if origin not in v.origins:
                    v.origins.append(origin)
                if r.label not in v.sets:
                    v.sets.append(r.label)
                    v.raws[r.label] = edit
                if key not in seen_in_set:
                    seen_in_set.add(key)
                    per_set[r.label].append(v)

    if stats.get("split_codons"):
        warn(f"{stats['split_codons']} codon(s) split by an intron were skipped in codon mode")
    if not variants:
        die("no variants generated (empty region or all-N reference?)")

    annot = Annotator(tx, ref)
    for v in variants.values():
        annot.annotate(v)

    ordered = sorted(variants.values(), key=lambda v: (v.pos, len(v.ref), v.ref, v.alt))
    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)
    vcf_path = args.out + ".vcf"
    tsv_path = args.out + ".tsv"
    command = "vcfart " + " ".join(shlex.quote(a) for a in sys.argv[1:])
    write_vcf(vcf_path, ordered, tx, ref, args.build or "", command,
              sites_only=args.sites_only, gt=args.gt, sample=args.sample_name,
              id_style=args.id, all_contigs=args.all_contigs)

    rows = []
    fasta_records = []
    for r in regions:
        olig = {}
        if args.oligo:
            ws, we = r.start - args.oligo_pad, r.end + args.oligo_pad
            wt = ref.fetch(chrom, ws, we)

            def orient(seq):
                s = seq if (args.oligo_strand == "plus" or tx.plus) else revcomp(seq)
                return args.adapter5 + s + args.adapter3

            fasta_records.append((f"{r.label}|WT", orient(wt)))
            for v in per_set[r.label]:
                m = apply_edit(wt, ws, v.raws.get(r.label, v.raw))
                if m is None:
                    m = apply_edit(wt, ws, v.edit)
                if m is None:
                    continue
                oid = f"{r.label}|{v.hgvs_c}"
                olig[v.key] = (oid, orient(m))
                fasta_records.append(olig[v.key])
        for v in sorted(per_set[r.label], key=lambda v: (v.pos, len(v.ref), v.ref, v.alt)):
            rows.append((r.label, v, olig.get(v.key)))
    write_tsv(tsv_path, rows, tx, oligo=args.oligo)
    outputs = [vcf_path, tsv_path]
    if args.oligo:
        fa_path = args.out + ".oligos.fa"
        write_fasta(fa_path, fasta_records)
        outputs.append(fa_path)

    by_type = Counter(v.mut_type for v in ordered)
    by_csq = Counter(v.consequence.split("&")[0] if v.consequence else "intergenic"
                     for v in ordered)
    type_col = {"SNV": LILAC, "MNV": PURPLE, "DEL": AMBER, "INS": GREEN, "DELINS": AMBER}
    log(f"💨 {raw_count:,} raw edits → " + paint(f"{len(ordered):,} unique variants", BOLD), "generate")
    log("🧪 " + "   ".join(paint(k, BOLD, type_col.get(k, GREY)) + f" {n:,}"
                          for k, n in by_type.most_common()), "types")
    _bar_chart(by_csq, "consequence")
    for i, p in enumerate(outputs):
        log(f"💾 {p}", "write" if i == 0 else "")
    if not utils.QUIET:
        print("\n  " + paint("💨 뿌뿡!", BOLD, LILAC) + "  " + paint(f"{len(ordered):,} variants", BOLD)
              + paint(" → ", GREY) + vcf_path + "\n", file=sys.stderr)


# ========================================================================= info
def cmd_info(args):
    ref = open_reference(args.fasta) if args.fasta else None
    tx = _load_tx(args, ref)
    tty = utils.COLOR and sys.stdout.isatty()

    def out(text, *codes):
        print(paint(text, *codes) if tty else text)
    out(tx.describe(), BOLD)
    if tx.coding:
        out(f"CDS: c.1 = {tx.chrom}:{tx.n_to_g(tx.cds_start_n)}, "
            f"c.{tx.cds_len} (stop end) = {tx.chrom}:{tx.n_to_g(tx.cds_end_n)}, "
            f"{tx.cds_len // 3 - 1} aa", GREY)
    print()
    out("feature\tgenomic_start\tgenomic_end\tlength\tfrom\tto", GREY)
    rows = tx.exon_table()
    for i, row in enumerate(rows):
        k, s, e, length, a, b = row
        out(f"exon{k}\t{s}\t{e}\t{length}\t{a}\t{b}", BOLD, PURPLE)
        if i < len(rows) - 1:
            nk = rows[i + 1]
            lo = min(e, nk[2]) + 1
            hi = max(s, nk[1]) - 1
            ga, gb = (lo, hi) if tx.plus else (hi, lo)
            out(f"intron{k}\t{lo}\t{hi}\t{hi - lo + 1}\t{tx.prefix}{tx.pos(ga)}\t"
                f"{tx.prefix}{tx.pos(gb)}", DIM)
    if ref is not None:
        utils.QUIET = False
        _sanity(tx, ref)


def cmd_fetch(args):
    try:
        fetch(args.build, args.track, args.dir, args.insecure)
    except (RuntimeError, ValueError) as e:
        die(str(e))


# ========================================================================= main
def build_parser():
    p = argparse.ArgumentParser(
        prog="vcfart",
        description="VCFART: generate artificial saturation-mutagenesis VCFs from a transcript "
                    "and an HGVS c./g. range (MAVE / SGE / DMS / SpliceAI scans).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="examples:\n"
               "  vcfart generate -f hg38.fa -a hg38_refGene.txt -t NM_007294.4 \\\n"
               "      -r int16=c.5075-100_c.5075-1 -r ex17=c.5075_c.5152 "
               "-r int17=c.5152+1_c.5152+100 -o BRCA1_ex17\n"
               "  vcfart generate -f hg38.fa -b hg38 -t NM_000546.6 -r exon5:20 "
               "-m snv,del,codon -o TP53_ex5\n"
               "  vcfart info -a hg38_refGene.txt -t NM_007294")
    p.add_argument("-V", "--version", action="version", version=f"vcfart {__version__}")
    sub = p.add_subparsers(dest="cmd")

    def common(sp, need_fasta):
        sp.add_argument("-f", "--fasta", required=need_fasta,
                        help="reference FASTA (hg19/hg38; .fai used or built)")
        sp.add_argument("-a", "--annotation",
                        help="UCSC refGene/ncbiRefSeq table (ANNOVAR humandb/*_refGene.txt ok) or GTF")
        sp.add_argument("-b", "--build", help="genome build label for the VCF header "
                                              "(hg19/hg38); also selects cached annotation")
        sp.add_argument("-t", "--transcript", required=True,
                        help="transcript accession, e.g. NM_007294.4 (version optional) or ENST...")
        sp.add_argument("--chrom", help="pick the copy on this chromosome if the transcript maps "
                                        "to several loci")
        sp.add_argument("-q", "--quiet", action="store_true", help="no banner / progress")
        sp.add_argument("--no-color", action="store_true",
                        help="plain output (also: NO_COLOR=1; colour is off when piped)")

    g = sub.add_parser("generate", aliases=["gen", "fart"], help="generate the variant VCF",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    common(g, True)
    g.add_argument("-r", "--region", action="append", required=True,
                   help="target range, repeatable, optional LABEL= prefix: c.5075-100_c.5075-1 | "
                        "g.43063774_43064051 | chr17:43063774-43064051 | exon17[:5flank[:3flank]] "
                        "| intron16 | cds")
    g.add_argument("-o", "--out", required=True, help="output prefix (writes PREFIX.vcf, PREFIX.tsv)")
    g.add_argument("-m", "--mut", default="snv",
                   help="mutation types, comma-separated: snv,del,ins,codon")
    g.add_argument("--del-sizes", default="1", help="deletion lengths for -m del, e.g. 1,2,3")
    g.add_argument("--ins-sizes", default="1",
                   help="insertion lengths for -m ins (all 4^k sequences per gap)")
    g.add_argument("--codon-set", choices=["aa", "all", "nnk"], default="aa",
                   help="codon mode: one codon per amino acid (+Ter, +synonymous), all 63, or NNK")
    g.add_argument("--codon-pick", choices=["usage", "min-edit"], default="usage",
                   help="codon choice for --codon-set aa: human preferred codon or fewest edits")
    g.add_argument("--exonic-only", action="store_true",
                   help="only exonic positions (cDNA-style); default keeps introns inside the range")
    g.add_argument("--sites-only", action="store_true", help="omit FORMAT/sample columns")
    g.add_argument("--gt", default="0/1", help="genotype written in the dummy sample column")
    g.add_argument("--sample-name", default="VCFART", help="dummy sample column name")
    g.add_argument("--id", choices=["hgvs", "index", "none"], default="hgvs",
                   help="what to put in the VCF ID column")
    g.add_argument("--all-contigs", action="store_true",
                   help="write every FASTA contig to the VCF header (default: only the used one)")
    g.add_argument("--oligo", action="store_true",
                   help="also write mutant sequences (TSV columns + PREFIX.oligos.fa)")
    g.add_argument("--oligo-pad", type=int, default=0,
                   help="extend each region by N bp on both sides for the oligo sequence")
    g.add_argument("--oligo-strand", choices=["transcript", "plus"], default="transcript",
                   help="orientation of oligo sequences")
    g.add_argument("--adapter5", default="", help="constant sequence prepended to each oligo")
    g.add_argument("--adapter3", default="", help="constant sequence appended to each oligo")
    g.set_defaults(func=cmd_generate)

    i = sub.add_parser("info", help="show exon/intron table with c. and genomic coordinates")
    common(i, False)
    i.set_defaults(func=cmd_info)

    f = sub.add_parser("fetch-annotation", help="download a UCSC transcript table into ~/.vcfart")
    f.add_argument("-b", "--build", required=True, help="hg19 or hg38")
    f.add_argument("--track", choices=TRACKS, default="ncbiRefSeq")
    f.add_argument("--dir", help="output directory (default ~/.vcfart or $VCFART_HOME)")
    f.add_argument("--insecure", action="store_true",
                   help="skip TLS certificate verification (SSL-inspecting proxies)")
    f.add_argument("-q", "--quiet", action="store_true")
    f.add_argument("--no-color", action="store_true")
    f.set_defaults(func=cmd_fetch)
    return p


def main(argv=None):
    if hasattr(signal, "SIGPIPE"):
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "cmd", None):
        parser.print_help()
        return 0
    utils.QUIET = bool(getattr(args, "quiet", False)) or args.cmd == "info"
    utils.COLOR = utils.color_wanted(getattr(args, "no_color", False))
    if args.cmd in ("generate", "gen", "fart") and not utils.QUIET:
        print(banner(), file=sys.stderr)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
