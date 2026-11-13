"""Transcript models and HGVS c./n. <-> genomic coordinate conversion.

Supported annotation inputs
---------------------------
* UCSC genePred tables: ``refGene.txt(.gz)``, ``ncbiRefSeq.txt(.gz)``,
  ``ncbiRefSeqSelect.txt(.gz)`` -- with or without the leading ``bin`` column.
  ANNOVAR ``humandb/hg38_refGene.txt`` / ``hg38_refGeneWithVer.txt`` /
  ``hg19_refGene.txt`` are the same format and work as-is.
* GTF: UCSC ``hg38.ncbiRefSeq.gtf.gz``, GENCODE / Ensembl GTF (ENST IDs).
"""

import gzip
import re

from .sequtils import revcomp
from .utils import warn

_ACC_RE = re.compile(r"^(?P<acc>[A-Za-z]+_?\d+)(?:\.(?P<ver>\d+))?(?:_(?P<dup>\d+))?$")


def split_accession(name):
    m = _ACC_RE.match(name)
    if m:
        return m.group("acc"), m.group("ver")
    head, dot, tail = name.rpartition(".")
    if dot and tail.isdigit():
        return head, tail
    return name, None


class Transcript:
    """A spliced transcript on the genome.

    Coordinates: genomic positions are 1-based inclusive. ``n`` is the
    1-based position along the mature transcript (exons only, 5'->3').
    """

    def __init__(self, name, chrom, strand, exons, cds_start=None, cds_end=None,
                 gene="", source=""):
        if strand not in "+-" or not exons:
            raise ValueError(f"bad transcript {name}")
        self.name = name
        self.gene = gene or name
        self.chrom = chrom
        self.strand = strand
        self.plus = strand == "+"
        self.source = source
        self.exons_g = sorted(exons)
        self.tx_start = self.exons_g[0][0]
        self.tx_end = self.exons_g[-1][1]
        self.exons = self.exons_g if self.plus else list(reversed(self.exons_g))
        self.exon_n = []
        n = 1
        for s, e in self.exons:
            self.exon_n.append((n, n + (e - s)))
            n += e - s + 1
        self.length = n - 1
        self.coding = cds_start is not None and cds_end is not None and cds_end >= cds_start
        self.cds_start = cds_start if self.coding else None
        self.cds_end = cds_end if self.coding else None
        if self.coding:
            first = cds_start if self.plus else cds_end
            last = cds_end if self.plus else cds_start
            self.cds_start_n = self._exonic_n(first)
            self.cds_end_n = self._exonic_n(last)
            if self.cds_start_n is None or self.cds_end_n is None:
                raise ValueError(f"{name}: CDS boundary not inside an exon")
            self.cds_len = self.cds_end_n - self.cds_start_n + 1
        self._tx_seq = None

    # --------------------------------------------------------------- helpers
    @property
    def prefix(self):
        return "c." if self.coding else "n."

    @property
    def n_exons(self):
        return len(self.exons)

    def _exonic_n(self, g):
        for i, (s, e) in enumerate(self.exons):
            if s <= g <= e:
                n0 = self.exon_n[i][0]
                return n0 + (g - s if self.plus else e - g)
        return None

    def g_to_n(self, g):
        """genomic -> (n, offset). offset != 0 only for intronic positions.

        Upstream of the transcript gives n <= 0, downstream gives n > length.
        """
        n = self._exonic_n(g)
        if n is not None:
            return n, 0
        if self.plus:
            if g < self.tx_start:
                return 1 - (self.tx_start - g), 0
            if g > self.tx_end:
                return self.length + (g - self.tx_end), 0
        else:
            if g > self.tx_end:
                return 1 - (g - self.tx_end), 0
            if g < self.tx_start:
                return self.length + (self.tx_start - g), 0
        for k in range(len(self.exons) - 1):
            s1, e1 = self.exons[k]
            s2, e2 = self.exons[k + 1]
            if self.plus and e1 < g < s2:
                d5, d3 = g - e1, s2 - g
            elif not self.plus and e2 < g < s1:
                d5, d3 = s1 - g, g - e2
            else:
                continue
            if d5 <= d3:
                return self.exon_n[k][1], d5
            return self.exon_n[k + 1][0], -d3
        raise ValueError(f"cannot place {g} on {self.name}")

    def n_to_g(self, n, offset=0):
        if 1 <= n <= self.length:
            for i, (n0, n1) in enumerate(self.exon_n):
                if n0 <= n <= n1:
                    s, e = self.exons[i]
                    g = s + (n - n0) if self.plus else e - (n - n0)
                    break
        elif n < 1:
            g = self.tx_start - (1 - n) if self.plus else self.tx_end + (1 - n)
        else:
            g = self.tx_end + (n - self.length) if self.plus else self.tx_start - (n - self.length)
        return g + offset if self.plus else g - offset

    def n_to_pos(self, n):
        """transcript position -> c./n. numbering string (no prefix, no offset)."""
        if self.coding:
            if n < self.cds_start_n:
                return f"-{self.cds_start_n - n}"
            if n <= self.cds_end_n:
                return str(n - self.cds_start_n + 1)
            return f"*{n - self.cds_end_n}"
        if n < 1:
            return f"-{1 - n}"
        if n <= self.length:
            return str(n)
        return f"*{n - self.length}"

    def pos(self, g):
        """genomic position -> HGVS position string, e.g. '5075-1', '*12', '-30+2'."""
        n, off = self.g_to_n(g)
        s = self.n_to_pos(n)
        if off > 0:
            s += f"+{off}"
        elif off < 0:
            s += str(off)
        return s

    _POS_RE = re.compile(r"^(?:([cngr])\.)?([-*]?)(\d+)(?:([+-])(\d+))?$")

    def parse_pos(self, text):
        """HGVS position ('c.5075-1', '-30', '*5+1', 'n.12') -> genomic position."""
        t = text.strip().replace(" ", "")
        m = self._POS_RE.match(t)
        if not m:
            raise ValueError(f"cannot parse position '{text}'")
        kind, sign, num, osign, onum = m.groups()
        num = int(num)
        if kind == "c" and not self.coding:
            raise ValueError(f"{self.name} is non-coding; use n. positions")
        if kind == "g":
            raise ValueError("use g.START_END for genomic ranges")
        if self.coding and kind != "n":
            if sign == "" and num == 0:
                raise ValueError("c.0 does not exist (c.-1 is followed by c.1)")
            if sign == "-":
                n = self.cds_start_n - num
            elif sign == "*":
                n = self.cds_end_n + num
            else:
                n = self.cds_start_n + num - 1
        else:
            if sign == "-":
                n = 1 - num
            elif sign == "*":
                n = self.length + num
            else:
                n = num
        off = 0
        if osign:
            off = int(onum) * (1 if osign == "+" else -1)
            ok = any((off > 0 and n == n1) or (off < 0 and n == n0)
                     for n0, n1 in self.exon_n)
            if not ok:
                warn(f"'{text}': intronic offset is not anchored on an exon boundary "
                     "(non-standard HGVS) -- computed anyway")
        return self.n_to_g(n, off)

    def region(self, g):
        """('exon', k) / ('intron', k) / ('upstream', 0) / ('downstream', 0); k is 1-based."""
        for i, (s, e) in enumerate(self.exons):
            if s <= g <= e:
                return "exon", i + 1
        if self.plus:
            if g < self.tx_start:
                return "upstream", 0
            if g > self.tx_end:
                return "downstream", 0
        else:
            if g > self.tx_end:
                return "upstream", 0
            if g < self.tx_start:
                return "downstream", 0
        for k in range(len(self.exons) - 1):
            lo = min(self.exons[k][1], self.exons[k + 1][1])
            hi = max(self.exons[k][0], self.exons[k + 1][0])
            if lo < g < hi:
                return "intron", k + 1
        return "intergenic", 0

    def region_label(self, g):
        kind, k = self.region(g)
        return f"{kind}{k}" if k else kind

    # ------------------------------------------------------------- sequence
    def transcript_seq(self, ref):
        if self._tx_seq is None:
            parts = [ref.fetch(self.chrom, s, e) for s, e in self.exons]
            if self.plus:
                self._tx_seq = "".join(parts)
            else:
                self._tx_seq = "".join(revcomp(p) for p in parts)
        return self._tx_seq

    def cds_seq(self, ref):
        if not self.coding:
            return ""
        return self.transcript_seq(ref)[self.cds_start_n - 1:self.cds_end_n]

    def describe(self):
        span = f"{self.chrom}:{self.tx_start}-{self.tx_end}"
        kind = f"coding, CDS {self.cds_len} bp" if self.coding else "non-coding"
        return (f"{self.name} ({self.gene}) {span} strand {self.strand}, "
                f"{self.n_exons} exons, {self.length} nt, {kind}")

    def exon_table(self):
        rows = []
        for i, ((s, e), (n0, n1)) in enumerate(zip(self.exons, self.exon_n)):
            first = self.n_to_pos(n0)
            last = self.n_to_pos(n1)
            rows.append((i + 1, s, e, e - s + 1, f"{self.prefix}{first}", f"{self.prefix}{last}"))
        return rows


# ====================================================================== parsing
def _open(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _name_matches(name, acc, ver):
    a, v = split_accession(name)
    if a != acc:
        return False
    return True


def load_genepred(path, query):
    """Return transcripts whose name (accession) or gene symbol matches ``query``."""
    acc, ver = split_accession(query)
    hits = []
    with _open(path) as fh:
        for line in fh:
            if acc not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            o = 1 if (f[0].isdigit() and len(f) > 3 and f[3] in ("+", "-")) else 0
            name, chrom, strand = f[o], f[o + 1], f[o + 2]
            gene = f[o + 11] if len(f) > o + 11 else ""
            if not (_name_matches(name, acc, ver) or gene == query):
                continue
            starts = [int(x) for x in f[o + 8].strip(",").split(",") if x]
            ends = [int(x) for x in f[o + 9].strip(",").split(",") if x]
            exons = [(s + 1, e) for s, e in zip(starts, ends)]
            cs, ce = int(f[o + 5]), int(f[o + 6])
            cds = (cs + 1, ce) if ce > cs else (None, None)
            hits.append(Transcript(name, chrom, strand, exons, cds[0], cds[1], gene,
                                   source=str(path)))
    return hits


_ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def load_gtf(path, query):
    acc, ver = split_accession(query)
    recs = {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or acc not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feat = f[2]
            if feat not in ("exon", "CDS", "start_codon", "stop_codon"):
                continue
            attrs = dict(_ATTR_RE.findall(f[8]))
            tid = attrs.get("transcript_id", "")
            gene = attrs.get("gene_name") or attrs.get("gene_id", "")
            if not (_name_matches(tid, acc, ver) or gene == query):
                continue
            key = (tid, f[0], f[6])
            r = recs.setdefault(key, {"exons": [], "cds": [], "gene": gene})
            s, e = int(f[3]), int(f[4])
            if feat == "exon":
                r["exons"].append((s, e))
            else:
                r["cds"].append((s, e))
    hits = []
    for (tid, chrom, strand), r in recs.items():
        if not r["exons"]:
            continue
        cds = (min(s for s, _ in r["cds"]), max(e for _, e in r["cds"])) if r["cds"] else (None, None)
        hits.append(Transcript(tid, chrom, strand, sorted(set(r["exons"])), cds[0], cds[1],
                               r["gene"], source=str(path)))
    return hits


def load_transcripts(path, query):
    p = str(path).lower()
    if p.endswith((".gtf", ".gtf.gz", ".gff", ".gff.gz")):
        return load_gtf(path, query)
    return load_genepred(path, query)


def choose_transcript(hits, query, ref=None, chrom=None):
    """Pick one transcript from candidate hits; returns (Transcript, notes)."""
    acc, ver = split_accession(query)
    notes = []
    if not hits:
        raise LookupError(f"'{query}' not found in the annotation file")

    by_acc = [t for t in hits if split_accession(t.name)[0] == acc]
    if not by_acc:
        accs = sorted({split_accession(t.name)[0] for t in hits})
        if len(accs) > 1:
            raise LookupError(
                f"'{query}' is a gene symbol with several transcripts; choose one with -t: "
                + ", ".join(sorted({t.name for t in hits})[:40]))
        notes.append(f"gene symbol {query} -> only transcript {hits[0].name}")
        acc, ver = accs[0], None
    else:
        hits = by_acc

    def _v(t):
        return split_accession(t.name)[1]

    exact = [t for t in hits if ver is not None and _v(t) == ver]
    if exact:
        hits = exact
    else:
        vers = sorted({_v(t) for t in hits}, key=lambda v: int(v or 0))
        if ver is not None:
            if vers == [None]:
                notes.append(f"annotation has no version numbers; assuming {acc} == {query}")
            else:
                have = "/".join(str(v) for v in vers if v)
                notes.append(f"{query} not in annotation (has {acc}.{have}); "
                             f"using {acc}.{vers[-1]}")
        top = vers[-1]
        hits = [t for t in hits if _v(t) == top]

    if chrom:
        hits = [t for t in hits if t.chrom == chrom or t.chrom.replace("chr", "") == chrom.replace("chr", "")]
        if not hits:
            raise LookupError(f"{query} not annotated on {chrom}")

    primary = [t for t in hits if "_" not in t.chrom]
    if primary and len(primary) < len(hits):
        notes.append("ignored copies on alt/fix/random contigs: "
                     + ", ".join(sorted({t.chrom for t in hits if '_' in t.chrom})))
        hits = primary

    if ref is not None:
        present = [t for t in hits if ref.resolve(t.chrom)]
        if not present:
            raise LookupError(f"contig {hits[0].chrom} of {query} is not in the FASTA")
        hits = present

    uniq = {}
    for t in hits:
        uniq.setdefault((t.chrom, t.strand, tuple(t.exons_g), t.cds_start, t.cds_end), t)
    hits = list(uniq.values())
    if len(hits) > 1:
        locs = ", ".join(f"{t.name}@{t.chrom}:{t.tx_start}-{t.tx_end}" for t in hits)
        raise LookupError(f"{query} maps to multiple loci ({locs}); pick one with --chrom")
    return hits[0], notes
