"""Target-region parsing.

Accepted forms for ``-r/--region`` (optionally prefixed with ``LABEL=``):

  c.5075-100_c.5075-1     HGVS c. range (second 'c.' optional: c.5075-100_5075-1)
  c.-30_c.*20             UTRs / beyond the transcript are fine
  n.10_n.200              non-coding transcripts
  c.5075                  single position
  g.43063774_43064051     genomic (transcript chromosome)
  chr17:43063774-43064051 genomic with chromosome
  exon17                  whole exon 17 (transcript order)
  exon17:100              exon 17 + 100 bp intron on both sides
  exon17:100:50           exon 17 + 100 bp 5' flank + 50 bp 3' flank (transcript orientation)
  intron16                whole intron 16
  cds                     c.1 .. last base of the stop codon (combine with --exonic-only)
"""

import re
from dataclasses import dataclass


@dataclass
class Region:
    label: str
    chrom: str
    start: int
    end: int
    spec: str
    hgvs: str = ""

    @property
    def length(self):
        return self.end - self.start + 1


_EXON_RE = re.compile(r"^(?:exon|ex|e)(\d+)(?::(\d+)(?::(\d+))?)?$", re.I)
_INTRON_RE = re.compile(r"^(?:intron|int|ivs|i)(\d+)$", re.I)
_GEN_RE = re.compile(r"^(?:(?P<chrom>[^:\s]+):)?(?:g\.)?(?P<a>[\d,]+)(?:[-_](?:g\.)?(?P<b>[\d,]+))?$")


def parse_region(spec, tx, ref, index):
    label = None
    body = spec.strip()
    if "=" in body:
        label, body = body.split("=", 1)
        label, body = label.strip(), body.strip()
    label = label or f"R{index}"

    if body.lower() == "cds":
        if not tx.coding:
            raise ValueError(f"{tx.name} is non-coding (no CDS)")
        a, b = tx.n_to_g(tx.cds_start_n), tx.n_to_g(tx.cds_end_n)
        return _finish(label, tx, ref, min(a, b), max(a, b), spec)

    m = _EXON_RE.match(body)
    if m:
        k = int(m.group(1))
        if not 1 <= k <= tx.n_exons:
            raise ValueError(f"{tx.name} has {tx.n_exons} exons (asked exon {k})")
        f5 = int(m.group(2) or 0)
        f3 = int(m.group(3)) if m.group(3) is not None else f5
        s, e = tx.exons[k - 1]
        if tx.plus:
            start, end = s - f5, e + f3
        else:
            start, end = s - f3, e + f5
        return _finish(label, tx, ref, start, end, spec)

    m = _INTRON_RE.match(body)
    if m:
        k = int(m.group(1))
        if not 1 <= k < tx.n_exons:
            raise ValueError(f"{tx.name} has {tx.n_exons - 1} introns (asked intron {k})")
        a = tx.exons[k - 1]
        b = tx.exons[k]
        start = min(a[1], b[1]) + 1
        end = max(a[0], b[0]) - 1
        return _finish(label, tx, ref, start, end, spec)

    low = body.lower()
    if low.startswith(("c.", "n.")):
        kind = body[:2]
        parts = body[2:].split("_")
        if len(parts) > 2:
            raise ValueError(f"bad range '{body}'")
        gs = [tx.parse_pos(p if p.lower().startswith(("c.", "n.")) else kind + p) for p in parts]
        return _finish(label, tx, ref, min(gs), max(gs), spec)

    m = _GEN_RE.match(body)
    if m and (m.group("chrom") or low.startswith("g.")):
        chrom = m.group("chrom")
        if chrom and ref.resolve(chrom) != ref.resolve(tx.chrom):
            raise ValueError(f"region {body} is on {chrom}, transcript is on {tx.chrom}")
        a = int(m.group("a").replace(",", ""))
        b = int((m.group("b") or m.group("a")).replace(",", ""))
        return _finish(label, tx, ref, min(a, b), max(a, b), spec)

    raise ValueError(
        f"cannot parse region '{spec}'. Use c.X_Y, g.X_Y, chr:X-Y, exonN[:flank], intronN")


def _finish(label, tx, ref, start, end, spec):
    chrom = ref.resolve(tx.chrom)
    if start < 1 or end > ref.length(chrom):
        raise ValueError(f"region {spec} ({start}-{end}) is outside {chrom}")
    a, b = (start, end) if tx.plus else (end, start)
    hg = f"{tx.prefix}{tx.pos(a)}_{tx.pos(b)}" if start != end else f"{tx.prefix}{tx.pos(a)}"
    return Region(label, chrom, start, end, spec, hg)
