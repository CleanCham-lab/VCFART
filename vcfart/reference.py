"""Indexed FASTA access without external dependencies (samtools-compatible .fai)."""

import os
import sys

from .utils import log, warn


class FastaError(Exception):
    pass


class FastaReference:
    """Random access to a (plain-text) FASTA file through a .fai index.

    * uses an existing ``<fasta>.fai`` (samtools faidx) if present
    * otherwise builds the index once (and tries to save it next to the FASTA)
    * bgzip-compressed FASTA (.gz) is supported only if ``pysam`` is installed
    """

    CHUNK = 1 << 16

    def __init__(self, path):
        self.path = os.path.abspath(os.path.expanduser(path))
        if not os.path.exists(self.path):
            raise FastaError(f"FASTA not found: {self.path}")
        self._pysam = None
        self._cache = {}
        if self.path.endswith((".gz", ".bgz")):
            try:
                import pysam
            except ImportError:
                raise FastaError(
                    "Compressed FASTA needs pysam (pip install pysam) -- "
                    "or use an uncompressed FASTA (gunzip hg38.fa.gz)."
                )
            self._pysam = pysam.FastaFile(self.path)
            self.index = {
                n: (l, None, None, None)
                for n, l in zip(self._pysam.references, self._pysam.lengths)
            }
            self.order = list(self._pysam.references)
            return
        self.index, self.order = self._load_index()
        self._fh = open(self.path, "rb")

    # ------------------------------------------------------------------ index
    def _load_index(self):
        fai = self.path + ".fai"
        if os.path.exists(fai):
            index, order = {}, []
            with open(fai) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 5:
                        continue
                    index[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
                    order.append(f[0])
            return index, order
        log(f"🗂️  no .fai index for {self.path}; building one (one-time, ~1 min for hg38)", "reference")
        index, order = self._build_index()
        try:
            with open(fai, "w") as out:
                for name in order:
                    l, off, lb, lw = index[name]
                    out.write(f"{name}\t{l}\t{off}\t{lb}\t{lw}\n")
            log(f"🗂️  saved {fai}", "reference")
        except OSError:
            warn(f"Could not write {fai} (read-only?). Index kept in memory only; "
                 "consider `samtools faidx` once.")
        return index, order

    def _build_index(self):
        index, order = {}, []
        name = None
        length = offset = lb = lw = 0
        short_line_seen = False
        pos = 0
        with open(self.path, "rb") as fh:
            for raw in fh:
                n = len(raw)
                if raw.startswith(b">"):
                    if name is not None:
                        index[name] = (length, offset, lb, lw)
                    name = raw[1:].split()[0].decode()
                    order.append(name)
                    length = lb = lw = 0
                    offset = pos + n
                    short_line_seen = False
                else:
                    seq = raw.rstrip(b"\r\n")
                    if seq:
                        if lb == 0:
                            lb, lw = len(seq), n
                        elif short_line_seen:
                            raise FastaError(
                                f"Irregular line lengths in sequence {name}; "
                                "please reformat or use samtools faidx."
                            )
                        elif len(seq) != lb:
                            short_line_seen = True
                        length += len(seq)
                pos += n
            if name is not None:
                index[name] = (length, offset, lb, lw)
        return index, order

    # --------------------------------------------------------------- queries
    def resolve(self, chrom):
        """Return the FASTA's own name for ``chrom`` (chr-prefix / MT aliasing)."""
        if chrom in self.index:
            return chrom
        cands = []
        if chrom.startswith("chr"):
            core = chrom[3:]
            cands.append(core)
            if core == "M":
                cands.append("MT")
        else:
            cands.append("chr" + chrom)
            if chrom == "MT":
                cands.append("chrM")
        for c in cands:
            if c in self.index:
                return c
        return None

    def length(self, chrom):
        return self.index[chrom][0]

    def _read(self, chrom, start0, end0):
        length, off, lb, lw = self.index[chrom]
        start0 = max(0, start0)
        end0 = min(length, end0)
        if start0 >= end0:
            return ""
        if self._pysam is not None:
            return self._pysam.fetch(chrom, start0, end0).upper()
        b0 = off + (start0 // lb) * lw + start0 % lb
        b1 = off + ((end0 - 1) // lb) * lw + (end0 - 1) % lb + 1
        self._fh.seek(b0)
        data = self._fh.read(b1 - b0)
        return data.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()

    def fetch(self, chrom, start, end):
        """1-based inclusive [start, end]; returns '' if empty. Clamped to contig."""
        if end < start:
            return ""
        c = self.CHUNK
        s0, e0 = start - 1, end
        out = []
        for k in range(max(0, s0) // c, (e0 - 1) // c + 1):
            key = (chrom, k)
            chunk = self._cache.get(key)
            if chunk is None:
                if len(self._cache) > 512:
                    self._cache.clear()
                chunk = self._read(chrom, k * c, (k + 1) * c)
                self._cache[key] = chunk
            lo = max(s0, k * c) - k * c
            hi = min(e0, (k + 1) * c) - k * c
            out.append(chunk[lo:hi])
        return "".join(out)

    def base(self, chrom, pos):
        s = self.fetch(chrom, pos, pos)
        return s if s else "N"


def open_reference(path):
    try:
        return FastaReference(path)
    except FastaError as e:
        sys.exit(f"[vcfart] ERROR: {e}")
