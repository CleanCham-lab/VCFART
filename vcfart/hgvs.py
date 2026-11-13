"""HGVS c./n./p. descriptions and Sequence-Ontology consequences."""

from .sequtils import aa3, revcomp, translate, translate_codon
from .variants import ref_of, shift

SEVERITY = [
    "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant",
    "stop_lost", "start_lost", "inframe_insertion", "inframe_deletion", "missense_variant",
    "protein_altering_variant", "splice_donor_5th_base_variant", "splice_region_variant",
    "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
    "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
    "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant",
    "intron_variant", "upstream_gene_variant", "downstream_gene_variant",
]
_RANK = {t: i for i, t in enumerate(SEVERITY)}


def _common_prefix(a, b):
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


class Annotator:
    def __init__(self, tx, ref):
        self.tx = tx
        self.ref = ref
        self.chrom = ref.resolve(tx.chrom)
        self._wt_prot = None

    # ------------------------------------------------------------------ c.
    def hgvs_edit(self, edit):
        """3'-shift along the transcript (HGVS 3' rule, in genomic sequence)."""
        return shift(edit, self.ref, 1 if self.tx.plus else -1)

    def _t(self, seq):
        return seq if self.tx.plus else revcomp(seq)

    def _rng(self, x, y):
        P = self.tx.pos
        if x == y:
            return P(x)
        return f"{P(x)}_{P(y)}" if self.tx.plus else f"{P(y)}_{P(x)}"

    def describe_c(self, e3):
        tx, pre = self.tx, self.tx.prefix
        r = ref_of(e3, self.ref)
        a = e3.alt
        if r and a:
            if len(r) == 1 and len(a) == 1:
                return f"{pre}{tx.pos(e3.start)}{self._t(r)}>{self._t(a)}"
            return f"{pre}{self._rng(e3.start, e3.end)}delins{self._t(a)}"
        if r:
            return f"{pre}{self._rng(e3.start, e3.end)}del"
        k, s = len(a), e3.start
        if tx.plus and self.ref.fetch(self.chrom, s - k, s - 1) == a:
            return f"{pre}{self._rng(s - k, s - 1)}dup"
        if not tx.plus and self.ref.fetch(self.chrom, s, s + k - 1) == a:
            return f"{pre}{self._rng(s, s + k - 1)}dup"
        flank = f"{tx.pos(s - 1)}_{tx.pos(s)}" if tx.plus else f"{tx.pos(s)}_{tx.pos(s - 1)}"
        return f"{pre}{flank}ins{self._t(a)}"

    # --------------------------------------------------------- consequences
    def _base_terms(self, g):
        tx = self.tx
        terms = set()
        kind, k = tx.region(g)
        if kind == "intron":
            _, off = tx.g_to_n(g)
            terms.add("intron_variant")
            if off in (1, 2):
                terms.add("splice_donor_variant")
            if off in (-1, -2):
                terms.add("splice_acceptor_variant")
            if off == 5:
                terms.add("splice_donor_5th_base_variant")
            if 3 <= off <= 6:
                terms.add("splice_donor_region_variant")
            if -17 <= off <= -3:
                terms.add("splice_polypyrimidine_tract_variant")
            if 3 <= abs(off) <= 8:
                terms.add("splice_region_variant")
        elif kind == "exon":
            n, _ = tx.g_to_n(g)
            n0, n1 = tx.exon_n[k - 1]
            if (k > 1 and n - n0 < 3) or (k < tx.n_exons and n1 - n < 3):
                terms.add("splice_region_variant")
            if tx.coding:
                if n < tx.cds_start_n:
                    terms.add("5_prime_UTR_variant")
                elif n > tx.cds_end_n:
                    terms.add("3_prime_UTR_variant")
                else:
                    terms.add("CDS")
            else:
                terms.add("non_coding_transcript_exon_variant")
        elif kind == "upstream":
            terms.add("upstream_gene_variant")
        elif kind == "downstream":
            terms.add("downstream_gene_variant")
        return terms

    def annotate(self, v):
        tx = self.tx
        e3 = self.hgvs_edit(v.edit)
        v.hgvs_c = self.describe_c(e3)
        if e3.ref_len > 0:
            bases = list(range(e3.start, e3.end + 1))
            terms = set().union(*(self._base_terms(g) for g in bases))
        else:
            bases = [e3.start - 1, e3.start]
            t1, t2 = self._base_terms(bases[0]), self._base_terms(bases[1])
            terms = t1 | t2
            for site in ("splice_donor_variant", "splice_acceptor_variant"):
                if not (site in t1 and site in t2):
                    terms.discard(site)
            if "CDS" in terms and not ("CDS" in t1 and "CDS" in t2):
                terms.discard("CDS")
        hgvs_p = ""
        if "CDS" in terms:
            terms.discard("CDS")
            hgvs_p, cterms = self.coding_effect(e3)
            terms |= cterms
        v.hgvs_p = hgvs_p
        v.consequence = "&".join(sorted(terms, key=lambda t: _RANK.get(t, 99)))
        labels = []
        ordered = bases if tx.plus else list(reversed(bases))
        for g in ordered:
            lab = tx.region_label(g)
            if lab not in labels:
                labels.append(lab)
        v.region = "-".join(labels)
        return v

    # ------------------------------------------------------------------ p.
    @property
    def wt_prot(self):
        if self._wt_prot is None:
            self._wt_prot = translate(self.tx.cds_seq(self.ref))
        return self._wt_prot

    def coding_effect(self, e3):
        tx = self.tx
        cs, ce = tx.cds_start_n, tx.cds_end_n
        if e3.ref_len > 0:
            (ns, os_), (ne, oe) = tx.g_to_n(e3.start), tx.g_to_n(e3.end)
            if os_ or oe or abs(ns - ne) != e3.end - e3.start:
                return "p.?", {"coding_sequence_variant"}
            lo, hi = min(ns, ne), max(ns, ne)
        else:
            (n1, o1), (n2, o2) = tx.g_to_n(e3.start - 1), tx.g_to_n(e3.start)
            if o1 or o2 or abs(n1 - n2) != 1:
                return "p.?", {"coding_sequence_variant"}
            lo = max(n1, n2)
            hi = lo - 1
        if lo < cs or hi > ce:
            terms = {"coding_sequence_variant"}
            if lo < cs:
                terms.add("start_lost")
                return "p.Met1?", terms
            terms.add("stop_lost")
            return "p.?", terms

        tseq = tx.transcript_seq(self.ref)
        s = tseq[cs - 1:]
        i0, i1 = lo - cs, hi - cs
        alt = self._t(e3.alt)
        mut = s[:i0] + alt + s[i1 + 1:]
        delta = len(alt) - (i1 - i0 + 1)
        W_all = self.wt_prot
        ci = i0 // 3

        if ci == 0 and W_all[:1] == "M" and translate_codon(mut[0:3]) != "M":
            return "p.Met1?", {"start_lost"}

        if delta == 0:
            return self._substitution(W_all, mut, ci, i1 // 3)

        W = W_all[ci:]
        M = translate(mut[ci * 3:], stop_at_stop=True)
        p = _common_prefix(W, M)
        pos = ci + p + 1
        if "X" in W[p:p + 1] + M[p:p + 1]:
            return "p.?", {"coding_sequence_variant"}

        if delta % 3:
            terms = {"frameshift_variant"}
            if p >= len(M) or p >= len(W):
                return "p.?", terms
            if W[p] == "*":
                terms.add("stop_lost")
                return self._ext(pos, M[p:]), terms
            if M[p] == "*":
                terms.add("stop_gained")
                return f"p.{aa3(W[p])}{pos}Ter", terms
            j = M.find("*", p)
            n = f"{j - p + 1}" if j >= 0 else "?"
            return f"p.{aa3(W[p])}{pos}{aa3(M[p])}fsTer{n}", terms

        terms = {"inframe_insertion" if delta > 0 else "inframe_deletion"}
        if p >= len(W):
            return "p.?", terms
        if p < len(M) and M[p] == "*" and W[p] != "*":
            terms.add("stop_gained")
            return f"p.{aa3(W[p])}{pos}Ter", terms

        expected = len(W) + delta // 3
        stop0 = tx.cds_len - 3                      # 0-based index of the stop codon in s
        if e3.ref_len > 0:
            stop_hit = i0 <= stop0 + 2 and i1 >= stop0
        else:
            stop_hit = stop0 < i0 <= stop0 + 2
        if stop_hit and not (M.endswith("*") and len(M) == expected):
            terms.add("stop_lost")
            k = len(W_all) - 1 - ci                 # WT Ter index within the tails
            return self._ext(len(W_all), M[k:] if k < len(M) else ""), terms
        if "*" not in M:
            return "p.?", terms

        if M.endswith("*") and len(M) < expected:
            terms.add("stop_gained")
            q = i1 // 3 - ci
            wd = W[p:q + 1] if q >= p else ""
            md = M[p:]
            return self._delins(W_all, ci + p, wd, md), terms

        smax = min(len(W), len(M)) - p
        t = 0
        while t < smax and W[len(W) - 1 - t] == M[len(M) - 1 - t]:
            t += 1
        wd, md = W[p:len(W) - t], M[p:len(M) - t]
        A = ci + p  # absolute 0-based index of first changed residue
        if not wd and not md:
            return "p.=", {"synonymous_variant"}
        if not md:
            if len(wd) == 1:
                return f"p.{aa3(wd)}{A + 1}del", terms
            return f"p.{aa3(wd[0])}{A + 1}_{aa3(wd[-1])}{A + len(wd)}del", terms
        if not wd:
            k = len(md)
            if A - k >= 0 and W_all[A - k:A] == md:
                if k == 1:
                    return f"p.{aa3(md)}{A}dup", terms
                return f"p.{aa3(md[0])}{A - k + 1}_{aa3(md[-1])}{A}dup", terms
            if A == 0:
                return "p.?", terms
            return f"p.{aa3(W_all[A - 1])}{A}_{aa3(W_all[A])}{A + 1}ins{aa3(md)}", terms
        return self._delins(W_all, A, wd, md), terms

    @staticmethod
    def _ext(pos, M_from):
        if not M_from:
            return f"p.Ter{pos}?"
        j = M_from.find("*")
        n = f"{j + 1}" if j >= 0 else "?"
        return f"p.Ter{pos}{aa3(M_from[0])}extTer{n}"

    @staticmethod
    def _delins(W_all, A, wd, md):
        if not wd:
            return f"p.{aa3(W_all[A - 1])}{A}_{aa3(W_all[A])}{A + 1}ins{aa3(md)}"
        if len(wd) == 1:
            first = f"{aa3(wd)}{A + 1}"
        else:
            first = f"{aa3(wd[0])}{A + 1}_{aa3(wd[-1])}{A + len(wd)}"
        return f"p.{first}delins{aa3(md)}"

    def _substitution(self, W_all, mut, ci, cj):
        w = W_all[ci:cj + 1]
        m = translate(mut[ci * 3:(cj + 1) * 3])
        if "X" in w or "X" in m:
            return "p.?", {"coding_sequence_variant"}
        if w == m:
            if "*" in w and len(w) == 1:
                return f"p.Ter{ci + 1}=", {"stop_retained_variant"}
            if len(w) == 1:
                return f"p.{aa3(w)}{ci + 1}=", {"synonymous_variant"}
            return "p.=", {"synonymous_variant"}
        p = _common_prefix(w, m)
        t = 0
        while t < len(w) - p and w[len(w) - 1 - t] == m[len(m) - 1 - t]:
            t += 1
        wd, md = w[p:len(w) - t], m[p:len(m) - t]
        A = ci + p
        if "*" in wd and wd[0] == "*":
            M = translate(mut[A * 3:], stop_at_stop=True)
            return self._ext(A + 1, M), {"stop_lost"}
        if md[0] == "*":
            return f"p.{aa3(wd[0])}{A + 1}Ter", {"stop_gained"}
        if len(wd) == 1:
            return f"p.{aa3(wd)}{A + 1}{aa3(md)}", {"missense_variant"}
        terms = {"stop_gained"} if "*" in md else {"missense_variant"}
        return self._delins(W_all, A, wd, md), terms
