"""Output writers: VCF 4.2, TSV table, oligo-library FASTA."""

import datetime

from . import __version__

INFO_DEFS = [
    ("VA_TX", "1", "String", "Transcript used by VCFART"),
    ("VA_GENE", "1", "String", "Gene symbol"),
    ("VA_HGVSc", "1", "String", "HGVS c./n. description (3'-shifted on the transcript)"),
    ("VA_HGVSp", "1", "String", "HGVS p. description (predicted; '=' encoded as %3D)"),
    ("VA_CSQ", "1", "String", "Sequence Ontology consequence(s), '&'-separated, most severe first"),
    ("VA_REGION", "1", "String", "Exon/intron location on the transcript (e.g. intron8, exon9)"),
    ("VA_TYPE", "1", "String", "Variant class: SNV, MNV, DEL, INS, DELINS"),
    ("VA_ORIGIN", ".", "String", "Generator(s) that produced this variant: snv, delN, insN, codon"),
    ("VA_SET", ".", "String", "Target region label(s) the variant was generated from"),
]


def _info_escape(v):
    return (str(v).replace("%", "%25").replace(";", "%3B").replace("=", "%3D")
            .replace(",", "%2C").replace(" ", "%20"))


def write_vcf(path, variants, tx, ref, build, command, sites_only=False, gt="0/1",
              sample="VCFART", id_style="hgvs", all_contigs=False):
    chrom = ref.resolve(tx.chrom)
    contigs = ref.order if all_contigs else [chrom]
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##fileDate={datetime.date.today():%Y%m%d}\n")
        out.write(f"##source=VCFART v{__version__} (artificial saturation-mutagenesis VCF)\n")
        out.write(f"##reference=file://{ref.path}\n")
        asm = f",assembly={build}" if build else ""
        for c in contigs:
            out.write(f"##contig=<ID={c},length={ref.length(c)}{asm}>\n")
        out.write(f"##VCFART_transcript={tx.name} {tx.gene} {tx.chrom}:{tx.tx_start}-{tx.tx_end} "
                  f"strand={tx.strand} annotation={tx.source}\n")
        out.write(f"##VCFART_command={command}\n")
        for i, n, t, d in INFO_DEFS:
            out.write(f'##INFO=<ID={i},Number={n},Type={t},Description="{d}">\n')
        cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        if not sites_only:
            out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            cols += ["FORMAT", sample]
        out.write("\t".join(cols) + "\n")
        for i, v in enumerate(variants, 1):
            if id_style == "hgvs":
                vid = v.hgvs_c
            elif id_style == "index":
                vid = f"VA{i:06d}"
            else:
                vid = "."
            info = [
                f"VA_TX={tx.name}", f"VA_GENE={_info_escape(tx.gene)}",
                f"VA_HGVSc={_info_escape(v.hgvs_c)}",
            ]
            if v.hgvs_p:
                info.append(f"VA_HGVSp={_info_escape(v.hgvs_p)}")
            info += [
                f"VA_CSQ={v.consequence or 'intergenic_variant'}",
                f"VA_REGION={v.region}", f"VA_TYPE={v.mut_type}",
                "VA_ORIGIN=" + ",".join(_info_escape(o) for o in v.origins),
                "VA_SET=" + ",".join(_info_escape(s) for s in v.sets),
            ]
            row = [v.chrom, str(v.pos), vid, v.ref, v.alt, ".", ".", ";".join(info)]
            if not sites_only:
                row += ["GT", gt]
            out.write("\t".join(row) + "\n")


TSV_COLS = ["set", "chrom", "pos", "ref", "alt", "type", "origin", "transcript", "gene",
            "hgvs_c", "hgvs_p", "consequence", "region"]


def write_tsv(path, rows, tx, oligo=False):
    cols = TSV_COLS + (["oligo_id", "oligo_len", "oligo_seq"] if oligo else [])
    with open(path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for set_label, v, olig in rows:
            r = [set_label, v.chrom, str(v.pos), v.ref, v.alt, v.mut_type, ",".join(v.origins),
                 tx.name, tx.gene, f"{tx.name}:{v.hgvs_c}", v.hgvs_p or ".",
                 v.consequence or ".", v.region]
            if oligo:
                oid, seq = olig if olig else (".", "")
                r += [oid, str(len(seq)) if seq else ".", seq or "."]
            out.write("\t".join(r) + "\n")


def write_fasta(path, records, width=0):
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            if width:
                for i in range(0, len(seq), width):
                    out.write(seq[i:i + width] + "\n")
            else:
                out.write(seq + "\n")
