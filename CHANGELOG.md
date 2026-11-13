# Changelog

## 1.0.0 — 2026-11-13
First stable public release.
- `vcfart generate`: SNV / deletion / insertion / codon (DMS) saturation from HGVS c./n./g. ranges,
  exon/intron/cds shortcuts, `--exonic-only`, left-normalized de-duplicated VCF 4.2, TSV with HGVS c./p.
  + SO consequences, optional oligo-library FASTA with adapters.
- `vcfart info`: exon/intron table in c. and genomic coordinates + build sanity checks.
- `vcfart fetch-annotation`: UCSC refGene / ncbiRefSeq download into `~/.vcfart`.
- MultiQC-style coloured log with a consequence bar chart (`--no-color`, `NO_COLOR`).
- Reads UCSC genePred (incl. ANNOVAR `humandb/*_refGene.txt`) and GTF; no third-party dependencies.
- Validated on hg38: reproduces the original in-house MAVE scripts (1,332 SNVs) and matches ANNOVAR
  refGene c./p. for all 6,669 CDS SNVs of DRC1 (`examples/compare_annovar.py`).
- Korean README (`README.ko.md`) and a homepage (`docs/`, GitHub Pages).
