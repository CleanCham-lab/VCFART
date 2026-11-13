#!/usr/bin/env bash
# VCFART demo on a tiny synthetic genome (minus-strand gene DEMO1, transcript NM_DEMO.1)
set -euo pipefail
cd "$(dirname "$0")"

# 1) look at the transcript structure (c. and genomic coordinates)
vcfart info -f demo.fa -a demo_refGene.txt -t NM_DEMO.1

# 2) exon 3 with 50 bp of flanking intron on each side: all SNVs + 1-bp del/ins
vcfart generate -f demo.fa -a demo_refGene.txt -t NM_DEMO.1 \
    -r ex3=exon3:50 -m snv,del,ins -o out/demo_ex3 --oligo

# 3) the same region written as explicit HGVS ranges (intron / exon / intron)
vcfart generate -q -f demo.fa -a demo_refGene.txt -t NM_DEMO.1 \
    -r int2=c.181-50_c.181-1 -r ex3=c.181_c.342 -r int3=c.342+1_c.342+50 -o out/demo_split

# 4) deep-mutational-scanning style codon swaps for codons 61-70
vcfart generate -q -f demo.fa -a demo_refGene.txt -t NM_DEMO.1 \
    -r c.181_c.210 -m codon --codon-set aa --codon-pick min-edit -o out/demo_dms
