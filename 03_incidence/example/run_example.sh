#!/bin/bash
# Variant-level and gene-level incidence for all ClinVar P/LP variants of
# ATAD3A (AR), ATP6AP1 and BCAP31 (XLR).
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

python ../01_annotate_variant_incidence.py -i example.outfilt.vcf -o example.incidence.vcf

python ../02_variant_table.py -i example.incidence.vcf \
    -g example_genes.tsv --aso example_aso_organs.tsv -o example.variant_table.tsv

python ../03_gene_incidence_table.py -i example.incidence.vcf \
    -g example_genes.tsv --aso-organs example_aso_organs.tsv \
    -o example.gene_incidence.tsv \
    --all-transcripts-output example.gene_incidence.alltrx.tsv \
    --aso-output example.gene_incidence.ASO_target.tsv

for f in example.incidence.vcf example.variant_table.tsv example.gene_incidence.tsv \
         example.gene_incidence.alltrx.tsv example.gene_incidence.ASO_target.tsv; do
    if diff -q "$f" "expected/$f" > /dev/null; then
        echo "$f matches expected/$f"
    else
        echo "$f differs from expected/$f"
    fi
done
