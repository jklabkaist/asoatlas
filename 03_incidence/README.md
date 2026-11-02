# 03. Variant-level and gene-level incidence

Turns the cohort allele frequencies of the P/LP variants into disease incidence: per
variant-transcript pair first, then per gene for the whole pathogenic set and for each
ASO-amenable subset. Every estimate carries a 95% Jeffreys credible interval.

Input: the VCF from step 02 (ClinVar P/LP + VEP + ASO amenability), carrying the
All of Us and UK Biobank allele counts joined in step 01.

```bash
# 1. incidence of every variant-transcript pair, written into the CSQ field
python 01_annotate_variant_incidence.py -i <step02_output>.vcf -o <step02_output>.bytrx.vcf

# 2. variant table: one row per variant-transcript pair of the gene set
python 02_variant_table.py \
    -i <step02_output>.bytrx.vcf \
    -g <AR_XLR_genes>.tsv \
    --aso <gene_target_organs>.tsv \
    -o variant_table.tsv

# 3. gene table: all pairs, one transcript per gene, and the ASO-target subset
python 03_gene_incidence_table.py \
    -i <step02_output>.bytrx.vcf \
    -g <AR_XLR_genes>.tsv \
    --aso-organs <gene_target_organs>.tsv \
    -o gene_incidence.tsv \
    --all-transcripts-output gene_incidence.alltrx.tsv \
    --aso-output gene_incidence.ASO_target.tsv
```

`-g` is a TSV with an `HGNC_num` column: the gene universe to tabulate (the AR / XLR
genes of the paper). Genes with no observed variant stay in the gene table with zero
incidence. The target-organ TSV (`HGNC_num`, `GT_Target_Organs`) fills the
`ASO_target_organs` column and defines the ASO-target subset, which also drops the genes
that already have an approved therapy. `--mane` and `--refflat` default to the annotation
files bundled with [`aso_amenability`](../aso_amenability). The incidence model and the
credible intervals live in [`incidence.py`](incidence.py), shared by steps 1 and 3.

## Unit of analysis: the variant-transcript pair

A variant can be annotated to several transcripts, and those transcripts may belong to
different, overlapping genes. Since the incidence of a variant depends on its **gene's**
pathogenic allele pool, the same variant has a different incidence in each gene, and
every quantity is keyed on `SYMBOL:Feature`. A variant annotated to two genes of the
gene set is counted in both, each time against that gene's own pool.

A VCF INFO field holds one value per variant and cannot represent this, so step 1 writes
the per-pair values **inside the CSQ field** (18 fields appended, AoU and UKBB). Steps 2
and 3 read them from the row's own transcript. The scalar INFO fields are kept only for
backward compatibility — they carry the primary pair (MANE Select first, then
alphabetical) plus a packed `<cohort>_NoBiHypo_bygene` field — and must not be used for
gene-aware analysis.

gnomAD is used only to flag hypomorphic alleles; its aggregation over heterogeneous
cohorts does not fit this incidence model, so no gnomAD incidence is reported.

## How the incidence is computed

| | AR genes | XLR genes (chrX) |
| --- | --- | --- |
| incidence of a variant set | `(sum AF)^2` minus the combinations without that set | `0.5 * sum(female AF)` of that set |
| hypomorphic alleles | bi-hypomorphic combinations are subtracted (>= 2 homozygotes in any cohort) | hypomorphic variants are dropped (>= 2 hemizygotes / female homozygotes) |

`PR` is the incidence, `EN` the expected number of patients among the newborns of one
year (3.7M for All of Us / US, 690K for UK Biobank / UK), `PRd` the fraction of the
gene's cases the variant set explains. A variant is used for a cohort only when it is
called in >= 30% of that cohort and is not flagged `Ambiguous` by step 01.

The interval propagates the Jeffreys posterior of the focal allele count,
`q ~ Beta(AC + 1/2, AN - AC + 1/2)`, through the same incidence functions; the gene
denominator of `PRd` is recomputed at the same draw. An allele pool with no observed
allele is left-censored, so its lower limit is 0 (one-sided interval).

## Step 1 fields

Per cohort, appended to every CSQ entry and mirrored into the scalar INFO fields:

| field | meaning |
| --- | --- |
| `PR_NoBiHypo` | incidence of the cases carrying this variant, in this gene |
| `EN_NoBiHypo` | expected patients per year |
| `PRd_NoBiHypo` | fraction of the gene's cases carrying this variant |
| `*_CIlo`, `*_CIhi` | 95% Jeffreys credible interval of each of the three |

INFO also keeps `AN_percent` (call rate) and `nhomalt_ge_2` (hypomorphic), both
gene-independent.

## Step 3 columns

One row per gene (or per gene-transcript pair with `--all-transcripts-output`), with
`<cohort>_` prefixed columns and each value followed by its `_CIlo` / `_CIhi` pair:

| variant set | columns |
| --- | --- |
| all pathogenic | `patho_PR_NoBiHypo`, `patho_EN_NoBiHypo`, allele / variant counts |
| all amenable (Probably + Possibly) | `amen_PR_NoBiHypo`, `amen_EN_NoBiHypo`, `amen_PRd_NoBiHypo` |
| Probably amenable | `prob_*` |
| recurrent amenable (ASO with EN >= 0.1 patients) | `amen_*_recurr`, `prob_*_recurr`, with the number of variants and of ASOs |
| the N highest-EN ASO groups | `*_top{N}` (`--top-n`, default 10) |

An "ASO" is one distinct (amenability call, splicing molecular consequence) group, i.e.
the variants a single ASO would address; where the molecular consequence is empty the
ClinVar variant name is used, making the variant its own group.

## Caveats

- `ASO_amenability` is a per-transcript call, so a variant can be amenable in one gene
  and not in another (FHL1 and GRIA3, whose MANE Select and MANE Plus Clinical
  transcripts disagree). Each pair is scored on its own call.
- Nine variants that experimental evidence showed to be unamenable are overwritten to
  `Unlikely amenable` in steps 2 and 3.

## Example

```bash
cd example && ./run_example.sh
```

All P/LP variants of ATAD3A (AR), ATP6AP1 and BCAP31 (XLR); the outputs are compared
with `example/expected/`. Note that the call-rate filter uses the largest AN found in
the input file, so a subset VCF like this one gives slightly different numbers than the
genome-wide run.
