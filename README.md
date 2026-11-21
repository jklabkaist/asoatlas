# ASOatlas — ASO amenability and incidence of ClinVar pathogenic variants

Which pathogenic variants could be corrected by an antisense oligonucleotide (ASO), and
how many patients they account for.

The pipeline keeps the ClinVar variants with a pathogenic or likely pathogenic germline
submission, predicts for every variant-transcript pair whether an ASO could correct its
effect, and turns All of Us and UK Biobank allele frequencies into disease incidence,
expected patients per year and the share of a gene's patients an ASO would reach — each
with a 95% Bayesian credible interval.

## Pipeline

| | step | output |
| --- | --- | --- |
| 1 | [`01_clinvar_pathogenic_variant_filteration`](01_clinvar_pathogenic_variant_filteration) | VCF of P/LP variants, ambiguous ones flagged |
| 2 | [`02_ASO_amenability`](02_ASO_amenability) | `Probably` / `Possibly` / `Unlikely amenable` per variant-transcript pair |
| 3 | [`03_incidence`](03_incidence) | incidence, patients per year and %amenability, per variant and per gene |

[`aso_amenability`](aso_amenability) is step 2 as an importable module, usable on its own.
Each step folder has a README with its fields, columns and options.

## Install

```bash
pip install -r aso_amenability/requirements.txt          # python >= 3.9
./aso_amenability/setup_data.sh <reference file directory>
```

Step 1 additionally needs Ensembl VEP with a GRCh38 offline cache and the REVEL, NMD and
EVE plugins. A GPU is optional — SpliceAI also runs on CPU, more slowly.

Reference files are not distributed here; `setup_data.sh` lists each one and where it
comes from (hg38.2bit, the SpliceAI models, MANE v1.0, the GENCODE v41 MANE refFlat and
the LaBranchoR branchpoints). All of Us and UK Biobank allele counts have to be produced
inside their own controlled-access environments and joined onto the VCF, as described at
the bottom of [`01_variant_summary_to_vcf.py`](01_clinvar_pathogenic_variant_filteration/01_variant_summary_to_vcf.py).

## Usage

```bash
IN=variant_summary_2025-07.GRCh38       # ClinVar release
XML=ClinVarRCVRelease_2025-07.xml
GENES=union_AR_XLR.tsv                  # TSV with an HGNC_num column
ORGANS=union_AR_XLR_gt_gene.tsv         # TSV with HGNC_num and GT_Target_Organs

# 1. ClinVar P/LP variant set
cd 01_clinvar_pathogenic_variant_filteration
python 01_variant_summary_to_vcf.py -i variant_summary_2025-07.txt -o $IN.sorted.vcf
#   ... join the gnomAD / All of Us / UK Biobank allele counts here ...
sbatch 02_run_vep.sh $IN.sorted.joined.vcf                      # calls 03_curate_vep.py
python 04_rcv_xml_to_tsv.py -i $XML -o clinvar.scv.tsv
python 05_annotate_low_penetrance.py -i clinvar.scv.tsv -o clinvar.scv.lowP.tsv
python 06_filter_pathogenic_variants.py -i $IN.sorted.joined.VEP.MANE.vcf \
                                        -s clinvar.scv.lowP.tsv -o $IN.cs_plp.vcf
python 07_filter_af_outlier.py -i $IN.cs_plp.vcf -o $IN.cs_plp.outfilt.vcf
python 08_merge_duplicate_alleles.py -i $IN.cs_plp.outfilt.vcf -o $IN.cs_plp.outfilt.dedup.vcf

# 2. ASO amenability
cd ../02_ASO_amenability
python annotate_vcf.py -i $IN.cs_plp.outfilt.dedup.vcf -o $IN.amenability.vcf -b 500
#   on a GPU node:  sbatch slurm_annotate_vcf.sh $IN.cs_plp.outfilt.dedup.vcf

# 3. Incidence
cd ../03_incidence
python 01_annotate_variant_incidence.py -i $IN.amenability.vcf -o $IN.incidence.vcf
python 02_variant_table.py -i $IN.incidence.vcf -g $GENES --aso $ORGANS -o variant_table.tsv
python 03_gene_incidence_table.py -i $IN.incidence.vcf -g $GENES --aso-organs $ORGANS \
       -o gene_incidence.tsv \
       --all-transcripts-output gene_incidence.alltrx.tsv \
       --aso-output gene_incidence.ASO_target.tsv
```

## Examples

Both ship with their inputs and expected outputs, and compare the two:

```bash
cd 02_ASO_amenability/example && ./run_example.sh   # 6 variants -> amenability VCF
cd 03_incidence/example       && ./run_example.sh   # 63 variants of 3 genes -> incidence
```

## Amenability as a module

```python
from aso_amenability import annotate_variants

vepD = {...}          # one VEP CSQ entry as a {field: value} dict
resL = annotate_variants([('chr1', '2406956', 'C', 'T', 'ENST00000447513.7', vepD)])

resL[0]['ASO_amenability']    # 'Probably amenable' | 'Possibly amenable' | 'Unlikely amenable'
resL[0]['STATUS_mis_splice']  # the predicted mis-splicing behind the call
```

`aso_amenability/vendor/` holds the two in-house packages it needs — `jkspliceai`, the
modified SpliceAI API, in full, and the part of `jklib` it uses — importable without
being installed.

## Notes

- Every quantity is computed per variant-transcript pair, so a variant annotated to two
  overlapping genes is scored against each gene's own allele pool.
- ClinVar records that describe the same allele under different VariationIDs are merged
  in step 1, so each allele enters the incidence once.
- The gene set is matched to the VEP annotations by HGNC ID and Ensembl transcript ID, not
  by gene symbol, since ClinVar, MANE v1.0 and VEP can name the same gene differently.
- A variant is `Probably amenable` when it gains mis-splicing without severe coding or
  canonical-splice damage, `Possibly amenable` when it causes exon skipping or retention.
- Incidence uses `(sum AF)^2` for autosomal recessive genes and `0.5 * sum(female AF)`
  for X-linked ones, excluding combinations of two hypomorphic alleles, and is scaled to
  3.7M US and 690K UK annual births.
- SpliceAI scores can differ in the third decimal between GPU and CPU, and the credible
  intervals in their last digit between SciPy versions; neither changes a call.

## Code authorship and AI assistance

The analysis, the pipeline and the underlying code were designed and written by the
authors. Claude Opus 5 (Anthropic), used through Claude Code between October and
November 2026, assisted with preparing this release: reorganising the working scripts
into the three steps above, replacing hard-coded paths with command-line arguments,
writing the step READMEs and the runnable examples, and reviewing the reorganised code
against the original scripts for behavioural equivalence. It was not used to design the
method or to generate the results. Every figure and table in the paper comes from the
code as deposited here, and the authors are fully responsible for it.
