# ASOatlas — ASO amenability and incidence of ClinVar pathogenic variants

Code for estimating which pathogenic variants are amenable to antisense oligonucleotide
(ASO) treatment, and how many patients they account for.

The pipeline starts from a ClinVar release, keeps the variants with a pathogenic or
likely pathogenic germline submission, predicts for every variant-transcript pair
whether an ASO could correct its effect, and converts population allele frequencies
(All of Us, UK Biobank) into disease incidence, expected patients per year and the
fraction of a gene's patients that a variant or an ASO would reach — each with a 95%
Bayesian credible interval.

```
01_clinvar_pathogenic_variant_filteration/   ClinVar P/LP variant set
02_ASO_amenability/                          ASO amenability annotation   (runner)
aso_amenability/                             ASO amenability             (module)
03_incidence/                                variant- and gene-level incidence (runner)
incidence/                                   incidence + credible intervals (module)
```

The numbered folders are the pipeline steps, run in order. The two unnumbered folders
are importable modules: `aso_amenability` is the amenability engine, `incidence` the
incidence engine; both are usable on their own.

## Requirements

- Python >= 3.9, with the packages in [`aso_amenability/requirements.txt`](aso_amenability/requirements.txt)
  (TensorFlow, the SpliceAI package, numpy / pandas / scipy, pytabix, twobitreader, PyYAML)
- Ensembl VEP with an offline GRCh38 cache and the REVEL, NMD and EVE plugins (step 01)
- A GPU is optional; SpliceAI runs on CPU, more slowly

Reference files are not distributed with this repository. `aso_amenability/setup_data.sh`
links them into place and lists where each one comes from:

| file | source |
| --- | --- |
| `hg38.2bit` | UCSC goldenPath |
| SpliceAI models `spliceai{1..5}.h5` | the SpliceAI package (Illumina) |
| `MANE.GRCh38.v1.0.summary.txt` | NCBI MANE release 1.0 |
| GENCODE v41 MANE refFlat | derived from the GENCODE v41 GTF |
| LaBranchoR branchpoints (GENCODE v40) | LaBranchoR predictions |

```bash
pip install -r aso_amenability/requirements.txt
./aso_amenability/setup_data.sh <directory holding the reference files>
```

Allele counts from All of Us and the UK Biobank must be generated inside their own
controlled-access environments and joined onto the VCF before step 03; the join step is
described in [`01_clinvar_pathogenic_variant_filteration/01_variant_summary_to_vcf.py`](01_clinvar_pathogenic_variant_filteration/01_variant_summary_to_vcf.py).

## Running the pipeline

```bash
IN=variant_summary_2025-07.GRCh38
CLINVAR_XML=ClinVarRCVRelease_2025-07.xml
GENES=union_AR_XLR.tsv                 # TSV with an HGNC_num column
ORGANS=union_AR_XLR_gt_gene.tsv        # TSV with HGNC_num and GT_Target_Organs

# ── 01. ClinVar P/LP variant set ───────────────────────────────────────────────
cd 01_clinvar_pathogenic_variant_filteration

python 01_variant_summary_to_vcf.py -i variant_summary_2025-07.txt -o $IN.sorted.vcf
#   ... left-join the gnomAD / UK Biobank / All of Us allele counts here ...
sbatch 02_run_vep.sh $IN.sorted.joined.vcf                     # calls 03_curate_vep.py
python 04_rcv_xml_to_tsv.py        -i $CLINVAR_XML -o clinvar.scv.tsv
python 05_annotate_low_penetrance.py -i clinvar.scv.tsv -o clinvar.scv.lowP.tsv
python 06_filter_pathogenic_variants.py -i $IN.sorted.joined.VEP.MANE.vcf \
                                        -s clinvar.scv.lowP.tsv -o $IN.cs_plp.vcf
python 07_filter_af_outlier.py -i $IN.cs_plp.vcf -o $IN.cs_plp.outfilt.vcf

# ── 02. ASO amenability ────────────────────────────────────────────────────────
cd ../02_ASO_amenability

python annotate_vcf.py -i $IN.cs_plp.outfilt.vcf -o $IN.amenability.vcf -b 500
# or on a GPU node:  sbatch slurm_annotate_vcf.sh $IN.cs_plp.outfilt.vcf

# ── 03. Incidence ──────────────────────────────────────────────────────────────
cd ../03_incidence

python 01_annotate_variant_incidence.py -i $IN.amenability.vcf -o $IN.incidence.vcf
python 02_variant_table.py   -i $IN.incidence.vcf -g $GENES --aso $ORGANS -o variant_table.tsv
python 03_gene_incidence_table.py -i $IN.incidence.vcf -g $GENES --aso-organs $ORGANS \
       -o gene_incidence.tsv \
       --all-transcripts-output gene_incidence.alltrx.tsv \
       --aso-output gene_incidence.ASO_target.tsv
```

Each step folder has a README with its own fields, columns and options:
[01](01_clinvar_pathogenic_variant_filteration/README.md),
[02](02_ASO_amenability/README.md),
[03](03_incidence/README.md).

## Using the modules directly

ASO amenability of a variant, given its VEP annotation:

```python
from aso_amenability import annotate_variants

vepD = {...}          # one VEP CSQ entry as a {field: value} dict
resL = annotate_variants([('chr1', '2406956', 'C', 'T', 'ENST00000447513.7', vepD)])

resL[0]['ASO_amenability']    # 'Probably amenable' | 'Possibly amenable' | 'Unlikely amenable'
resL[0]['STATUS_mis_splice']  # the predicted mis-splicing behind the call
```

Incidence of any variant set of a gene:

```python
from incidence import incidence, expected_patients, amenable_fraction

afL  = [2.1e-4, 5.0e-5, 1.2e-5]   # allele frequency of every P/LP variant of the gene
mild = [2]                        # hypomorphic variants
amen = [0, 1]                     # the variant set to evaluate

ir_gene = incidence(afL, mild_idxL=mild)         # all pathogenic variants
ir_set  = incidence(afL, amen, mild_idxL=mild)   # the set
expected_patients(ir_set, 3_700_000)             # patients per year
amenable_fraction(ir_set, ir_gene)               # fraction of the gene's patients
```

Pass `xlr=True` for chrX genes, where `afL` holds female allele frequencies.

## Examples

Both examples ship with their inputs and expected outputs, and compare the two:

```bash
cd 02_ASO_amenability/example && ./run_example.sh   # 6 variants -> amenability VCF
cd 03_incidence/example       && ./run_example.sh   # 63 variants of 3 genes -> incidence
```

The amenability example annotates two `Probably amenable`, two `Possibly amenable` and
two `Unlikely amenable` variants. The incidence example covers one autosomal recessive
gene (ATAD3A) and two X-linked genes (ATP6AP1, BCAP31), including a gene with no
observed variant, and produces the variant table and all three gene tables.

SpliceAI scores can differ in the third decimal between GPU and CPU, and Jeffreys
intervals in their last digit between SciPy versions; neither changes a call.

## Method summary

- **Pathogenic set** — a variant is kept when at least one germline submission calls it
  pathogenic or likely pathogenic. It is flagged `Ambiguous`, and excluded from the
  incidence estimates, when benign submissions outnumber the pathogenic ones, when VUS
  submissions predominate, when a submission describes low penetrance (from the ClinVar
  terms or from the submission comments), or when its cohort frequency exceeds three
  times the 99th percentile of its gene.
- **Amenability** — SpliceAI (5 kb context, transcript-aware), MaxEntScan and LaBranchoR
  predictions are combined with the VEP consequence. A variant is `Probably amenable`
  when it gains mis-splicing without severe coding or canonical-splice damage,
  `Possibly amenable` when it causes exon skipping or retention, and `Unlikely amenable`
  otherwise.
- **Incidence** — autosomal recessive genes use `(sum AF)^2` minus the combinations
  without the variant set, X-linked genes `0.5 * sum(female AF)`; combinations of two
  hypomorphic alleles are excluded. Expected patients use 3.7M US and 640K UK annual
  births. Credible intervals propagate the Jeffreys posterior of the focal allele count
  through the same functions, recomputing the denominator at the same draw.
- **ASO grouping** — variants sharing an (amenability call, splicing consequence) key are
  treated as addressable by one ASO; a *recurrent* ASO is one expected in >= 0.1 US
  patients per year.

Every quantity is computed per variant-transcript pair, so a variant annotated to two
overlapping genes is scored against each gene's own allele pool.

## Bundled dependencies

`aso_amenability/vendor/` contains the two in-house packages the amenability step needs:
`jkspliceai` (the modified SpliceAI API) in full, and the part of `jklib` it uses. They
are importable as `jklib` and `jkspliceai` without being installed.
