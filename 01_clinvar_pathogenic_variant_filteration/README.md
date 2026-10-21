# 01. ClinVar pathogenic variant filtration

Input files (ClinVar FTP, 2025-07 release):

- `variant_summary_2025-07.txt`
- `ClinVarRCVRelease_2025-07.xml`

Requirements: Python >= 3.9, GNU `sort`, Ensembl VEP (singularity image, offline cache, REVEL / NMD / EVE plugins).

## Run order

```bash
# 1. variant_summary -> GRCh38-only, position-sorted VCF
python 01_variant_summary_to_vcf.py \
    -i variant_summary_2025-07.txt \
    -o variant_summary_2025-07.GRCh38.sorted.vcf

# ... left-join gnomAD v4.1 joint / UK Biobank / All of Us INFO fields onto the VCF
#     (not distributed here; see the note at the bottom of 01_variant_summary_to_vcf.py)

# 2. VEP annotation; step 3 is called at the end of the script
sbatch 02_run_vep.sh variant_summary_2025-07.GRCh38.sorted.<joined>.vcf

# 3. drop non-variant rows, keep non-intergenic MANE VEP entries (standalone use)
python 03_curate_vep.py -i variant_summary_2025-07.GRCh38.sorted.<joined>.VEP.noskip.vcf

# 4. RCV XML -> SCV-level TSV (one row per RCV x trait x SCV)
python 04_rcv_xml_to_tsv.py \
    -i ClinVarRCVRelease_2025-07.xml \
    -o ClinVarRCVRelease_2025-07.comment2.tsv

# 5. flag low-penetrance submissions from the submission comments
python 05_annotate_low_penetrance.py \
    -i ClinVarRCVRelease_2025-07.comment2.tsv \
    -o ClinVarRCVRelease_2025-07.comment2.lowP.tsv \
    -k lowP.keywords.txt

# 6. keep variants with >= 1 P/LP germline submission, annotate ambiguity flags
python 06_filter_pathogenic_variants.py \
    -i <annotated>.vcf \
    -s ClinVarRCVRelease_2025-07.comment2.lowP.tsv \
    -o <annotated>.cs_plp.vcf

# 7. flag per-gene allele-frequency outliers and add them to the ambiguous set
python 07_filter_af_outlier.py \
    -i <annotated>.cs_plp.vcf \
    -o <annotated>.cs_plp.outfilt.vcf
```

`02_run_vep.sh` reads its reference paths from the environment:
`VEP_SIF`, `VEP_CACHE_DIR`, `VEP_PLUGIN_DIR`, `VEP_FASTA_DIR`, `REVEL_DIR`
(optional: `VEP_FASTA`, `REVEL_TSV`, `EVE_VCF`, `FORK`, `CONDA_BASE`, `CONDA_ENV`).

## INFO fields added by steps 6-7

| field | meaning |
| --- | --- |
| `subs_P`, `subs_LP`, `subs_P/LP` | number of submissions classified P, LP, and their sum |
| `subs_B`, `subs_LB`, `subs_B/LB` | number of submissions classified B, LB, and their sum |
| `subs_uncertain`, `subs_association`, `subs_no_class` | number of VUS / association / unclassified submissions |
| `subs_low_penetrance` | low-penetrance submissions (ClinVar low-penetrance terms + comment-derived flag from step 5) |
| `Low_penetrance_variant` | `subs_low_penetrance > 0` |
| `Conflicting_but_likely_non_pathogenic` | `subs_P/LP <= subs_B/LB` |
| `Predominant_VUS` | `subs_B/LB > 0` and `subs_uncertain + subs_B/LB > subs_P/LP` |
| `Ambiguous` | any of the three flags above; step 7 adds the frequency-outlier variants |
| `over3x99tile_AoU`, `over3x99tile_UKB`, `over3x99tile_gnomAD` | cohort AF above 3 x the 99th percentile of the P/LP variants of the same gene (`-` for variants already ambiguous) |

Submissions are counted after restricting to traits of type `Disease`, `Finding` or
`PhenotypeInstruction` and dropping somatic / cancer-susceptibility traits.

Step 7 computes the per-gene percentile over the non-ambiguous variants of the input and
re-flags a frequency outlier as `Ambiguous` only when the variant also carries a B/LB
submission or >= 2 VUS submissions.

Step 6 rewrites its fields in place, so it can be re-run after any annotation step. In the
paper run it was applied once to the sorted VCF before the frequency joins, to cut the
number of variants carried through VEP.
