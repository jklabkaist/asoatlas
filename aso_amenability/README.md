# aso_amenability

ASO amenability of a variant-transcript pair, from SpliceAI (5 kb context, modified API),
MaxEntScan, LaBranchoR and the VEP annotation of the variant.

```python
from aso_amenability import annotate_variants

recordL = [('chr1', '2406956', 'C', 'T', 'ENST00000447513.7', vepD)]   # vepD: VEP CSQ fields
resL = annotate_variants(recordL, batch_size=500)

resL[0]['ASO_amenability']   # 'Probably amenable' | 'Possibly amenable' | 'Unlikely amenable'
```

`../02_ASO_amenability/` holds the VCF batch runner and a runnable example.

## Layout

| file | content |
| --- | --- |
| `annotate.py` | `annotate_variants()` / `annotate_variant()`: SpliceAI batch + amenability of each record |
| `utils.py` | amenability calculation: effect aggregation, classification, output fields |
| `coding.py` | coding effect (pLOF consequence, REVEL, coding-sequence indel) |
| `splicing.py` | splicing effect (canonical splice-site damage, mis-splicing, splicing molecular consequence) |
| `config.py` | reference file paths (`ASO_AMEN_REFFLAT`, `ASO_AMEN_LABRANCHOR` override them) |
| `vendor/jkspliceai` | SpliceAI API modified for transcript-aware, masked, MaxEntScan-paired scoring (bundled in full, including its MaxEntScan matrices) |
| `vendor/jklib` | the part of jklib this package needs: `genome.locus` / `genome.loadBlatOutput` (refFlat parsing, 2bit sequence), `utils.margin` / `utils.overlap`, `utils.basic.rc`, and the `path.yaml` config |

`__init__.py` puts `vendor/` on `sys.path`, so the vendored `jklib` / `jkspliceai` keep importing
each other by name; neither has to be installed separately. The SpliceAI model is loaded on the
first call, not at import.

MaxEntScan is taken from `jkspliceai.MaxEntScan` only (`MaxEntScan2`, used by `spliceAI_batch(MES=True)`
and read back from the SpliceAI output text), so jklib's own `splice/` package is not bundled.

## Install

```bash
pip install -r requirements.txt       # third-party packages only
./setup_data.sh <source_dir>          # links the reference files into vendor/ (--copy to copy)
```

`setup_data.sh` lists the reference files and where they come from; they are not distributed
with this repository. Running on GPU needs `CUDA_VISIBLE_DEVICES` and, for TensorFlow 2.9,
`TF_FORCE_GPU_ALLOW_GROWTH=true`.

## Amenability call

| field | meaning |
| --- | --- |
| `STATUS_pLOF` | severe coding consequence (stop gained / frameshift / start lost) |
| `STATUS_missense` | REVEL > 0.5, or an indel inside the coding sequence |
| `STATUS_canon_splice` | damage to a canonical splice site (SpliceAI loss + MaxEntScan + branchpoint) |
| `STATUS_mis_splice` | gain of mis-splicing, or exon skipping / intron retention |
| `spliceAI_molecular_consequence` | predicted splicing alteration and its effect on the ORF |
| `ASO_amenability` | the call below |
| `Affected` | the variant changes the protein product by any of the routes above |

- `Probably amenable`: gain of mis-splicing with no severe coding or canonical-splice damage
- `Possibly amenable`: exon skipping / retention, or a gain of mis-splicing with moderate damage
- `Unlikely amenable`: severe coding or canonical-splice damage, or no predicted mis-splicing
