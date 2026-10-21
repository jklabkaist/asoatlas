# incidence

Incidence of an arbitrary variant set, and its Jeffreys credible interval.

```python
from incidence import incidence, expected_patients, amenable_fraction

afL  = [2.1e-4, 5.0e-5, 1.2e-5]     # allele frequency of every P/LP variant of the gene
mild = [2]                          # hypomorphic variants
amen = [0, 1]                       # the variant set to evaluate

ir_gene = incidence(afL, mild_idxL=mild)                  # all pathogenic variants
ir_set  = incidence(afL, amen, mild_idxL=mild)            # the variant set
expected_patients(ir_set, 3_700_000)                      # patients per year
amenable_fraction(ir_set, ir_gene)                        # fraction of the gene's cases
```

Pass `xlr=True` for chrX genes, where `afL` holds female allele frequencies and the
incidence is linear in AF.

| file | content |
| --- | --- |
| `rates.py` | `IR_pathogenic` / `IR_amenable` (+ `_XLR`), and the `incidence()` wrapper |
| `ci.py` | Jeffreys Beta posterior: `beta_quantiles`, `pathogenic_pool_ci`, `variant_set_ci` |
| `fields.py` | VCF INFO conventions shared by both levels (cohort keys, hypomorphic rules, call rate, newborns per cohort) |

[`../03_incidence`](../03_incidence) applies this to a VCF, per variant and per gene.
