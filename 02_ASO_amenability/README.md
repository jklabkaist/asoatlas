# 02. ASO amenability annotation

Adds the ASO amenability fields of [`aso_amenability`](../aso_amenability) to the `CSQ` field of a
VEP-annotated VCF: one entry per variant-transcript pair, 33 fields appended to the VEP fields.

Input: the VCF produced by step 01 (VEP-annotated, MANE transcripts, P/LP variants).

```bash
# single VCF
python annotate_vcf.py -i <input>.vcf -o <input>.amenability.vcf -b 500

# several VCFs, one path per line
python annotate_vcf.py -I vcf_list.txt -b 500

# on a GPU node
sbatch slurm_annotate_vcf.sh <input>.vcf
```

`-b` is the number of variant-transcript pairs sent to SpliceAI at a time (500 on an A6000).
`CUDA_VISIBLE_DEVICES` selects the GPU; leave it empty to run on CPU.

## Appended CSQ fields

| field | content |
| --- | --- |
| `spliceAI_5k_transcriptID` | transcript the SpliceAI scores were computed on |
| `spliceAI_5k_raw_{relpos,score,MES}_{AG,AL,DG,DL}` | top-5 acceptor/donor gain/loss predictions, with the MaxEntScan score of each site |
| `spliceAI_5k_masked_{relpos,score,MES}_{AG,AL,DG,DL}` | the same, after masking predictions that agree with the annotated splice sites |
| `spliceAI_molecular_consequence` | predicted splicing alteration (pseudoexon / exitron / exon skipping / truncation / retention) and its effect on the ORF |
| `spliceAI_damage2orf_by_splicing` | `ORFdamaged` when the alteration gains a PTC or shifts the frame |
| `STATUS_pLOF`, `STATUS_missense`, `STATUS_canon_splice`, `STATUS_mis_splice` | the four effect evaluations behind the call |
| `ASO_amenability` | `Probably amenable` / `Possibly amenable` / `Unlikely amenable` |
| `Affected` | the variant changes the protein product |

## Example

```bash
cd example && ./run_example.sh          # example.vcf -> example.amenability.vcf
python annotate_api_example.py          # same six variants through the module API
```

`run_example.sh` compares its output with `example/expected/example.amenability.vcf`, the annotated
VCF these six variants give (two `Probably amenable`, two `Possibly amenable`, two `Unlikely amenable`).
SpliceAI scores can differ in the third decimal between GPU and CPU; the amenability calls do not.
