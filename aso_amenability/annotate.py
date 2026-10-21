"""ASO amenability annotation of variant-transcript pairs."""

from jkspliceai import spliceAI_batch
from jkspliceai.utils import filter_res_dict, get_resD_masked, outform_vcf_info_refalt

from .splicing import h
from .utils import amenability_result, coding_only_result, log

MAX_SPLICEAI_INDEL = 4500


def transcript_record(tid):
    try:
        return h[tid][0]
    except (KeyError, IndexError):
        raise KeyError(f"transcript is not in the MANE refFlat annotation: {tid}")


def _chunks(idxL, size):
    for i in range(0, len(idxL), size):
        yield idxL[i:i + size]


def annotate_variants(records, batch_size=500, top_k=5, delta_th=0.01):
    """Evaluate ASO amenability for each (chrom, pos, ref, alt, transcript_id, vepD) record.

    vepD is the VEP annotation of that variant-transcript pair as a {field: value} dict;
    'Consequence', 'REVEL', 'EXON', 'CDS_position' and 'Protein_position' are used.

    Returns one result dict per record, keyed by `utils.amenability_headerL`.
    """

    recordL = list(records)
    resL = [None] * len(recordL)

    sai_idxL = []
    for i, (chrom, pos, ref, alt, tid, vepD) in enumerate(recordL):
        # SpliceAI is not evaluated for variants longer than its sequence context
        if len(ref) > MAX_SPLICEAI_INDEL or len(alt) > MAX_SPLICEAI_INDEL:
            resL[i] = coding_only_result(chrom, pos, ref, alt, vepD)
        else:
            sai_idxL.append(i)

    for idxL in _chunks(sai_idxL, batch_size):

        varL = [f"{recordL[i][0]}-{recordL[i][1]}-{recordL[i][2]}-{recordL[i][3]}" for i in idxL]
        tidL = [recordL[i][4] for i in idxL]

        resD = spliceAI_batch(varL, tidL=tidL, verbose=False, MES=True)
        resD = filter_res_dict(varL, resD, tidL=tidL, delta_th=delta_th)
        maskD = get_resD_masked(varL, resD, tidL, delta_th=delta_th, top_k=top_k)
        sai_outD = outform_vcf_info_refalt(varL, resD, tidL, delta_th=delta_th, top_k=top_k, sep='|')

        for i, var, tid in zip(idxL, varL, tidL):
            chrom, pos, ref, alt, _, vepD = recordL[i]
            key = f'{var};{tid}'

            if key not in sai_outD:
                log(f" --- Skipped: {key}")
                resL[i] = coding_only_result(chrom, pos, ref, alt, vepD)
                continue

            resL[i] = amenability_result(chrom, pos, ref, alt, vepD,
                                         transcript_record(tid), sai_outD[key], maskD[key])

    return resL


def annotate_variant(chrom, pos, ref, alt, tid, vepD, top_k=5, delta_th=0.01):
    return annotate_variants([(chrom, pos, ref, alt, tid, vepD)], top_k=top_k, delta_th=delta_th)[0]
