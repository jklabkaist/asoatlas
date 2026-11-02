"""Incidence of a variant set, and its Jeffreys credible interval.

Used by the variant-level (01) and gene-level (03) scripts of this folder: the VCF
INFO conventions both read, the AR / XLR incidence model, and the credible intervals.
"""

import numpy as np
from scipy.stats import beta as _beta_dist

# ── VCF INFO conventions ──────────────────────────────────────────────────────
COHORT_NEWBORNS = {'AoU_CDRv8_AF': 3_700_000,   # All of Us: annual US newborns
                   'ukbb_AF': 690_000}          # UK Biobank: annual UK newborns

CALL_RATE_MIN = 0.3                             # minimum AN / max(AN) of a cohort

AUTOSOMAL_NHOMALT_KEYS = ['gnomADv4.1_nhomalt_joint', 'ukbb_nhomalt', 'AoU_CDRv8_nhomalt']


def newborns(af_key, default=0):
    for k, n in COHORT_NEWBORNS.items():
        if k.split('_')[0] in af_key:
            return n
    return default


def is_cohort(af_key):
    return ('AoU' in af_key) or ('ukbb' in af_key)


def parse_info(info_txt):
    return {x.split('=')[0]: x.split('=')[1] for x in info_txt.split(';')}


def safe_int(infoD, key):
    """int value of a count field, or None when the field is absent / negative."""
    v = int(infoD.get(key, -1))
    return v if v >= 0 else None


def is_hypomorphic(infoD):
    """Autosomal: >= 2 homozygotes observed in any cohort."""
    return any([int(infoD.get(k, -1)) >= 2 for k in AUTOSOMAL_NHOMALT_KEYS
                if int(infoD.get(k, -1)) >= 0])


def is_xlr_hypomorphic(infoD):
    """chrX: >= 2 hemizygous males (gnomAD XY / UKBB / AoU) or female homozygotes."""
    gnomad_hemi = safe_int(infoD, 'gnomADv4.1_nhomalt_joint_XY')
    gnomad_fhom = safe_int(infoD, 'gnomADv4.1_nhomalt_joint_XX')
    gnomad_hypo = (gnomad_hemi is not None and gnomad_hemi >= 2) or \
                  (gnomad_fhom is not None and gnomad_fhom >= 2)

    def _eff_count_hypo(hemi_k, fhom_k, uhom_k):
        # UKBB / AoU carry no nhomalt field: male_hemi + 0.5*female_hom + 0.5*unknown_hom
        hemi = safe_int(infoD, hemi_k)
        if hemi is None:
            return False  # source data unavailable
        fhom = safe_int(infoD, fhom_k) or 0
        uhom = safe_int(infoD, uhom_k) or 0
        return (hemi + 0.5 * fhom + 0.5 * uhom) >= 2

    ukbb_hypo = _eff_count_hypo('ukbb_AC_hemi_male', 'ukbb_AC_hom_female', 'ukbb_AC_hom_unknown')
    aou_hypo = _eff_count_hypo('AoU_CDRv8_AC_hemi_male', 'AoU_CDRv8_AC_hom_female', 'AoU_CDRv8_AC_hom_unknown')

    return gnomad_hypo or ukbb_hypo or aou_hypo


def xlr_female_af_key(af_key):
    if is_cohort(af_key):
        return af_key.replace('_AF', '_AF_female')   # AoU_CDRv8_AF_female / ukbb_AF_female
    return af_key + '_XX'                            # gnomADv4.1_AF_joint_XX


# ── incidence of a variant set ────────────────────────────────────────────────
def IR_pathogenic(afL, mild_idxL=[]):
    # P(at least two pathogenic alleles), excluding bi-hypomorphic combinations
    return sum(afL)**2 - (sum([afL[j] for j in mild_idxL])**2)


def IR_amenable(afL, idxL, mild_idxL=[], mild_amen=True):
    # P(at least one allele of the variant set idxL, with a second pathogenic allele)
    n = len(afL)
    idx_nonamenL = [i for i in range(n) if i not in idxL]

    if len(mild_idxL) and mild_amen:
        idx_nonamen_or_mildL = list(set(idx_nonamenL + mild_idxL))
        idx_nonamen_sevL = [i for i in idx_nonamenL if i not in mild_idxL]
        idx_amen_mildL = [i for i in idxL if i in mild_idxL]

        return sum(afL)**2 - (sum([afL[i] for i in idx_nonamen_or_mildL])**2) + \
            (2 * sum([afL[i] for i in idx_nonamen_sevL]) * sum([afL[i] for i in idx_amen_mildL]))

    elif len(mild_idxL) and not mild_amen:  # omit mild_amenable only cases
        idx_nonamen_or_mildL = list(set(idx_nonamenL + mild_idxL))

        return sum(afL)**2 - (sum([afL[i] for i in idx_nonamen_or_mildL])**2)

    else:
        return sum(afL)**2 - (sum([afL[i] for i in idx_nonamenL])**2)


def IR_pathogenic_XLR(afL, mild_idxL=[]):
    # afL holds female allele frequencies; linear in AF (hemizygous males)
    mild_set = set(mild_idxL)
    return sum(0.5 * af for i, af in enumerate(afL) if i not in mild_set)


def IR_amenable_XLR(afL, idxL, mild_idxL=[]):
    mild_set = set(mild_idxL)
    amen_set = set(idxL)
    return sum(0.5 * af for i, af in enumerate(afL) if i in amen_set and i not in mild_set)


def incidence(afL, idxL=None, mild_idxL=(), xlr=False, mild_amen=True):
    """Incidence of the variant set `idxL` (all pathogenic variants when None).

    afL       : allele frequency of every pathogenic variant of the gene
                (female allele frequency when xlr)
    idxL      : indices of the variant set to evaluate
    mild_idxL : indices of the hypomorphic (mild) variants
    """
    mild_idxL = list(mild_idxL)

    if idxL is None:
        return IR_pathogenic_XLR(afL, mild_idxL) if xlr else IR_pathogenic(afL, mild_idxL)

    if xlr:
        return IR_amenable_XLR(afL, list(idxL), mild_idxL)
    return IR_amenable(afL, list(idxL), mild_idxL, mild_amen=mild_amen)


def expected_patients(ir, n_newborns):
    return ir * n_newborns


def amenable_fraction(ir_set, ir_patho):
    return 0 if ir_patho == 0 else ir_set / ir_patho


# ── credible intervals ────────────────────────────────────────────────────────
JEFFREYS_A, JEFFREYS_B = 0.5, 0.5
CI_LO, CI_HI = 0.025, 0.975


def beta_quantiles(ac, an, lo=CI_LO, hi=CI_HI, a=JEFFREYS_A, b=JEFFREYS_B):
    """Posterior quantiles of q ~ Beta(AC+a, AN-AC+b); (None, None) if AN invalid."""
    try:
        ac = float(ac)
        an = float(an)
    except Exception:
        return None, None
    if not np.isfinite(an) or an <= 0 or not np.isfinite(ac) or ac < 0:
        return None, None
    aa = ac + a
    bb = max(an - ac, 0.0) + b
    return float(_beta_dist.ppf(lo, aa, bb)), float(_beta_dist.ppf(hi, aa, bb))


def scale_pool(afL, pool_set, f):
    return [af * f if i in pool_set else af for i, af in enumerate(afL)]


def pathogenic_pool_ci(afL, acL, mild_idxL, IRp, n_newborns):
    """CI of the pathogenic incidence, propagating the non-hypomorphic pool count."""
    nonhypo = [i for i in range(len(afL)) if i not in set(mild_idxL)]
    af_pool = sum(afL[i] for i in nonhypo)
    ac_pool = sum(acL[i] for i in nonhypo)
    if af_pool <= 0 or ac_pool <= 0:
        return (0.0, 0.0, "0.000", "0.000")

    qlo, qhi = beta_quantiles(ac_pool, ac_pool / af_pool)
    if qlo is None:
        return (0.0, 0.0, "0.000", "0.000")

    ps = set(nonhypo)
    pr_lo = IRp(scale_pool(afL, ps, qlo / af_pool), mild_idxL)
    pr_hi = IRp(scale_pool(afL, ps, qhi / af_pool), mild_idxL)
    return (pr_lo, pr_hi, f"{n_newborns*pr_lo:.3f}", f"{n_newborns*pr_hi:.3f}")


def variant_set_ci(afL, acL, pool, mild_idxL, IRp, IRa, n_newborns, patho_pt, AN_gene, is_xlr):
    """CI of one variant-set scenario: (PR_lo, PR_hi, EN_lo, EN_hi, PRd_lo, PRd_hi)."""
    pool = list(pool)
    af_pool = sum(afL[i] for i in pool) if pool else 0.0
    ac_pool = sum(acL[i] for i in pool) if pool else 0.0

    if af_pool > 0 and ac_pool > 0:
        qlo, qhi = beta_quantiles(ac_pool, ac_pool / af_pool)
        if qlo is None:
            return (0.0, 0.0, "0.000", "0.000", 0.0, 0.0)
        ps = set(pool)
        a_lo = scale_pool(afL, ps, qlo / af_pool)
        a_hi = scale_pool(afL, ps, qhi / af_pool)
        pr_lo = IRa(a_lo, pool, mild_idxL)
        pr_hi = IRa(a_hi, pool, mild_idxL)
        pgl = IRp(a_lo, mild_idxL)
        pgh = IRp(a_hi, mild_idxL)
        prd_lo = 0 if pgl == 0 else pr_lo / pgl
        prd_hi = 0 if pgh == 0 else pr_hi / pgh
    else:
        # zero-count rule: one-sided [0, upper] at the gene's effective allele number
        pr_lo, prd_lo = 0.0, 0.0
        if AN_gene and AN_gene > 0:
            _, qz = beta_quantiles(0.0, AN_gene)
            if is_xlr:
                pr_hi = 0.5 * qz
            else:
                Qtot = sum(afL)
                pr_hi = qz * (2 * Qtot - qz) if (2 * Qtot - qz) > 0 else qz * qz
        else:
            pr_hi = 0.0
        prd_hi = 0.0 if (patho_pt + pr_hi) <= 0 else pr_hi / (patho_pt + pr_hi)

    return (pr_lo, pr_hi, f"{n_newborns*pr_lo:.3f}", f"{n_newborns*pr_hi:.3f}", prd_lo, prd_hi)
