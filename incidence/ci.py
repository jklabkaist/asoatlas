"""Bayesian (Jeffreys Beta-binomial) 95% credible intervals for the incidences.

The Jeffreys posterior of the focal allele count,
    q ~ Beta(AC + 1/2, AN - AC + 1/2),
is propagated through the same deterministic incidence functions. Every quantity
is monotone increasing in the focal frequency, so the interval is the function
evaluated at the 2.5% / 97.5% posterior quantiles (pushforward, no Monte-Carlo
error). For the amenable fraction (PRd) the denominator is recomputed at the SAME
draw, which preserves the numerator-denominator dependence.

ZERO-COUNT RULE: an allele pool with no observed allele is left-censored, not an
exact zero. Its lower limit is forced to 0 instead of the strictly positive 2.5%
Jeffreys quantile, so the interval is one-sided [0, upper].
"""

import numpy as np
from scipy.stats import beta as _beta_dist

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
