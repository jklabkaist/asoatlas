"""Incidence of an arbitrary variant set, from cohort allele frequencies.

AR genes: incidence = P(two pathogenic alleles), i.e. (sum AF)^2, with the
bi-hypomorphic combinations subtracted. XLR genes: incidence = 0.5 * sum of the
female allele frequencies (hemizygous males).
"""


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
