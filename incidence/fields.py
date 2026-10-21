"""VCF INFO conventions shared by the variant-level and gene-level incidence steps."""

COHORT_NEWBORNS = {'AoU_CDRv8_AF': 3_700_000,   # All of Us: annual US newborns
                   'ukbb_AF': 640_000}          # UK Biobank: annual UK newborns

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
