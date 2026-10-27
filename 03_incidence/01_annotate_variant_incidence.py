#!/usr/bin/env python3
"""Per variant-transcript incidence.

For every (gene, transcript) pair a variant is annotated to, the incidence of the
cases carrying it (PR), the expected number of patients (EN), the fraction of the
gene's cases it explains (PRd) and their Jeffreys credible intervals.

A VCF INFO field holds one value per variant and cannot carry a different incidence
for each gene, so the per-pair values are written into the CSQ field, one set per VEP
annotation. The scalar INFO fields are kept for backward compatibility: they hold the
values of the primary pair, plus a packed `<cohort>_NoBiHypo_bygene` field when the
variant belongs to more than one pair.
"""

import argparse
import collections

from incidence import (CALL_RATE_MIN, IR_amenable, IR_pathogenic, beta_quantiles,
                       is_cohort, is_hypomorphic, is_xlr_hypomorphic, newborns,
                       parse_info, xlr_female_af_key)

outline = lambda x: '\t'.join(map(str, x)) + '\n'

gnomad_selectedL = ['AoU_CDRv8_AF',
                    'ukbb_AF',
                    'gnomADv4.1_AF_joint',
                    'gnomADv4.1_AF_joint_afr',
                    'gnomADv4.1_AF_joint_amr',
                    'gnomADv4.1_AF_joint_fin',
                    'gnomADv4.1_AF_joint_nfe',
                    'gnomADv4.1_AF_joint_eas',
                    'gnomADv4.1_AF_joint_sas']

# chrX is evaluated on the female allele frequency of these three cohorts only
xlr_main_keys = ['AoU_CDRv8_AF', 'ukbb_AF', 'gnomADv4.1_AF_joint']

# gnomAD aggregates heterogeneous cohorts and is not used for incidence; only AoU and
# UKBB values are carried per transcript into the CSQ field
csq_keyL = ['AoU_CDRv8_AF', 'ukbb_AF']

gene_fieldL = ('PR_NoBiHypo', 'EN_NoBiHypo', 'PRd_NoBiHypo',
               'PR_NoBiHypo_CIlo', 'PR_NoBiHypo_CIhi',
               'EN_NoBiHypo_CIlo', 'EN_NoBiHypo_CIhi',
               'PRd_NoBiHypo_CIlo', 'PRd_NoBiHypo_CIhi')


def csq_field_names(keyL=csq_keyL, fieldL=gene_fieldL):
    return [k.replace('AF', f) for k in keyL for f in fieldL]


def blank_fields(resD, k):
    for field in gene_fieldL:
        resD[k.replace('AF', field)] = ''


def genekeys_of(vepDL):
    """The variant-transcript pairs a variant is annotated to, in CSQ order."""
    return list(dict.fromkeys(f"{d['SYMBOL']}:{d['Feature']}" for d in vepDL))


def primary_genekey(vepDL):
    """Pair backing the scalar INFO fields: MANE Select first, then alphabetical.

    ClinVar's own GeneSymbol is not used - it is frequently not the HGNC symbol and
    cannot be matched to VEP's SYMBOL in general.
    """
    maneL = [f"{d['SYMBOL']}:{d['Feature']}" for d in vepDL if d.get('MANE_SELECT', '')]
    return sorted(set(maneL))[0] if maneL else sorted(genekeys_of(vepDL))[0]


def get_pathoD(input_vcf):
    """Per (gene, transcript) pool of pathogenic variants, and each cohort's maximum AN."""

    most_freq = collections.defaultdict(lambda: {k: {'varL': []} for k in gnomad_selectedL})
    xlr_pathoD = collections.defaultdict(lambda: {k: {'varL': []} for k in xlr_main_keys})
    AN_maxD = {k: 0 for k in gnomad_selectedL}

    for line in open(input_vcf):

        if line.startswith("##INFO=<ID=CSQ"):
            vephL = line.strip(' \r\n').split("Format: ")[-1].replace('">', '').split('|')
            continue

        if line.startswith("#"):
            continue

        tokL = line.strip(' \r\n').split('\t')

        varstr = f'{tokL[0]}-{tokL[1]}-{tokL[3]}-{tokL[4]}'
        chrom = tokL[0]
        info_txt, vep_txt = tokL[7].split(';CSQ=')
        infoD = parse_info(info_txt)
        vepDL = [{k: v for k, v in zip(vephL, vep.split('|'))} for vep in vep_txt.split(',')]

        if infoD['Ambiguous'] == 'True':
            continue

        genekeyL = genekeys_of(vepDL)
        isHypo = is_hypomorphic(infoD)

        for k in gnomad_selectedL:
            AN_maxD[k] = max(int(infoD[k.replace('AF', 'AN')]), AN_maxD[k])
            if chrom == 'chrX':
                continue  # the AR pool is not used for chrX
            if float(infoD[k]) > 0:
                for genekey in genekeyL:
                    most_freq[genekey][k]['varL'].append({'varstr': varstr,
                                                          'AF': float(infoD[k]),
                                                          'AC': int(infoD[k.replace('AF', 'AC')]),
                                                          'AN': int(infoD[k.replace('AF', 'AN')]),
                                                          'nhomalt': int(infoD[k.replace('AF', 'nhomalt')]),
                                                          'isHypo': isHypo})

        if chrom == 'chrX':
            isHypo_XLR = is_xlr_hypomorphic(infoD)
            for k in xlr_main_keys:
                female_af = float(infoD.get(xlr_female_af_key(k), 0))
                if female_af > 0:
                    for genekey in genekeyL:
                        xlr_pathoD[genekey][k]['varL'].append({'varstr': varstr,
                                                               'female_AF': female_af,
                                                               'AN': int(infoD[k.replace('AF', 'AN')]),
                                                               'isHypo_XLR': isHypo_XLR})

    return most_freq, AN_maxD, xlr_pathoD


def xlr_incidence(infoD, genekey, k, varstr, xlr_pathoD, AN_maxD, isHypo_XLR):
    """chrX: incidence = 0.5 * female allele frequency."""

    resD = {}

    if is_cohort(k):   # call-rate filter (total AN as proxy)
        an = int(infoD[k.replace('AF', 'AN')])
        an_pct = an / AN_maxD[k] if min(an, AN_maxD[k]) > 0 else -1
        if an_pct < CALL_RATE_MIN:
            blank_fields(resD, k)
            return resD

    if infoD['Ambiguous'] == 'True' or isHypo_XLR:
        blank_fields(resD, k)
        return resD

    female_af = float(infoD.get(xlr_female_af_key(k), 0))
    pr = 0.5 * female_af
    nb = newborns(k)

    resD[k.replace('AF', 'PR_NoBiHypo')] = pr
    if nb:
        resD[k.replace('AF', 'EN_NoBiHypo')] = f"{pr * nb:.3f}"

    # PRd: this variant's PR over the gene's non-hypomorphic female pool, restricted to
    # the variants that pass the same call-rate filter
    xlr_varL = xlr_pathoD[genekey][k]['varL']
    if is_cohort(k):
        xlr_varL = [v for v in xlr_varL if v['AN'] > 0 and v['AN'] / AN_maxD[k] >= CALL_RATE_MIN]
    gene_nonhypo_sum = sum(v['female_AF'] for v in xlr_varL if not v['isHypo_XLR'])
    gene_pr = 0.5 * gene_nonhypo_sum
    resD[k.replace('AF', 'PRd_NoBiHypo')] = -1 if gene_pr == 0 else pr / gene_pr

    # CI: propagate this variant's own female allele count
    fem_k = xlr_female_af_key(k)
    fem_ac = int(infoD.get(fem_k.replace('AF', 'AC'), -1))
    fem_an = int(infoD.get(fem_k.replace('AF', 'AN'), -1))
    qlo, qhi = beta_quantiles(fem_ac, fem_an)

    if qlo is None:
        for field in gene_fieldL[3:]:
            resD[k.replace('AF', field)] = ''
        return resD

    qlo_eff = 0.0 if fem_ac == 0 else qlo
    pr_lo, pr_hi = 0.5 * qlo_eff, 0.5 * qhi
    resD[k.replace('AF', 'PR_NoBiHypo_CIlo')] = pr_lo
    resD[k.replace('AF', 'PR_NoBiHypo_CIhi')] = pr_hi
    resD[k.replace('AF', 'EN_NoBiHypo_CIlo')] = f"{pr_lo*nb:.3f}"
    resD[k.replace('AF', 'EN_NoBiHypo_CIhi')] = f"{pr_hi*nb:.3f}"

    if gene_pr == 0:
        resD[k.replace('AF', 'PRd_NoBiHypo_CIlo')] = -1
        resD[k.replace('AF', 'PRd_NoBiHypo_CIhi')] = -1
    else:
        sum_rest = gene_nonhypo_sum - female_af
        gp_lo = 0.5 * (sum_rest + qlo_eff)
        gp_hi = 0.5 * (sum_rest + qhi)
        resD[k.replace('AF', 'PRd_NoBiHypo_CIlo')] = -1 if gp_lo <= 0 else pr_lo / gp_lo
        resD[k.replace('AF', 'PRd_NoBiHypo_CIhi')] = -1 if gp_hi <= 0 else pr_hi / gp_hi

    return resD


def ar_incidence(infoD, genekey, k, varstr, pathoD, AN_maxD):
    """Autosomes: incidence = P(two pathogenic alleles)."""

    resD = {}

    pathovarL = pathoD[genekey][k]['varL']
    if is_cohort(k):   # keep only variants called in >= 30% of the cohort
        pathovarL = [var for var in pathovarL if var['AN'] > 0 and var['AN'] / AN_maxD[k] >= CALL_RATE_MIN]

    patho_afL = [d['AF'] for d in pathovarL]
    currvar_idxL = [i for i, d in enumerate(pathovarL) if d['varstr'] == varstr]
    mildvar_idxL = [i for i, d in enumerate(pathovarL) if d['isHypo']]

    ir_pnh = IR_pathogenic(patho_afL, mildvar_idxL)
    ir_vnh = IR_amenable(patho_afL, currvar_idxL, mildvar_idxL)

    if infoD['Ambiguous'] == 'True' or \
            (is_cohort(k) and infoD[k.replace('AF', 'AN_percent')] < CALL_RATE_MIN):
        blank_fields(resD, k)
        return resD

    nb = newborns(k)
    resD[k.replace('AF', 'PR_NoBiHypo')] = ir_vnh
    if nb:
        resD[k.replace('AF', 'EN_NoBiHypo')] = f"{ir_vnh * nb:.3f}"
    resD[k.replace('AF', 'PRd_NoBiHypo')] = 0 if ir_pnh == 0 else ir_vnh / ir_pnh

    # CI: propagate this variant's own allele count, recomputing the gene denominator
    # at the same draw
    foc_ac = int(infoD.get(k.replace('AF', 'AC'), -1))
    foc_an = int(infoD.get(k.replace('AF', 'AN'), -1))
    qlo, qhi = beta_quantiles(foc_ac, foc_an)

    if (qlo is None) or (len(currvar_idxL) == 0):
        for field in gene_fieldL[3:]:
            resD[k.replace('AF', field)] = ''
        return resD

    ci = currvar_idxL[0]
    qlo_eff = 0.0 if foc_ac == 0 else qlo
    af_lo = list(patho_afL)
    af_lo[ci] = qlo_eff
    af_hi = list(patho_afL)
    af_hi[ci] = qhi

    prv_lo = IR_amenable(af_lo, currvar_idxL, mildvar_idxL)
    prv_hi = IR_amenable(af_hi, currvar_idxL, mildvar_idxL)
    pg_lo = IR_pathogenic(af_lo, mildvar_idxL)
    pg_hi = IR_pathogenic(af_hi, mildvar_idxL)

    resD[k.replace('AF', 'PR_NoBiHypo_CIlo')] = prv_lo
    resD[k.replace('AF', 'PR_NoBiHypo_CIhi')] = prv_hi
    if nb:
        resD[k.replace('AF', 'EN_NoBiHypo_CIlo')] = f"{prv_lo*nb:.3f}"
        resD[k.replace('AF', 'EN_NoBiHypo_CIhi')] = f"{prv_hi*nb:.3f}"
    resD[k.replace('AF', 'PRd_NoBiHypo_CIlo')] = 0 if pg_lo == 0 else prv_lo / pg_lo
    resD[k.replace('AF', 'PRd_NoBiHypo_CIhi')] = 0 if pg_hi == 0 else prv_hi / pg_hi

    return resD


def annotate_variant_incidence(input_vcf, output_vcf):

    pathoD, AN_maxD, xlr_pathoD = get_pathoD(input_vcf)

    out = open(output_vcf, 'w')
    csq_extraL = csq_field_names()

    for line in open(input_vcf):

        if line.startswith("##INFO=<ID=CSQ"):
            vephL = line.strip(' \r\n').split("Format: ")[-1].replace('">', '').split('|')
            line = line.rstrip('\r\n')
            assert line.endswith('">'), 'unexpected CSQ header format'
            line = line[:-2] + '|' + '|'.join(csq_extraL) + '">\n'

        if line.startswith("#"):
            out.write(line)
            continue

        tokL = line.strip(' \r\n').split('\t')

        varstr = f'{tokL[0]}-{tokL[1]}-{tokL[3]}-{tokL[4]}'
        chrom = tokL[0]
        info_txt, vep_txt = tokL[7].split(';CSQ=')
        infoD = parse_info(info_txt)
        vepDL = [{k: v for k, v in zip(vephL, vep.split('|'))} for vep in vep_txt.split(',')]

        genekeyL = genekeys_of(vepDL)
        isHypo = is_hypomorphic(infoD)
        isHypo_XLR = is_xlr_hypomorphic(infoD)

        # gene-independent fields, written once
        infoD['nhomalt_ge_2'] = isHypo_XLR if chrom == 'chrX' else isHypo
        for k in gnomad_selectedL:
            if is_cohort(k):
                an = int(infoD[k.replace('AF', 'AN')])
                infoD[k.replace('AF', 'AN_percent')] = an / AN_maxD[k] if min(an, AN_maxD[k]) > 0 else -1

        resD = {}
        for genekey in genekeyL:
            res = {}
            if chrom == 'chrX':
                for k in gnomad_selectedL:
                    if k not in xlr_main_keys:   # gnomAD subpopulations: not evaluated for XLR
                        blank_fields(res, k)
                for k in xlr_main_keys:
                    res.update(xlr_incidence(infoD, genekey, k, varstr, xlr_pathoD, AN_maxD, isHypo_XLR))
            else:
                for k in gnomad_selectedL:
                    res.update(ar_incidence(infoD, genekey, k, varstr, pathoD, AN_maxD))
            resD[genekey] = res

        # scalar INFO fields: the primary pair, with every pair packed into _bygene
        #   GENE:TX|PR|EN|PRd|PRlo|PRhi|ENlo|ENhi|PRdlo|PRdhi & GENE2:TX2|...
        #   ';' and ',' are INFO separators, so '&' and '|' are used instead
        infoD.update(resD.get(primary_genekey(vepDL), {}))
        if len(genekeyL) > 1:
            for k in gnomad_selectedL:
                packed = []
                for gk in genekeyL:
                    vals = [str(resD[gk].get(k.replace('AF', f), '')) for f in gene_fieldL]
                    if any(v != '' for v in vals):
                        packed.append('|'.join([gk] + vals))
                if packed:
                    infoD[k.replace('AF', 'NoBiHypo_bygene')] = '&'.join(packed)

        # per-transcript values appended to each VEP annotation
        out_entryL = []
        for j, entry in enumerate(vep_txt.split(',')):
            res = resD.get(f"{vepDL[j]['SYMBOL']}:{vepDL[j]['Feature']}", {}) if j < len(vepDL) else {}
            valL = [str(res.get(k.replace('AF', f), '')) for k in csq_keyL for f in gene_fieldL]
            out_entryL.append(entry + '|' + '|'.join(valL))
        vep_txt = ','.join(out_entryL)

        info_txt = ';'.join([f'{k}={v}' for k, v in infoD.items()])
        tokL[7] = ';CSQ='.join([info_txt, vep_txt])
        out.write(outline(tokL))

    out.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', action="store", required=True)
    parser.add_argument('-o', '--output', dest='output', action="store", default="")
    args = parser.parse_args()

    infile = args.input
    outfile = args.output if args.output else infile.replace(".vcf", ".PRd_variant2.bytrx.vcf")

    print("START!")
    annotate_variant_incidence(infile, outfile)
    print("DONE!")
