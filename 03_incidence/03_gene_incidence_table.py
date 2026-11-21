#!/usr/bin/env python3
"""Per-gene incidence table.

For every gene, the incidence (PR), expected number of patients (EN), amenable
fraction (PRd) and Jeffreys credible interval of the whole pathogenic variant
set and of each ASO-amenable subset:

    patho        all non-hypomorphic pathogenic variants (total disease incidence)
    amen         Probably + Possibly amenable
    prob         Probably amenable
    amen_recurr  amenable variants of a recurrent ASO (EN >= 0.1 patients)
    prob_recurr  the Probably amenable part of those
    amen_topN    the N highest-EN amenable ASO groups
    prob_topN    the N highest-EN Probably amenable ASO groups

An "ASO" is one distinct (ASO_amenability, splicing molecular consequence) group.
"""

import argparse
import collections
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import aso_amenability  # noqa: F401  (puts the bundled jklib on sys.path)
from jklib.genome import loadBlatOutput

from incidence import (IR_amenable, IR_amenable_XLR, IR_pathogenic, IR_pathogenic_XLR,
                       is_hypomorphic, is_xlr_hypomorphic, newborns, parse_info,
                       pathogenic_pool_ci, variant_set_ci, xlr_female_af_key)

TARGET_TISSUES = {'CNS', 'eye', 'liver'}

# genes with an approved / established therapy, excluded from the ASO-target universe
EXC_EXISTING_THERAPY = set("""GAA F8 F11 VWF SLC22A5 ACADM CACNA1F F7 ACSF3 GPR143 BTD UGT1A1 ACADVL
KIAA0586 L1CAM MMACHC CNGB3 ADAMTS13 CYP27A1 BRWD3 G6PC1 RPE65 FAH TRPM1 MME LIPA PHKB
PLG PTS SLC19A3 AKR1D1 F9 C6 POLD3 MSTO1 C7 GRHPR LMF1""".split())

# amenability overwritten to "Unlikely amenable" from experimental evidence
unlikely_confirmedL = """NM_001360.3(DHCR7):c.964-1G>C
NM_006715.4(MAN2C1):c.2303G>A (p.Arg768Gln)
NM_207346.3(TSEN54):c.919G>T (p.Ala307Ser)
NM_001277062.2(MFF):c.-40-842dup
NM_001244710.2(GFPT1):c.*22C>A
NM_176787.5(PIGN):c.2679C>G (p.Ser893Arg)
NM_020320.5(RARS2):c.-2A>G
NM_018706.7(DHTKD1):c.2185G>A (p.Gly729Arg)
NM_001492.6(GDF1):c.401C>G (p.Ser134Trp)""".split('\n')

AMENABLE_CALLS = ('Probably amenable', 'Possibly amenable')
RECURRENT_EN_MIN = 0.1


def load_mane_by_hgnc(path):
    maneD = collections.defaultdict(list)
    header = None
    for line in open(path):
        if line.startswith("#"):
            header = line.rstrip().split('\t')
            continue
        tokD = {k.replace("#", ""): v for k, v in zip(header, line.strip(' \r\n').split('\t'))}
        maneD[tokD['HGNC_ID'].replace("HGNC:", "")].append(tokD)
    return maneD


def load_mane_by_transcript(path):
    manetD = collections.defaultdict(list)
    header = None
    for line in open(path):
        if line.startswith("#"):
            header = line.rstrip().split('\t')
            continue
        tokD = {k.replace("#", ""): v for k, v in zip(header, line.strip(' \r\n').split('\t'))}
        manetD[tokD['Ensembl_nuc']].append(tokD)
    return manetD


def load_aso_target_organs(path):
    """HGNC number -> target organs of the ASO-treatable diseases of that gene."""
    gdD_ASO = collections.defaultdict(list)
    if not path:
        return gdD_ASO

    df = pd.read_csv(path, sep='\t', dtype=str)
    for _, row in df.iterrows():
        hgnc_key = str(row.get('HGNC_num', '')).strip()
        if not hgnc_key or hgnc_key == 'nan':
            continue
        for tissue in str(row.get('GT_Target_Organs', '')).split('|'):
            tissue = tissue.strip()
            if tissue not in TARGET_TISSUES:
                continue
            gdD_ASO[hgnc_key].append({'Disease_classification': tissue,
                                      'Disorder_Name': row.get('Active_Diseases', ''),
                                      'Disorder_OrphaCode': row.get('disease_id', ''),
                                      'DisorderType_name': row.get('DisorderType_name', '')})
    return gdD_ASO


def load_gene_preset(path, maneD):
    """MANE transcripts of the genes to tabulate, keyed by the unversioned Ensembl ID."""
    hgncS = set()
    with open(path) as f:
        hL = f.readline().strip(' \r\n').split('\t')
        for line in f:
            tokD = dict(zip(hL, line.strip(' \r\n').split('\t')))
            hgncS.add(tokD['HGNC_num'])

    presetD = {}
    for hgnc in hgncS:
        for td in maneD.get(hgnc, []):
            presetD[td['Ensembl_nuc'].split('.')[0]] = {
                'tid': td['Ensembl_nuc'], 'hgnc': hgnc, 'mane_symbol': td['symbol'],
                'mane_key': f"{td['symbol']}:{td['Ensembl_nuc']}"}
    return presetD


def gene_field(vepD, k, field):
    """Per-transcript value written into the CSQ field by 01_annotate_variant_incidence.py.

    The scalar INFO field holds one value per variant, so for a variant annotated to
    several genes it carries only the primary pair's value; the CSQ entry being
    processed is the pair itself.
    """
    return vepD.get(k.replace('AF', field), '')


def sorted_indices(afL, idxL):
    l = [afL[i] for i in idxL]
    idx_sortedL = sorted(range(len(l)), key=lambda i: l[i], reverse=True)
    return [idxL[i] for i in idx_sortedL]


def flaten_list(ll):
    resL = []
    for l in ll:
        resL += l
    return resL


def parse_molcons(molcons):
    # use the first molecular consequence, i.e. the one with the biggest delta
    mcL = molcons.split('/')[0].split(':')
    if molcons == '':
        return None
    return ':'.join([mcL[0], mcL[1], mcL[3]])


def collect_gene_variants(input_vcf, presetD, cohortL, filter_ambiguous=True):
    """Per gene and cohort, the allele frequencies / counts of the pathogenic
    variants and the index lists of each amenable subset."""

    def _empty():
        return {k: {'AF_sum_patho': 0, 'AC_sum_patho': 0, 'NofVar_patho': 0, 'NofVar_ClinVar_patho': 0, 'NofVar_nonmild_patho': 0,
                    'AF_sum_amen': 0, 'AC_sum_amen': 0, 'NofVar_amen': 0, 'NofVar_ClinVar_amen': 0, 'NofVar_nonmild_amen': 0,
                    'AF_sum_prob': 0, 'AC_sum_prob': 0, 'NofVar_prob': 0, 'NofVar_ClinVar_prob': 0, 'NofVar_nonmild_prob': 0,
                    'AF_list_patho': [], 'AC_list_patho': [],
                    'AF_list_amen': [], 'AF_idxlist_amen': [], 'AF_list_prob': [], 'AF_idxlist_prob': [],
                    'AF_idxlist_mild': [], 'AF_idxlist_EN': [], 'AF_idxlist_prob_EN': [],
                    'AF_grplist_amen': [],
                    'ASO_set_ClinVar_amen': set(),
                    'ASO_set_ClinVar_prob': set()} for k in cohortL}

    most_freq = collections.defaultdict(_empty)
    chrX_genesetS = set()
    # keep genes without any observed variant
    [most_freq[e] for e in sorted(presetD, key=lambda e: presetD[e]['mane_key'])]
    vep_labelD = {}
    version_mismatchS = set()

    for line in open(input_vcf):

        if line.startswith("##INFO=<ID=CSQ"):
            vephL = line.strip(' \r\n').split("Format: ")[-1].replace('">', '').split('|')
            continue

        if line.startswith("#"):
            continue

        tokL = line.strip(' \r\n').split('\t')
        chrom = tokL[0]

        info_txt, vep_txt = tokL[7].split(';CSQ=')
        infoD = parse_info(info_txt)
        vepDL = [{k: v for k, v in zip(vephL, vep.split('|'))} for vep in vep_txt.split(',')]

        for vepD in vepDL:
            enst = vepD['Feature'].split('.')[0]
            if enst in presetD:
                vep_labelD.setdefault(enst, f"{vepD['SYMBOL']}:{vepD['Feature']}")
                if vepD['Feature'] != presetD[enst]['tid']:
                    version_mismatchS.add((presetD[enst]['tid'], vepD['Feature']))

        if filter_ambiguous and infoD['Ambiguous'] == 'True':
            continue

        isHypo = is_xlr_hypomorphic(infoD) if chrom == 'chrX' else is_hypomorphic(infoD)

        for vepD in vepDL:
            genekey = vepD['Feature'].split('.')[0]
            if genekey not in presetD:
                continue

            if chrom == 'chrX':
                chrX_genesetS.add(genekey)

            if infoD['Name'] in unlikely_confirmedL:
                vepD['ASO_amenability'] = "Unlikely amenable"

            mc = parse_molcons(vepD['spliceAI_molecular_consequence'])
            if mc is None:
                mc = infoD['Name']   # no splicing consequence: the variant is its own ASO group

            is_amen = vepD['ASO_amenability'] in AMENABLE_CALLS
            is_prob = vepD['ASO_amenability'] == 'Probably amenable'

            for k in cohortL:

                d = most_freq[genekey][k]
                d['NofVar_ClinVar_patho'] += 1

                if is_amen:
                    d['NofVar_ClinVar_amen'] += 1
                    d['ASO_set_ClinVar_amen'].add((vepD['ASO_amenability'], mc))
                if is_prob:
                    d['NofVar_ClinVar_prob'] += 1
                    d['ASO_set_ClinVar_prob'].add((vepD['ASO_amenability'], mc))

                # Inclusion gate: the per-variant PR_NoBiHypo of step 01, which is
                # blank exactly when the variant is not usable for this cohort
                # (ambiguous, call rate < 30%, or XLR-hypomorphic).
                if ('AoU' in k or 'ukbb' in k) and gene_field(vepD, k, 'PR_NoBiHypo') == '':
                    continue

                if chrom == 'chrX':
                    female_af_k = xlr_female_af_key(k)
                    raw_af = float(infoD.get(female_af_k, -1))
                    af_val = raw_af if raw_af > 0 else -1
                    ac_val = int(infoD.get(female_af_k.replace('AF', 'AC'), 0))
                else:
                    af_val = float(infoD[k])
                    ac_val = int(infoD[k.replace('AF', 'AC')])

                if af_val <= 0:
                    continue

                d['AF_list_patho'].append(af_val)
                d['AC_list_patho'].append(ac_val)
                d['AF_sum_patho'] += af_val
                d['AC_sum_patho'] += ac_val
                d['NofVar_patho'] += 1
                idx = len(d['AF_list_patho']) - 1

                if isHypo:
                    d['AF_idxlist_mild'].append(idx)
                else:
                    d['NofVar_nonmild_patho'] += 1

                en_str = gene_field(vepD, k, 'EN_NoBiHypo')
                en = float(en_str) if en_str != '' else -1

                if is_amen:
                    d['AF_grplist_amen'].append((idx, vepD['ASO_amenability'], en, mc))
                    if en >= RECURRENT_EN_MIN:
                        d['AF_idxlist_EN'].append((idx, vepD['ASO_amenability'], en, mc))
                    d['AF_idxlist_amen'].append(idx)
                    d['AF_list_amen'].append(af_val)
                    d['AF_sum_amen'] += af_val
                    d['AC_sum_amen'] += ac_val
                    d['NofVar_amen'] += 1
                    if not isHypo:
                        d['NofVar_nonmild_amen'] += 1

                if is_prob:
                    if en >= RECURRENT_EN_MIN:
                        d['AF_idxlist_prob_EN'].append((idx, vepD['ASO_amenability'], en, mc))
                    d['AF_idxlist_prob'].append(idx)
                    d['AF_list_prob'].append(af_val)
                    d['AF_sum_prob'] += af_val
                    d['AC_sum_prob'] += ac_val
                    d['NofVar_prob'] += 1
                    if not isHypo:
                        d['NofVar_nonmild_prob'] += 1

    # rows are labelled with VEP's symbol, which can differ from the MANE v1.0 one
    labelD = {e: vep_labelD.get(e, presetD[e]['mane_key']) for e in presetD}
    metaD = {labelD[e]: presetD[e] for e in presetD}
    most_freq = {labelD[e]: v for e, v in most_freq.items()}
    chrX_genesetS = {labelD[e] for e in chrX_genesetS}

    renamed = sorted((m['mane_key'], lab) for lab, m in metaD.items() if lab.split(':')[0] != m['mane_symbol'])
    print(f"  ENST version differs VEP/MANE : {len(version_mismatchS)} {sorted(version_mismatchS)[:5]}")
    print(f"  symbol differs VEP/MANE       : {len(renamed)} {renamed}")
    return most_freq, chrX_genesetS, metaD


def build_gene_table(most_freq_D, h, manetD, gdD_ASO, metaD, chrX_genesetS=set(),
                     cohortL=('AoU_CDRv8_AF', 'ukbb_AF'), n_amen=10, top_k=None,
                     one_trx_per_gene=True):

    gene_cols = ['Gene:Transcript', 'MANE_Status', 'ASO_target_organs',
                 'Transcript_length', 'NofExon', 'Exon_length', 'Intron_length', 'Intron_percent',
                 'ClinVar_patho_NofVar', 'ClinVar_amen_NofVar', 'ClinVar_amen_NofASO',
                 'ClinVar_prob_NofVar', 'ClinVar_prob_NofASO']

    data_col_suffix = [
        'patho_AC_sum', 'patho_AF_sum', 'patho_NofVar', 'patho_nonHypo_NofVar',
        'patho_PR_NoBiHypo', 'patho_EN_NoBiHypo',
        'amen_AC_sum', 'amen_AF_sum', 'amen_NofVar', 'amen_nonHypo_NofVar', 'amen_NofASO',
        'prob_AC_sum', 'prob_AF_sum', 'prob_NofVar', 'prob_nonHypo_NofVar', 'prob_NofASO',
        'amen_PR_NoBiHypo', 'amen_EN_NoBiHypo', 'amen_PRd_NoBiHypo',
        'prob_PR_NoBiHypo', 'prob_EN_NoBiHypo', 'prob_PRd_NoBiHypo',
        'amen_PR_NoBiHypo_recurr', 'amen_EN_NoBiHypo_recurr', 'amen_PRd_NoBiHypo_recurr', 'amen_nofVar', 'amen_nofASO',
        'prob_PR_NoBiHypo_recurr', 'prob_EN_NoBiHypo_recurr', 'prob_PRd_NoBiHypo_recurr', 'prob_nofVar', 'prob_nofASO',
        f'amen_PR_NoBiHypo_top{n_amen}', f'amen_EN_NoBiHypo_top{n_amen}', f'amen_PRd_NoBiHypo_top{n_amen}', f'amen_nofVar_top{n_amen}', f'amen_nofASO_top{n_amen}',
        f'amen_PR_NoBiHypo_top{n_amen}_prob', f'amen_EN_NoBiHypo_top{n_amen}_prob', f'amen_PRd_NoBiHypo_top{n_amen}_prob', f'amen_nofVar_top{n_amen}_prob', f'amen_nofASO_top{n_amen}_prob',
        f'amen_EN_NoBiHypo_top{n_amen}_list',
        f'prob_PR_NoBiHypo_top{n_amen}', f'prob_EN_NoBiHypo_top{n_amen}', f'prob_PRd_NoBiHypo_top{n_amen}', f'prob_EN_NoBiHypo_top{n_amen}_list',
    ]

    # credible-interval columns, in the same order as ci_vals below
    ci_suffix = ['patho_PR_NoBiHypo_CIlo', 'patho_PR_NoBiHypo_CIhi',
                 'patho_EN_NoBiHypo_CIlo', 'patho_EN_NoBiHypo_CIhi']
    for _base, _suf in [('amen', ''), ('prob', ''), ('amen', '_recurr'), ('prob', '_recurr'),
                        ('amen', f'_top{n_amen}'), ('amen', f'_top{n_amen}_prob'), ('prob', f'_top{n_amen}')]:
        ci_suffix += [f'{_base}_PR_NoBiHypo{_suf}_CIlo', f'{_base}_PR_NoBiHypo{_suf}_CIhi',
                      f'{_base}_EN_NoBiHypo{_suf}_CIlo', f'{_base}_EN_NoBiHypo{_suf}_CIhi',
                      f'{_base}_PRd_NoBiHypo{_suf}_CIlo', f'{_base}_PRd_NoBiHypo{_suf}_CIhi']
    data_col_suffix = data_col_suffix + ci_suffix

    gD = {k: [] for k in cohortL}
    gene_infoD = {}

    for genekey, d in most_freq_D.items():

        if not d:
            continue

        is_xlr = genekey in chrX_genesetS

        tid = metaD[genekey]['tid']
        tD = h[tid][0]
        intron_len = tD['txnLen'] - tD['exnLen']
        hgnc = manetD[tid][0]['HGNC_ID'].split(':')[-1]

        tissueL = []
        for gd in gdD_ASO[hgnc]:
            tissueL.extend(gd['Disease_classification'].split('|'))
        tissueL = sorted({t for t in tissueL if t in TARGET_TISSUES})

        ref_k = cohortL[0]
        gene_infoD[genekey] = (genekey,
                               manetD[tid][0]['MANE_status'],
                               ','.join(tissueL),
                               tD['txnLen'], len(tD['exnList']), tD['exnLen'], intron_len, intron_len / tD['txnLen'],
                               d[ref_k]['NofVar_ClinVar_patho'],
                               d[ref_k]['NofVar_ClinVar_amen'],
                               len(d[ref_k]['ASO_set_ClinVar_amen']),
                               d[ref_k]['NofVar_ClinVar_prob'],
                               len(d[ref_k]['ASO_set_ClinVar_prob']))

        for k in cohortL:

            # ---- ASO groups, keyed by (amenability call, molecular consequence) ----
            group_amen_allvar = collections.defaultdict(list)
            for tup in d[k]['AF_grplist_amen']:
                group_amen_allvar[(tup[1], tup[3])].append(tup)
            n_all_amen_aso = len(group_amen_allvar)
            n_all_prob_aso = len([kk for kk in group_amen_allvar if kk[0] == 'Probably amenable'])

            group_amen = collections.defaultdict(list)
            for tup in d[k]['AF_idxlist_EN']:
                group_amen[(tup[1], tup[3])].append(tup)
            n_recurr_amen_var = len(d[k]['AF_idxlist_EN'])
            n_recurr_amen_aso = len([tup for tup in group_amen.keys()])
            n_recurr_prob_var = len([tup for tup in d[k]['AF_idxlist_EN'] if tup[1] == 'Probably amenable'])
            n_recurr_prob_aso = len([tup for tup in group_amen.keys() if tup[0] == 'Probably amenable'])

            group_amen_sorted = sorted(list(group_amen.values()), key=lambda xL: sum([x[2] for x in xL]), reverse=True)
            group_prob_sorted = sorted(list([v for kk, v in group_amen.items() if kk[0] == 'Probably amenable']),
                                       key=lambda xL: sum([x[2] for x in xL]), reverse=True)

            idx_groupL_amen_all = flaten_list([[v[0] for v in vL] for vL in group_amen_sorted])
            idx_groupL_amen_prob_all = flaten_list([[v[0] for v in vL if v[1] == 'Probably amenable'] for vL in group_amen_sorted])

            group_amen_topnASO = group_amen_sorted[:n_amen]
            idx_groupL_amen_topnASO = flaten_list([[v[0] for v in vL] for vL in group_amen_topnASO])
            idx_groupL_amen_prob_topnASO = flaten_list([[v[0] for v in vL if v[1] == 'Probably amenable'] for vL in group_amen_topnASO])

            n_topn_amen_var = len(idx_groupL_amen_topnASO)
            n_topn_amen_aso = len(group_amen_topnASO)
            n_topn_prob_var = len(idx_groupL_amen_prob_topnASO)
            n_topn_prob_aso = len([vL for vL in group_amen_topnASO if vL[0][1] == 'Probably amenable'])
            grop_amen_topnASO_enlist = ';'.join([','.join([f'{tup[2]}:{tup[1]}' for tup in sorted(tupL, key=lambda x: x[2], reverse=True)]) for tupL in group_amen_topnASO])

            group_prob_topnASO = group_prob_sorted[:n_amen]
            idx_groupL_prob_topnASO = flaten_list([[v[0] for v in vL] for vL in group_prob_topnASO])
            grop_prob_topnASO_enlist = ';'.join([','.join([f'{tup[2]}:{tup[1]}' for tup in sorted(tupL, key=lambda x: x[2], reverse=True)]) for tupL in group_prob_topnASO])

            # ---- incidence of each variant set ----
            IRp = IR_pathogenic_XLR if is_xlr else IR_pathogenic
            IRa = IR_amenable_XLR if is_xlr else IR_amenable

            afL = d[k]['AF_list_patho']
            acL = d[k]['AC_list_patho']
            mild = d[k]['AF_idxlist_mild']

            scen_pools = [d[k]['AF_idxlist_amen'],        # amen
                          d[k]['AF_idxlist_prob'],        # prob
                          idx_groupL_amen_all,            # amen_recurr
                          idx_groupL_amen_prob_all,       # prob_recurr
                          idx_groupL_amen_topnASO,        # amen_top{n}
                          idx_groupL_amen_prob_topnASO,   # amen_top{n}_prob
                          idx_groupL_prob_topnASO]        # prob_top{n}

            ir_pnh = IRp(afL, mild)
            ir_scenL = [IRa(afL, pool, mild) for pool in scen_pools]
            (ir_anh, ir_prob_nh, ir_anh_rec_amen, ir_anh_rec_amen_prob,
             ir_anhtopn_amen, ir_anhtopn_amen_prob, ir_anhtopn_prob) = ir_scenL

            n_us = newborns(k, default=3_700_000)

            en_pnh = f"{n_us * ir_pnh:.3f}"
            (en_anh, en_prob_nh, en_anh_amen, en_anh_amen_prob,
             en_anhtopn_amen, en_anhtopn_amen_prob, en_anhtopn_prob) = [f"{n_us * ir:.3f}" for ir in ir_scenL]

            prd = lambda ir: 0 if ir_pnh == 0 else ir / ir_pnh

            # ---- credible intervals ----
            AN_gene = (d[k]['AC_sum_patho'] / d[k]['AF_sum_patho']) if d[k]['AF_sum_patho'] > 0 else None
            ci_vals = list(pathogenic_pool_ci(afL, acL, mild, IRp, n_us))
            for _pool in scen_pools:
                ci_vals += list(variant_set_ci(afL, acL, _pool, mild, IRp, IRa, n_us, ir_pnh, AN_gene, is_xlr))

            gD[k].append((
                genekey,
                d[k]['AC_sum_patho'], d[k]['AF_sum_patho'], d[k]['NofVar_patho'], d[k]['NofVar_nonmild_patho'], ir_pnh, en_pnh,
                d[k]['AC_sum_amen'], d[k]['AF_sum_amen'], d[k]['NofVar_amen'], d[k]['NofVar_nonmild_amen'], n_all_amen_aso,
                d[k]['AC_sum_prob'], d[k]['AF_sum_prob'], d[k]['NofVar_prob'], d[k]['NofVar_nonmild_prob'], n_all_prob_aso,
                ir_anh, en_anh, prd(ir_anh),
                ir_prob_nh, en_prob_nh, prd(ir_prob_nh),
                ir_anh_rec_amen, en_anh_amen, prd(ir_anh_rec_amen), n_recurr_amen_var, n_recurr_amen_aso,
                ir_anh_rec_amen_prob, en_anh_amen_prob, prd(ir_anh_rec_amen_prob), n_recurr_prob_var, n_recurr_prob_aso,
                ir_anhtopn_amen, en_anhtopn_amen, prd(ir_anhtopn_amen), n_topn_amen_var, n_topn_amen_aso,
                ir_anhtopn_amen_prob, en_anhtopn_amen_prob, prd(ir_anhtopn_amen_prob), n_topn_prob_var, n_topn_prob_aso, grop_amen_topnASO_enlist,
                ir_anhtopn_prob, en_anhtopn_prob, prd(ir_anhtopn_prob), grop_prob_topnASO_enlist,
            ) + tuple(ci_vals))

    # ---- one transcript per gene: MANE Select > MANE Plus Clinical > other ----
    def mane_priority(genekey):
        status = gene_infoD[genekey][1]
        if 'MANE Select' in status and 'Plus' not in status:
            return 0
        elif 'MANE Plus Clinical' in status:
            return 1
        return 2

    gene_hgnc_to_keys = collections.defaultdict(list)
    for genekey in gene_infoD:
        gene_hgnc_to_keys[metaD[genekey]['hgnc']].append(genekey)

    if one_trx_per_gene:
        selected_keys = set()
        for keys in gene_hgnc_to_keys.values():
            selected_keys.add(min(keys, key=mane_priority))
    else:
        selected_keys = set(gene_infoD)

    for k in cohortL:
        gD[k] = [row for row in gD[k] if row[0] in selected_keys]

    # ---- sort by the pathogenic AF sum ----
    ref_order = {row[0]: i for i, row in enumerate(sorted(
        [row for rows in gD.values() for row in rows], key=lambda x: x[2], reverse=True))}
    for k in cohortL:
        gD[k] = sorted(gD[k], key=lambda x: ref_order.get(x[0], 9999))

    dfL = []
    for k in cohortL:
        kp = k.replace('_AF', '')
        cols = ['Gene:Transcript'] + [f'{kp}_{c}' for c in data_col_suffix]
        dfL.append(pd.DataFrame(gD[k][:top_k], columns=cols).set_index('Gene:Transcript'))

    df_ = dfL[0]
    for x in dfL[1:]:
        df_ = df_.join(x, how='outer')

    gene_df = pd.DataFrame([gene_infoD[gk] for gk in selected_keys], columns=gene_cols).set_index('Gene:Transcript')
    df_ = gene_df.join(df_, how='right').reset_index()

    # gene columns first, then each value column followed by its CI pair
    orig_suffix = [c for c in data_col_suffix if not (c.endswith('_CIlo') or c.endswith('_CIhi'))]
    ci_value_suffixes = set(c[:-5] for c in data_col_suffix if c.endswith('_CIlo'))
    desired = []
    for k in cohortL:
        kp = k.replace('_AF', '')
        for suf in orig_suffix:
            col = f'{kp}_{suf}'
            if col in df_.columns:
                desired.append(col)
                if suf in ci_value_suffixes:
                    for tag in ('_CIlo', '_CIhi'):
                        cc = f'{kp}_{suf}{tag}'
                        if cc in df_.columns:
                            desired.append(cc)
    rest = [c for c in df_.columns if c not in gene_cols and c not in desired]
    df_ = df_[gene_cols + desired + rest].copy()
    df_['HGNC_ID'] = ['HGNC:' + metaD[g]['hgnc'] for g in df_['Gene:Transcript']]
    df_['MANE_v1.0_symbol'] = [metaD[g]['mane_symbol'] for g in df_['Gene:Transcript']]
    return df_


if __name__ == "__main__":

    default_annotations = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                       'aso_amenability', 'vendor', 'jkspliceai', 'annotations')

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', action="store", required=True,
                        help="VCF annotated by 01_annotate_variant_incidence.py")
    parser.add_argument('-o', '--output', dest='output', action="store", default="")
    parser.add_argument('-g', '--genes', dest='genes', action="store", required=True,
                        help="TSV with an HGNC_num column: the genes to tabulate (AR/XLR set)")
    parser.add_argument('--aso-organs', dest='aso_organs', action="store", default="",
                        help="TSV with HGNC_num and GT_Target_Organs columns")
    parser.add_argument('--aso-output', dest='aso_output', action="store", default="",
                        help="second TSV holding the ASO-target universe only")
    parser.add_argument('--all-transcripts-output', dest='alltrx_output', action="store", default="",
                        help="third TSV holding every gene-transcript pair")
    parser.add_argument('--mane', dest='mane', action="store",
                        default=os.path.join(default_annotations, 'MANE.GRCh38.v1.0.summary.txt'))
    parser.add_argument('--refflat', dest='refflat', action="store",
                        default=os.path.join(default_annotations, 'gencode.v41.annotation.gtf.refFlat.MANE.v1.arranged.txt'))
    parser.add_argument('-c', '--cohorts', dest='cohorts', action="store", default='AoU_CDRv8_AF,ukbb_AF')
    parser.add_argument('-n', '--top-n', dest='top_n', type=int, action="store", default=10)
    args = parser.parse_args()

    cohortL = args.cohorts.split(',')
    outfile = args.output if args.output else args.input.replace('.vcf', '.gene_incidence.tsv')

    print("START!")
    maneD = load_mane_by_hgnc(args.mane)
    manetD = load_mane_by_transcript(args.mane)
    h = loadBlatOutput(args.refflat)
    gdD_ASO = load_aso_target_organs(args.aso_organs)
    presetD = load_gene_preset(args.genes, maneD)
    print(f"{len(presetD)} MANE transcripts, {len(gdD_ASO)} genes with ASO target organs")

    most_freqD, chrX_genesetS, metaD = collect_gene_variants(args.input, presetD, cohortL)

    def _build(one_trx):
        return build_gene_table(most_freqD, h, manetD, gdD_ASO, metaD, chrX_genesetS=chrX_genesetS,
                                cohortL=cohortL, n_amen=args.top_n, one_trx_per_gene=one_trx)

    if args.alltrx_output:
        # every gene-transcript pair: a gene with both a MANE Select and a MANE Plus
        # Clinical transcript keeps both rows
        df_all = _build(False)
        df_all.to_csv(args.alltrx_output, sep='\t', header=True, index=False)
        print(f"all gene-transcript : {len(df_all)} rows -> {args.alltrx_output}")

    df_out = _build(True)
    df_out.to_csv(outfile, sep='\t', header=True, index=False)
    print(f"one trx per gene    : {len(df_out)} rows -> {outfile}")

    if args.aso_output:
        # ASO-target universe: an annotated target organ, minus the existing-therapy genes
        _sym = df_out['Gene:Transcript'].str.split(':').str[0]
        _organ = df_out['ASO_target_organs'].fillna('').astype(str).str.strip()
        _exc = _sym.isin(EXC_EXISTING_THERAPY) | df_out['MANE_v1.0_symbol'].isin(EXC_EXISTING_THERAPY)
        df_aso = df_out[(_organ != '') & ~_exc].copy()
        df_aso.to_csv(args.aso_output, sep='\t', header=True, index=False)
        print(f"  ASO target        : {len(df_aso)} rows -> {args.aso_output}")

    print("DONE!")
