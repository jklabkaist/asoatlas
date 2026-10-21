#!/usr/bin/env python3
"""Flatten the annotated VCF to one row per variant-transcript pair.

Restricted to the genes of a gene-set TSV, dropping `Ambiguous = True` records. The
incidence columns are filled from the CSQ values of the row's own transcript, not from
the scalar INFO fields, so a variant that belongs to two genes gets each gene's value.
"""

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path

outline = lambda x: '\t'.join(map(str, x)) + '\n'

TARGET_TISSUES = {'CNS', 'eye', 'liver'}

# amenability overwritten from experimental evidence
unlikely_confirmedL = """NM_001360.3(DHCR7):c.964-1G>C
NM_006715.4(MAN2C1):c.2303G>A (p.Arg768Gln)
NM_207346.3(TSEN54):c.919G>T (p.Ala307Ser)
NM_001277062.2(MFF):c.-40-842dup
NM_001244710.2(GFPT1):c.*22C>A
NM_176787.5(PIGN):c.2679C>G (p.Ser893Arg)
NM_020320.5(RARS2):c.-2A>G
NM_018706.7(DHTKD1):c.2185G>A (p.Gly729Arg)
NM_001492.6(GDF1):c.401C>G (p.Ser134Trp)""".split('\n')

# gnomAD aggregates heterogeneous cohorts; its allele counts are kept for hypomorphic
# calling but no gnomAD incidence is reported
drop_gnomad_incidenceL = ('gnomADv4.1_PR_NoBiHypo_joint',
                          'gnomADv4.1_EN_NoBiHypo_joint',
                          'gnomADv4.1_PRd_NoBiHypo_joint')

info_selectL = [
    'Name', 'ClinicalSignificance', 'VariationID', 'RCVaccession', 'PhenotypeList',
    'subs_P/LP', 'subs_P', 'subs_LP', 'subs_B/LB', 'subs_B', 'subs_LB',
    'subs_low_penetrance', 'subs_uncertain', 'subs_association', 'subs_no_class',
    'Low_penetrance_variant', 'Conflicting_but_likely_non_pathogenic', 'Predominant_VUS',
    'over3x99tile_AoU', 'over3x99tile_UKB', 'over3x99tile_gnomAD', 'Ambiguous',
] + [
    'nhomalt_ge_2',
    'AoU_CDRv8_AC', 'AoU_CDRv8_AN', 'AoU_CDRv8_AF', 'AoU_CDRv8_AN_percent',
    'AoU_CDRv8_nhomalt',
    'AoU_CDRv8_AC_male', 'AoU_CDRv8_AN_male', 'AoU_CDRv8_AF_male',
    'AoU_CDRv8_AC_hemi_male', 'AoU_CDRv8_AC_hom_male',
    'AoU_CDRv8_AC_female', 'AoU_CDRv8_AN_female', 'AoU_CDRv8_AF_female',
    'AoU_CDRv8_AC_hom_female', 'AoU_CDRv8_AC_hemi_female',
    'AoU_CDRv8_AC_unknown', 'AoU_CDRv8_AN_unknown',
    'AoU_CDRv8_AC_hom_unknown', 'AoU_CDRv8_AC_hemi_unknown',
    'AoU_CDRv8_PR_NoBiHypo',
    'AoU_CDRv8_PR_NoBiHypo_CIlo', 'AoU_CDRv8_PR_NoBiHypo_CIhi',
    'AoU_CDRv8_EN_NoBiHypo',
    'AoU_CDRv8_EN_NoBiHypo_CIlo', 'AoU_CDRv8_EN_NoBiHypo_CIhi',
    'AoU_CDRv8_PRd_NoBiHypo',
    'AoU_CDRv8_PRd_NoBiHypo_CIlo', 'AoU_CDRv8_PRd_NoBiHypo_CIhi',
] + [
    'ukbb_AC', 'ukbb_AN', 'ukbb_AF', 'ukbb_AN_percent',
    'ukbb_nhomalt',
    'ukbb_AC_male', 'ukbb_AN_male', 'ukbb_AF_male',
    'ukbb_AC_hemi_male', 'ukbb_AC_hom_male',
    'ukbb_AC_female', 'ukbb_AN_female', 'ukbb_AF_female',
    'ukbb_AC_hom_female', 'ukbb_AC_hemi_female',
    'ukbb_AC_unknown', 'ukbb_AN_unknown',
    'ukbb_AC_hom_unknown', 'ukbb_AC_hemi_unknown',
    'ukbb_PR_NoBiHypo',
    'ukbb_PR_NoBiHypo_CIlo', 'ukbb_PR_NoBiHypo_CIhi',
    'ukbb_EN_NoBiHypo',
    'ukbb_EN_NoBiHypo_CIlo', 'ukbb_EN_NoBiHypo_CIhi',
    'ukbb_PRd_NoBiHypo',
    'ukbb_PRd_NoBiHypo_CIlo', 'ukbb_PRd_NoBiHypo_CIhi',
] + [
    'gnomADv4.1_AC_joint', 'gnomADv4.1_AN_joint', 'gnomADv4.1_AF_joint',
    'gnomADv4.1_nhomalt_joint',
    'gnomADv4.1_AC_joint_XY', 'gnomADv4.1_AN_joint_XY', 'gnomADv4.1_AF_joint_XY',
    'gnomADv4.1_nhomalt_joint_XY',
    'gnomADv4.1_AC_joint_XX', 'gnomADv4.1_AN_joint_XX', 'gnomADv4.1_AF_joint_XX',
    'gnomADv4.1_nhomalt_joint_XX',
    'gnomADv4.1_PR_NoBiHypo_joint', 'gnomADv4.1_EN_NoBiHypo_joint',
    'gnomADv4.1_PRd_NoBiHypo_joint',
]
info_selectL = [c for c in info_selectL if c not in drop_gnomad_incidenceL]


def load_hgnc_set(tsv_path):
    """HGNC numbers of the genes to keep."""
    hgncS = set()
    with open(tsv_path) as f:
        header = f.readline().rstrip('\n').split('\t')
        if 'HGNC_num' in header:
            col = header.index('HGNC_num')
        elif 'HGNC' in header:
            col = header.index('HGNC')
        else:
            sys.exit("ERROR: gene TSV must have a column named 'HGNC_num' or 'HGNC'")
        for line in f:
            val = line.rstrip('\n').split('\t')[col].strip()
            if val:
                hgncS.add(val)
    return hgncS


def load_hgnc_to_tissues(tsv_path):
    """HGNC number -> target organs of the ASO-treatable diseases of that gene."""
    hgnc_to_tissues = defaultdict(set)
    with open(tsv_path) as f:
        header = f.readline().rstrip('\n').split('\t')
        if 'HGNC_num' in header:
            hgnc_col = header.index('HGNC_num')
        elif 'HGNC' in header:
            hgnc_col = header.index('HGNC')
        else:
            sys.exit("ERROR: ASO TSV must have a column named 'HGNC_num' or 'HGNC'")
        tissue_col = header.index('GT_Target_Organs')
        for line in f:
            tok = line.rstrip('\n').split('\t')
            hgnc_val = tok[hgnc_col].strip()
            tissue_val = tok[tissue_col].strip() if tissue_col < len(tok) else ''
            if not hgnc_val:
                continue
            for label in tissue_val.split('|'):
                label = label.strip()
                if label in TARGET_TISSUES:
                    hgnc_to_tissues[hgnc_val].add(label)
    return {k: '|'.join(sorted(v)) for k, v in hgnc_to_tissues.items()}


def open_vcf(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def flatten_vcf(infile, outfile, hgnc_filterS, hgnc_to_tissues):

    vephL = None
    csq_override, vep_keep_idx, vephL_out = {}, None, None
    n_vars_in = n_vars_out = n_rows_out = 0

    with open_vcf(infile) as fin, open(outfile, 'w') as fout:
        for line in fin:

            if line.startswith('##INFO=<ID=CSQ'):
                vephL = line.rstrip().split('"')[1].split('Format: ')[1].split('|')
                # the per-transcript fields appended by step 01 all carry 'NoBiHypo' in
                # their name; they overwrite the scalar incidence columns and are not
                # written again as CSQ columns, so the column count stays the same
                info_setS = set(info_selectL)
                csq_override = {c: j for j, c in enumerate(vephL)
                                if 'NoBiHypo' in c and c in info_setS}
                vep_keep_idx = [j for j, c in enumerate(vephL) if 'NoBiHypo' not in c]
                vephL_out = [vephL[j] for j in vep_keep_idx]
                print(f'  CSQ transcript-level fields : {sum(1 for c in vephL if "NoBiHypo" in c)}'
                      f' (overriding {len(csq_override)} table columns)')
                csq_hgnc_idx = vephL.index('HGNC_ID')
                csq_mane_sel_idx = vephL.index('MANE_SELECT')
                csq_amen_idx = vephL.index('ASO_amenability')
                continue

            if line.startswith('#CHROM'):
                if vephL is None:
                    sys.exit('ERROR: ##INFO=<ID=CSQ> header not found before #CHROM')
                vcfhL = line.rstrip().split('\t')[:-1]
                fout.write(outline(vcfhL + ['ASO_target_organs'] + info_selectL + vephL_out))
                continue

            if line.startswith('#'):
                continue

            n_vars_in += 1
            tokL = line.rstrip('\n').split('\t')

            try:
                info, vep_txt = tokL[7].split(';CSQ=')
                infoD = {e.split('=')[0]: e.split('=')[1] for e in info.split(';') if '=' in e}
                info_selectedL = [infoD.get(k, '') for k in info_selectL]
            except Exception as e:
                print(f'WARN: skipping malformed line: {e}', file=sys.stderr)
                continue

            if infoD.get('Ambiguous') == 'True':
                continue

            vepLL = [e.split('|') for e in vep_txt.split(',')]

            has_mane_select = any(csq_mane_sel_idx < len(v) and v[csq_mane_sel_idx] != '' for v in vepLL)
            if has_mane_select:
                vepLL = [v for v in vepLL if csq_mane_sel_idx < len(v) and v[csq_mane_sel_idx] != '']

            wrote = False
            for vepL in vepLL:
                if csq_hgnc_idx >= len(vepL):
                    continue

                if infoD['Name'] in unlikely_confirmedL:
                    vepL[csq_amen_idx] = 'Unlikely amenable'

                hgnc_id = vepL[csq_hgnc_idx].split(':')[-1].strip()
                if hgnc_id not in hgnc_filterS:
                    continue

                row_info = list(info_selectedL)
                for col, j in csq_override.items():
                    if j < len(vepL):
                        row_info[info_selectL.index(col)] = vepL[j]
                vepL_out = [vepL[j] if j < len(vepL) else '' for j in vep_keep_idx]

                fout.write(outline(tokL[:-1] + [hgnc_to_tissues.get(hgnc_id, '')] + row_info + vepL_out))
                n_rows_out += 1
                wrote = True

            if wrote:
                n_vars_out += 1

    print(f'  Variants in      : {n_vars_in:,}')
    print(f'  Variants written : {n_vars_out:,}')
    print(f'  Rows written     : {n_rows_out:,}')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='VCF from 01_annotate_variant_incidence.py')
    parser.add_argument('-g', '--genes', required=True, help='gene-set TSV with an HGNC_num column')
    parser.add_argument('--aso', default='', help='TSV with HGNC_num and GT_Target_Organs columns')
    parser.add_argument('-o', '--output', default='')
    args = parser.parse_args()

    inp = Path(args.input)
    outfile = args.output if args.output else \
        str(inp.parent / (inp.name.replace('.vcf.gz', '').replace('.vcf', '') + '.AR_XLR.tsv'))

    print(f'Gene filter : {args.genes}')
    hgnc_filterS = load_hgnc_set(args.genes)
    print(f'  {len(hgnc_filterS):,} unique HGNC IDs loaded')

    hgnc_to_tissues = load_hgnc_to_tissues(args.aso) if args.aso else {}
    if args.aso:
        print(f'ASO tissue  : {args.aso}')
        print(f'  {len(hgnc_to_tissues):,} genes with ASO tissue annotation'
              f' ({len(set(hgnc_to_tissues) & hgnc_filterS):,} in the gene filter)')

    print(f'Input VCF   : {args.input}')
    print(f'Output TSV  : {outfile}')
    flatten_vcf(args.input, outfile, hgnc_filterS, hgnc_to_tissues)
    print('Done.')
