#!/usr/bin/env python3
"""Keep ClinVar variants carrying at least one P/LP germline submission and
annotate submission-level counts / ambiguity flags into the VCF INFO field."""

import argparse
from collections import defaultdict

outline = lambda x: '\t'.join(map(str, x)) + '\n'

xml_PL = ['PATHOGENIC',
          'Pathogenic',
          'Pathogenic/Likely pathogenic',
          'pathogenic']

xml_LPL = ['LIKELY PATHOGENIC',
           'Likely Pathogenic',
           'Likely pathogenic',
           'likely Pathogenic',
           'likely pathogenic']

xml_patho_lowL = ['Pathogenic, low penetrance',
                  'Likely pathogenic, low penetrance',
                  'Established risk allele']

xml_BL = ['Benign',
          'Benign (Pseudodeficiency Allele)',
          'Benign/Likely benign',
          'benign']

xml_LBL = ['Likely Benign',
           'Likely benign',
           'likely benign',
           'no known pathogenicity']

xml_vusL = ['Uncertain Significance',
            'Uncertain risk allele',
            'Uncertain significance',
            'uncertain significance',
            'conflicting data from submitters']

xml_assoL = ['Affects',
             'Association',
             'association',
             'Likely risk allele',
             'Increased risk allele',
             'other',
             'drug response',
             'risk factor',
             'protective',
             'confers sensitivity']

xml_noclassL = ['not provided',
                'association not found']

xml_germlinemapD = {k: 'P' for k in xml_PL} | \
                   {k: 'LP' for k in xml_LPL} | \
                   {k: 'low_penetrance' for k in xml_patho_lowL} | \
                   {k: 'uncertain' for k in xml_vusL} | \
                   {k: 'association' for k in xml_assoL} | \
                   {k: 'B' for k in xml_BL} | \
                   {k: 'LB' for k in xml_LBL} | \
                   {k: 'no_class' for k in xml_noclassL}

ambig_hL = ['subs_P/LP',
            'subs_P',
            'subs_LP',
            'subs_low_penetrance',
            'subs_uncertain',
            'subs_association',
            'subs_B/LB',
            'subs_B', 'subs_LB',
            'subs_no_class',
            'Low_penetrance_variant', 'Conflicting_but_likely_non_pathogenic', 'Predominant_VUS',
            'Ambiguous']

scv_fieldL = ['TraitType', 'TraitName', 'SCV', 'SCV_GermlineClassification', 'low_penetrance']
proper_traitL = ['Disease', 'Finding', 'PhenotypeInstruction']


def load_variation_ids(infile):
    vidS = set()
    for line in open(infile):
        if line.startswith('#'):
            continue
        info = line.split('\t')[7]
        infoD = {x.split('=')[0]: x.split('=')[1] for x in info.split(';') if '=' in x}
        vidS.add(infoD['VariationID'])
    return vidS


def load_scv_tsv(infile, vidS=None):
    resD = defaultdict(list)
    with open(infile, encoding='utf-8', errors='replace') as f:
        hL = f.readline().strip(' \r\n').split('\t')
        vid_idx = hL.index('VariationID')
        field_idxL = [hL.index(k) for k in scv_fieldL]
        max_idx = max([vid_idx] + field_idxL)

        for line in f:
            tokL = line.rstrip('\r\n').split('\t')
            if len(tokL) <= max_idx:
                continue
            vid = tokL[vid_idx]
            if (vidS is not None) and (vid not in vidS):
                continue
            resD[vid].append({k: tokL[i] for k, i in zip(scv_fieldL, field_idxL)})
    return resD


def filter_cs_plp(infile, outfile, clinvar_scvD):

    out = open(outfile, 'w')
    for line in open(infile):
        if line.startswith('#'):
            out.write(line)
            continue

        tokL = line.strip(' \r\n').split('\t')
        infoD = {x.split('=')[0]: x.split('=')[1] for x in tokL[7].split(';') if '=' in x}

        dL = [d for d in clinvar_scvD[infoD['VariationID']]]

        # filter traitName
        properTraitdL = [d for d in dL if d['TraitType'] in proper_traitL]
        # filter cancer risk
        properTraitdL = [d for d in properTraitdL
                         if not (('cancer' in d['TraitName'].lower()) and
                                 ('hereditary' not in d['TraitName'].lower()) and
                                 ('somatic' in d['TraitName'].lower() or 'susceptibility' in d['TraitName'].lower()))
                         ]

        # aggregate traits for same SCV
        proper_scvD = defaultdict(list)
        for d in properTraitdL:
            proper_scvD[d['SCV']].append(d)

        # check SCV germlineClassification
        scvD = {k.replace('subs_', ''): [] for k in ambig_hL if k.startswith('subs')}
        lowp_cmt = 0
        for scv, dl_ in proper_scvD.items():
            d = dl_[0]
            if d['low_penetrance'] == "TRUE":
                lowp_cmt = lowp_cmt + 1
            cs = d['SCV_GermlineClassification']
            germline = xml_germlinemapD[cs]
            scvD[germline].append(d)

        if len(scvD['P']) + len(scvD['LP']) == 0:  # only variant with P/LP submission
            continue

        for ambig_h in ambig_hL:  # number of subs for each type
            if ambig_h.startswith('subs'):
                infoD[ambig_h] = len(scvD[ambig_h.replace('subs_', '')])

        infoD['subs_low_penetrance'] = int(infoD['subs_low_penetrance']) + lowp_cmt
        infoD['subs_P/LP'] = int(infoD['subs_P']) + int(infoD['subs_LP'])
        infoD['subs_B/LB'] = int(infoD['subs_B']) + int(infoD['subs_LB'])

        # ambiguity determination
        infoD['Low_penetrance_variant'] = int(infoD['subs_low_penetrance']) > 0
        infoD['Conflicting_but_likely_non_pathogenic'] = int(infoD['subs_P/LP']) <= int(infoD['subs_B/LB'])
        infoD['Predominant_VUS'] = int(infoD['subs_B/LB']) > 0 and \
            int(infoD['subs_uncertain']) + int(infoD['subs_B/LB']) > int(infoD['subs_P/LP'])
        infoD['Ambiguous'] = infoD['Low_penetrance_variant'] or \
            infoD['Conflicting_but_likely_non_pathogenic'] or infoD['Predominant_VUS']

        if ';CSQ=' in tokL[7]:
            info, vep = tokL[7].split(';CSQ=')
            existing_parts = [x for x in info.split(';') if x and x.split('=')[0] not in ambig_hL]
            new_fields = [f'{k}={infoD[k]}' for k in ambig_hL]
            tokL[7] = ';'.join(existing_parts + new_fields) + ';CSQ=' + vep
        else:
            existing_parts = [x for x in tokL[7].split(';') if x and x.split('=')[0] not in ambig_hL]
            new_fields = [f'{k}={infoD[k]}' for k in ambig_hL]
            tokL[7] = ';'.join(existing_parts + new_fields)

        out.write(outline(tokL))

    out.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', action="store", required=True)
    parser.add_argument('-s', '--scv', dest='scv', action="store", required=True,
                        help="SCV-level TSV with the low_penetrance column (step 05 output)")
    parser.add_argument('-o', '--output', dest='output', action="store", default="")
    args = parser.parse_args()

    infile = args.input
    outfile = args.output if args.output else infile.replace('.vcf', '.cs_plp.vcf')

    print("START!")
    vidS = load_variation_ids(infile)
    clinvar_scvD = load_scv_tsv(args.scv, vidS)
    filter_cs_plp(infile, outfile, clinvar_scvD)
    print("DONE!")
