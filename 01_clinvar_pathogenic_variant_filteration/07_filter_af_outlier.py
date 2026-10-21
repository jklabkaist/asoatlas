#!/usr/bin/env python3
"""Flag variants whose cohort allele frequency is a per-gene outlier (> 3 x 99th
percentile of the P/LP variants of that gene) and mark them ambiguous."""

import argparse
import collections

import numpy as np

outline = lambda x: '\t'.join(map(str, x)) + '\n'

ac_keyL = ['AoU_CDRv8_AC', 'gnomADv4.1_AC_joint', 'ukbb_AC']
outlier_hD = {'AoU_CDRv8_AC': ('over3x99tile_AoU', 'AoU_CDRv8_AF'),
              'ukbb_AC': ('over3x99tile_UKB', 'ukbb_AF'),
              'gnomADv4.1_AC_joint': ('over3x99tile_gnomAD', 'gnomADv4.1_AF_joint')}


def gene_af_percentiles(infile, pct=99):

    dataD = {ac_key: collections.defaultdict(list) for ac_key in ac_keyL}
    vephL = None

    for line in open(infile):
        if line.startswith("##INFO=<ID=CSQ"):
            vephL = line.strip(' \r\n').split("Format: ")[-1].replace('">', '').split('|')

        if line.startswith("#"):
            continue

        if 'Ambiguous=True' in line:
            continue

        tokL = line.strip(' \r\n').split('\t')
        info_txt, vep_txt = tokL[7].split(';CSQ=')
        infoD = {x.split('=')[0]: x.split('=')[1] for x in info_txt.split(';') if '=' in x}

        vepDL = [{k: v for k, v in zip(vephL, vep.split('|'))} for vep in vep_txt.split(',')]
        symbolL = list(set([vepD['SYMBOL'] for vepD in vepDL]))

        for ac_key in ac_keyL:
            if int(infoD[ac_key]) >= 1:
                for symbol in symbolL:
                    dataD[ac_key][symbol].append(float(infoD[ac_key.replace("AC", "AF")]))

    percentiles = {ac_key: collections.defaultdict(lambda: 1) for ac_key in ac_keyL}
    for ac_key, dD in dataD.items():
        for gene, afL in dD.items():
            percentiles[ac_key][gene] = np.percentile(afL, pct) if len(afL) > 0 else 0

    return percentiles, vephL


def filterout_outlier(infile, outfile, percentiles, vephL, fold=3):

    out = open(outfile, 'w')
    for line in open(infile):
        if line.startswith("#"):
            out.write(line)
            continue

        tokL = line.strip(' \r\n').split('\t')
        info_txt, vep_txt = tokL[7].split(';CSQ=')
        infoD = {x.split('=')[0]: x.split('=')[1] for x in info_txt.split(';') if '=' in x}

        if eval(infoD['Ambiguous']):
            for outlier_h, _ in outlier_hD.values():
                infoD[outlier_h] = '-'
        else:
            vepDL = [{k: v for k, v in zip(vephL, vep.split('|'))} for vep in vep_txt.split(',')]
            symbolL = list(set([vepD['SYMBOL'] for vepD in vepDL]))

            for ac_key, (outlier_h, af_key) in outlier_hD.items():
                infoD[outlier_h] = any([float(infoD[af_key]) > fold * percentiles[ac_key][symbol]
                                        for symbol in symbolL])

            if any([infoD[outlier_h] for outlier_h, _ in outlier_hD.values()]) and \
                    (int(infoD['subs_B/LB']) >= 1 or int(infoD['subs_uncertain']) >= 2):
                infoD['Ambiguous'] = True

        tokL[7] = ';'.join([f'{k}={v}' for k, v in infoD.items()]) + ';CSQ=' + vep_txt
        out.write(outline(tokL))

    out.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', action="store", required=True)
    parser.add_argument('-o', '--output', dest='output', action="store", default="")
    parser.add_argument('-p', '--percentile', dest='percentile', type=float, action="store", default=99)
    parser.add_argument('-f', '--fold', dest='fold', type=float, action="store", default=3)
    args = parser.parse_args()

    infile = args.input
    outfile = args.output if args.output else infile.replace('.vcf', '.outfilt.vcf')

    print("START!")
    percentiles, vephL = gene_af_percentiles(infile, pct=args.percentile)
    filterout_outlier(infile, outfile, percentiles, vephL, fold=args.fold)
    print("DONE!")
