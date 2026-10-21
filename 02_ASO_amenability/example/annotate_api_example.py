#!/usr/bin/env python3
"""ASO amenability through the module API, for the variant-transcript pairs of a VCF."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from aso_amenability import annotate_variants


def read_vep_vcf(path):
    vep_headerL, recordL = None, []

    for line in open(path):
        if line.startswith('##INFO=<ID=CSQ'):
            vep_headerL = line.strip(' \n\r').split('Format: ')[1][:-2].split('|')
        if line.startswith('#'):
            continue

        tokL = line.strip(' \n\r').split('\t')
        chrom, pos, ref, alt = tokL[0], tokL[1], tokL[3], tokL[4]

        for vep in tokL[7].split(';CSQ=')[1].split(','):
            vepD = dict(zip(vep_headerL, vep.split('|')))
            recordL.append((chrom, pos, ref, alt, vepD['Feature'], vepD))

    return recordL


if __name__ == '__main__':

    infile = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), 'example.vcf')

    recordL = read_vep_vcf(infile)
    resL = annotate_variants(recordL, batch_size=100)

    for (chrom, pos, ref, alt, tid, vepD), resD in zip(recordL, resL):
        print('\t'.join([chrom, pos, ref, alt, vepD['SYMBOL'], tid,
                         resD['ASO_amenability'], resD['STATUS_mis_splice'],
                         resD['spliceAI_molecular_consequence']]))
