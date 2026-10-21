#!/usr/bin/env python3
"""Add the ASO amenability fields to the CSQ field of a VEP-annotated VCF."""

import argparse
import collections
import datetime
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from aso_amenability import amenability_headerL, annotate_variants

outline = lambda x: '\t'.join(map(str, x)) + '\n'


def extend_csq_header(line, amen_headerL):
    left, fmt = line.strip(' \n\r').split('Format: ')
    vep_headerL = fmt[:-2].split('|')
    newline = 'Format: '.join([left, '|'.join(vep_headerL + amen_headerL)]) + '">' + '\n'
    return vep_headerL, newline


def write_batch(out, recordL, keyL, vep_headerL, batch_size, top_k):

    resL = annotate_variants(recordL, batch_size=batch_size, top_k=top_k)

    compD = collections.OrderedDict()
    for (key, vepD), resD in zip(keyL, resL):
        csq_txt = '|'.join([str(x) for x in [vepD[k] for k in vep_headerL] +
                                            [resD[k] for k in amenability_headerL]])
        compD.setdefault(key, []).append(csq_txt)

    for key, csq_txtL in compD.items():
        out.write(key + ';CSQ=' + ','.join(csq_txtL) + '\n')


def annotate_vcf(infile, outfile, batch_size=500, top_k=5):

    out = open(outfile, 'w')
    vep_headerL = None
    recordL, keyL = [], []

    for line in open(infile):

        if line.startswith("##INFO=<ID=CSQ"):
            vep_headerL, newline = extend_csq_header(line, amenability_headerL)
            out.write(newline)
            continue

        if line.startswith('#'):
            out.write(line)
            continue

        if len(recordL) >= batch_size:
            write_batch(out, recordL, keyL, vep_headerL, batch_size, top_k)
            recordL, keyL = [], []

        tokL = line.strip(' \n\r').split('\t')
        chrom, pos, ID, ref, alt, qual, filt = tokL[:7]

        try:
            info_txt, vep_txt = tokL[7].split(';CSQ=')
        except ValueError:
            print(f"[{datetime.datetime.now()}]  --- No CSQ: {chrom}-{pos}-{ref}-{alt}", file=sys.stderr)
            continue

        vepDL = [dict(zip(vep_headerL, vep.split('|'))) for vep in vep_txt.split(',')]
        if len(vepDL) <= 0:
            continue

        key = '\t'.join([chrom, pos, ID, ref, alt, qual, filt, info_txt])
        for vepD in vepDL:
            recordL.append((chrom, pos, ref, alt, vepD['Feature'], vepD))
            keyL.append((key, vepD))

    if len(recordL):
        write_batch(out, recordL, keyL, vep_headerL, batch_size, top_k)

    out.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", action='store', default=None)
    parser.add_argument("-I", "--input_list", dest="input_list", action='store', default=None,
                        help="text file listing one input VCF per line")
    parser.add_argument("-o", "--output", dest='output', action='store', default=None)
    parser.add_argument("-b", "--batch_size", dest="batch_size", type=int, action="store", default=500)
    parser.add_argument("-k", "--top_k", dest="top_k", type=int, action="store", default=5)
    args = parser.parse_args()

    if args.input:
        infileL = [args.input]
    elif args.input_list:
        infileL = [line.strip(' \r\n') for line in open(args.input_list) if line.strip()]
    else:
        parser.error("-i or -I is required")

    for infile in infileL:
        outfile = args.output if (args.output and args.input) else infile.replace('.vcf', '.amenability.vcf')
        print(f"[{datetime.datetime.now()}]  - Start: {infile}")
        annotate_vcf(infile, outfile, batch_size=args.batch_size, top_k=args.top_k)
        print(f"[{datetime.datetime.now()}]  - End: {outfile}")

    print("DONE!!")
