#!/usr/bin/env python3
"""Drop non-variant rows and keep MANE-transcript, non-intergenic VEP entries."""

import argparse

outline = lambda x: '\t'.join(map(str, x)) + '\n'

intergenicL = ['upstream_gene_variant', 'downstream_gene_variant', 'intergenic_variant']


def curate_vep(infile, outfile, outinter):

    out = open(outfile, 'w')
    out_inter = open(outinter, 'w')

    for line in open(infile):
        if line.startswith("##INFO=<ID=CSQ"):
            vephL = line.rstrip().split('"')[1].split('Format: ')[1].split('|')
            vep_consequence_idx = vephL.index("Consequence")
            vep_mane_select_idx = vephL.index("MANE_SELECT")
            vep_mane_plus_clinical_idx = vephL.index("MANE_PLUS_CLINICAL")

        if line.startswith("#"):
            out.write(line)
            out_inter.write(line)
            continue

        tokL = line.strip(' \r\n').split('\t')
        if tokL[3] == tokL[4]:  # non-variant
            continue

        try:
            info, vep_txt = tokL[7].split(';CSQ=')
        except ValueError:
            print(tokL)
            continue

        vepLL = [vep.split('|') for vep in vep_txt.split(',')]

        vepinterLL = []
        vepintraLL = []
        for vepL in vepLL:
            is_mane = (vepL[vep_mane_select_idx] != '') or (vepL[vep_mane_plus_clinical_idx] != '')
            if not is_mane:
                continue
            if vepL[vep_consequence_idx] not in intergenicL:
                vepintraLL.append(vepL)
            else:
                vepinterLL.append(vepL)

        if len(vepintraLL):
            tokL[7] = info + ';CSQ=' + ','.join(['|'.join(vepL) for vepL in vepintraLL])
            out.write(outline(tokL))

        if len(vepinterLL):
            tokL[7] = info + ';CSQ=' + ','.join(['|'.join(vepL) for vepL in vepinterLL])
            out_inter.write(outline(tokL))

    out.close()
    out_inter.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', action="store", required=True)
    parser.add_argument('-o', '--output', dest='output', action="store", default="")
    parser.add_argument('-s', '--skipped', dest='skipped', action="store", default="")
    args = parser.parse_args()

    infile = args.input
    outfile = args.output if args.output else infile.replace(".vcf", ".MANE.vcf")
    outinter = args.skipped if args.skipped else infile.replace(".vcf", ".MANE.skipped.vcf")

    print("START!")
    curate_vep(infile, outfile, outinter)
    print("DONE!")
