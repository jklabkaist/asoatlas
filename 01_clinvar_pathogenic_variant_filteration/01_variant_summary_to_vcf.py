#!/usr/bin/env python3
"""ClinVar variant_summary.txt -> position-sorted single-assembly VCF."""

import argparse
import os
import subprocess

INFO_KEYS = [
    'Type', 'Name', 'GeneID', 'GeneSymbol', 'HGNC_ID', 'ClinicalSignificance', 'ClinSigSimple', 'LastEvaluated',
    'RCVaccession', 'PhenotypeIDS', 'PhenotypeList', 'Origin', 'OriginSimple', 'Assembly', 'ChromosomeAccession',
    'Cytogenetic', 'ReviewStatus', 'NumberSubmitters', 'Guidelines', 'TestedInGTR', 'OtherIDs',
    'SubmitterCategories', 'VariationID', 'Chromosome', 'Start', 'Stop', 'ReferenceAllele', 'AlternateAllele',
    'PositionVCF', 'ReferenceAlleleVCF', 'AlternateAlleleVCF',
]

COORD_KEYS = ['Chromosome', 'PositionVCF', 'RS# (dbSNP)', 'ReferenceAlleleVCF', 'AlternateAlleleVCF']
VCF_HEADER = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
INVALID_COORD = {'na', '-1'}


def write_records(infile, outfile, assembly):
    seen = set()
    n_out = 0

    with open(infile) as f, open(outfile, 'w') as out:
        headerL = f.readline().strip(' \r\n').split('\t')

        for line in f:
            tokL = line.strip(' \r\n').split('\t')
            tokD = {k: v for k, v in zip(headerL, tokL)}

            if tokD['Assembly'] != assembly:
                continue

            allele_id = tokL[0]
            if allele_id in seen:
                continue
            seen.add(allele_id)

            coordL = [tokD[k] for k in COORD_KEYS]
            if INVALID_COORD & {coordL[0], coordL[1], coordL[3], coordL[4]}:
                continue

            info = ';'.join([f'{k}={tokD[k].replace(";", ",")}' for k in INFO_KEYS])
            out.write('chr' + '\t'.join(coordL + ['.', '.', info]) + '\n')
            n_out += 1

    return n_out


def sort_vcf(infile, outfile, tmpdir=None):
    cmd = ['sort', '-k1,1V', '-k2,2n']
    if tmpdir:
        cmd += ['-T', tmpdir]
    cmd.append(infile)

    with open(outfile, 'w') as out:
        out.write('\t'.join(VCF_HEADER) + '\n')
        out.flush()
        subprocess.run(cmd, stdout=out, env=dict(os.environ, LC_ALL='C'), check=True)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', action='store', required=True)
    parser.add_argument('-o', '--output', dest='output', action='store', default='')
    parser.add_argument('-a', '--assembly', dest='assembly', action='store', default='GRCh38')
    parser.add_argument('-T', '--tmpdir', dest='tmpdir', action='store', default=None)
    args = parser.parse_args()

    infile = args.input
    outfile = args.output
    if outfile == '':
        outfile = infile.replace('.txt', f'.{args.assembly}.sorted.vcf')

    unsorted_file = outfile + '.unsorted'

    print('START!')
    n_out = write_records(infile, unsorted_file, args.assembly)
    sort_vcf(unsorted_file, outfile, args.tmpdir)
    os.remove(unsorted_file)
    print(f'{n_out} records -> {outfile}')
    print('DONE!')

# The population-frequency INFO fields used downstream are left-joined onto this VCF
# before VEP annotation (join key: CHROM, POS, REF, ALT after left-trimming of REF/ALT):
#   - gnomAD v4.1 joint sites : AC/AN/AF/nhomalt, global and per-ancestry (and _XX/_XY on chrX/chrY)
#   - UK Biobank 500k exomes  : AC/AN/AF/nhomalt
#   - All of Us CDR v8        : AC/AN/AF/nhomalt
# Variants absent from a resource keep the sentinel value -1, and the join retains
# chr1-22, chrX and chrY only. The join code is not distributed here because the UK Biobank
# and All of Us allele counts are generated inside their controlled-access environments.
