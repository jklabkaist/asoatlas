#!/usr/bin/env python3
"""Merge ClinVar records that describe the same allele (same CHROM, POS, REF, ALT) under
different VariationIDs, so that the allele is counted once."""

import argparse
import sys

SUM_PREFIX = 'subs_'
SUM_FIELDS = {'NumberSubmitters'}
FLAG_FIELDS = {'Low_penetrance_variant', 'Conflicting_but_likely_non_pathogenic',
               'Predominant_VUS', 'Ambiguous'}
OUTLIER_FIELDS = ('over3x99tile_AoU', 'over3x99tile_UKB', 'over3x99tile_gnomAD')
CLINVAR_FIELDS = {'Type', 'Name', 'GeneID', 'GeneSymbol', 'HGNC_ID', 'ClinicalSignificance',
                  'ClinSigSimple', 'LastEvaluated', 'RS# (dbSNP)', 'nsv/esv (dbVar)',
                  'RCVaccession', 'PhenotypeIDS', 'PhenotypeList', 'Origin', 'OriginSimple',
                  'ReviewStatus', 'Guidelines', 'TestedInGTR', 'OtherIDs', 'SubmitterCategories',
                  'VariationID', 'AlleleID', 'NumberSubmitters'}

report_h = ['CHROM', 'POS', 'REF', 'ALT', 'n_records', 'VariationID_merged', 'VariationID_members',
            'Name', 'Name_merged', 'subsPLP_members', 'subsPLP_merged', 'subsBLB_members',
            'subsBLB_merged', 'subsVUS_members', 'subsVUS_merged', 'Ambiguous_members',
            'Ambiguous_merged']


def parse(line):
    t = line.rstrip('\n').split('\t')
    info, csq = t[7].split(';CSQ=')
    kv = [x.split('=', 1) for x in info.split(';')]
    return t, [k for k, _ in kv], dict(kv), csq


def union_tokens(vals):
    out = []
    for v in vals:
        for tok in v.split('|'):
            if tok not in out:
                out.append(tok)
    return '|'.join(out)


def primary_order(rec):
    _, _, d, _ = rec
    return (0 if d['Name'].startswith('NM_') else 1, -int(d['subs_P/LP']), int(d['VariationID']))


def merge(recs, report):
    recs = sorted(recs, key=primary_order)
    t0, keys0, d0, csq0 = recs[0]
    for t, keys, d, csq in recs[1:]:
        if csq != csq0:
            sys.exit(f'ABORT: CSQ differs within {t[0]}:{t[1]} {t[3]}>{t[4]}')
        if keys != keys0:
            sys.exit(f'ABORT: INFO keys differ within {t[0]}:{t[1]} {t[3]}>{t[4]}')
        for k in keys0:
            if d[k] != d0[k] and k not in CLINVAR_FIELDS and not k.startswith(SUM_PREFIX) \
                    and k not in FLAG_FIELDS and k not in OUTLIER_FIELDS:
                sys.exit(f'ABORT: non-ClinVar field {k} differs within {t[0]}:{t[1]} {t[3]}>{t[4]}')

    D = [r[2] for r in recs]
    m = dict(d0)
    for k in keys0:
        vals = [d[k] for d in D]
        if k.startswith(SUM_PREFIX) or k in SUM_FIELDS:
            m[k] = str(sum(int(v) for v in vals))
        elif k == 'Name':
            m[k] = vals[0]
        elif k == 'ClinicalSignificance' and set(vals) == {'Pathogenic', 'Likely pathogenic'}:
            m[k] = 'Pathogenic/Likely pathogenic'
        elif k in CLINVAR_FIELDS and len(set(vals)) > 1:
            m[k] = union_tokens(vals)

    # step 06 flags, from the summed submission counts
    plp, blb = int(m['subs_P/LP']), int(m['subs_B/LB'])
    vus, lowp = int(m['subs_uncertain']), int(m['subs_low_penetrance'])
    f_lowp = lowp > 0
    f_conf = plp <= blb
    f_pvus = blb > 0 and vus + blb > plp
    m['Low_penetrance_variant'] = str(f_lowp)
    m['Conflicting_but_likely_non_pathogenic'] = str(f_conf)
    m['Predominant_VUS'] = str(f_pvus)

    # step 07: the outlier flags depend on the allele frequency only
    if f_lowp or f_conf or f_pvus:
        for h in OUTLIER_FIELDS:
            m[h] = '-'
        amb = True
    else:
        outl = []
        for h in OUTLIER_FIELDS:
            vals = {d[h] for d in D}
            if '-' in vals:
                sys.exit(f'ABORT: outlier flag {h} was not evaluated for a member of '
                         f'{t0[0]}:{t0[1]}; re-run step 07 for this allele')
            if len(vals) > 1:
                sys.exit(f'ABORT: outlier flag {h} differs within {t0[0]}:{t0[1]}')
            m[h] = vals.pop()
            outl.append(m[h] == 'True')
        amb = any(outl) and (blb >= 1 or vus >= 2)
    m['Ambiguous'] = str(amb)
    m['Name_merged'] = union_tokens([d['Name'] for d in D])
    m['ClinVar_records_merged'] = str(len(recs))

    ids = [x for r in recs for x in r[0][2].split(';') if x not in ('-1', '.', '')]
    t = list(t0)
    t[2] = ';'.join(dict.fromkeys(ids)) or t0[2]
    t[7] = ';'.join(f'{k}={m[k]}' for k in keys0 + ['Name_merged', 'ClinVar_records_merged']) + ';CSQ=' + csq0

    report.write('\t'.join([t0[0], t0[1], t0[3], t0[4], str(len(recs)),
                            m['VariationID'], '|'.join(d['VariationID'] for d in D),
                            m['Name'], m['Name_merged'],
                            '|'.join(d['subs_P/LP'] for d in D), m['subs_P/LP'],
                            '|'.join(d['subs_B/LB'] for d in D), m['subs_B/LB'],
                            '|'.join(d['subs_uncertain'] for d in D), m['subs_uncertain'],
                            '|'.join(d['Ambiguous'] for d in D), m['Ambiguous']]) + '\n')
    return '\t'.join(t) + '\n'


def merge_duplicate_alleles(infile, outfile, reportfile):
    n_in = n_out = n_merged = 0
    block, block_key = [], None
    with open(infile) as f, open(outfile, 'w') as o, open(reportfile, 'w') as r:
        r.write('\t'.join(report_h) + '\n')

        def flush():
            nonlocal n_out, n_merged
            groups = {}
            for line in block:
                rec = parse(line)
                groups.setdefault((rec[0][3], rec[0][4]), []).append((line, rec))
            for g in groups.values():
                if len(g) == 1:
                    o.write(g[0][0])
                else:
                    o.write(merge([x[1] for x in g], r))
                    n_merged += len(g)
                n_out += 1

        for line in f:
            if line.startswith('##'):
                o.write(line)
                continue
            if line.startswith('#CHROM'):
                o.write('##INFO=<ID=Name_merged,Number=.,Type=String,Description="ClinVar names of all '
                        'records merged into this allele">\n')
                o.write('##INFO=<ID=ClinVar_records_merged,Number=1,Type=Integer,Description="Number of '
                        'ClinVar records (VariationIDs) merged into this allele">\n')
                o.write(line)
                continue
            n_in += 1
            key = tuple(line.split('\t', 2)[:2])
            if key != block_key and block:
                flush()
                block = []
            block_key = key
            block.append(line)
        if block:
            flush()
    print(f'records in : {n_in:,}')
    print(f'records out: {n_out:,} ({n_merged} records merged into {n_merged - (n_in - n_out)} alleles)')


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', action="store", required=True,
                        help="VCF from 07_filter_af_outlier.py, sorted by position")
    parser.add_argument('-o', '--output', dest='output', action="store", default="")
    parser.add_argument('-r', '--report', dest='report', action="store", default="",
                        help="TSV listing every merged allele")
    args = parser.parse_args()

    infile = args.input
    outfile = args.output if args.output else infile.replace('.vcf', '.dedup.vcf')
    reportfile = args.report if args.report else outfile.replace('.vcf', '.merged_alleles.tsv')

    print("START!")
    merge_duplicate_alleles(infile, outfile, reportfile)
    print("DONE!")
