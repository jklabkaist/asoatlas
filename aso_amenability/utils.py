"""ASO amenability calculation for a variant-transcript pair."""

import datetime
import sys

from jkspliceai.utils import parse_MES_text

from .coding import checkDamageCoding2
from .splicing import checkCanonSplice2, checkMisSplicing2, splicing_alteration

PROBABLY_AMENABLE = 'Probably amenable'
POSSIBLY_AMENABLE = 'Possibly amenable'
UNLIKELY_AMENABLE = 'Unlikely amenable'

intergenicL = ['upstream_gene_variant', 'downstream_gene_variant', 'intergenic_variant']

sai_headerL = ['spliceAI_5k_transcriptID']
for _prefix in ['spliceAI_5k_raw_', 'spliceAI_5k_masked_']:
    for _suf1 in ['relpos_', 'score_', 'MES_']:
        for _suf2 in ['AG', 'AL', 'DG', 'DL']:
            sai_headerL.append(_prefix + _suf1 + _suf2)

sai_cons_headerL = ['spliceAI_molecular_consequence', 'spliceAI_damage2orf_by_splicing']
amen_headerL = ['STATUS_pLOF', 'STATUS_missense', 'STATUS_canon_splice', 'STATUS_mis_splice',
                'ASO_amenability', 'Affected']

amenability_headerL = sai_headerL + sai_cons_headerL + amen_headerL


def log(msg):
    print(f"[{datetime.datetime.now()}] {msg}", file=sys.stderr)


def flatten_sai_text(sai_text, sep='|'):
    # [[(relpos, delta, MES_ref, MES_alt)] * [AG, AL, DG, DL]] for raw and masked scores
    tmpL = sai_text.split(sep)
    tid, tmpL = tmpL[0], tmpL[1:]
    relposL, deltaL, mesL, mrelposL, mdeltaL, mmesL = [[x.split('/') for x in tmpL[i * 4:(i + 1) * 4]]
                                                       for i in range(len(tmpL) // 4)]

    return ([[(r, d) + parse_MES_text(m) for r, d, m in zip(relpos, delta, mes)]
             for relpos, delta, mes in zip(relposL, deltaL, mesL)],
            [[(r, d) + parse_MES_text(m) for r, d, m in zip(relpos, delta, mes)]
             for relpos, delta, mes in zip(mrelposL, mdeltaL, mmesL)],
            tid)


def sai_delta(saiL):
    # keep the absolute change of the SpliceAI score, "ref>alt" -> |alt - ref|
    return [[(r, abs(float(d.split('>')[1]) - float(d.split('>')[0])) if '>' in d else 0, m1, m2)
             for r, d, m1, m2 in xL] for xL in saiL]


def canonical_splice_tag(dCanSpl):
    dCanSpl = [tag for tag in dCanSpl if not tag.startswith('No to little:')]
    tagL = [tag.split(":")[0] for tag in dCanSpl]

    if "Severe" in tagL:
        tag_canSpl = "Severe"
    elif "Moderate" in tagL:
        tag_canSpl = "Moderate"
    else:
        tag_canSpl = "No to little"

    return tag_canSpl, dCanSpl


def classify_amenability(tag_coding, tag_canSpl, missSpl, consequence=''):

    if tag_coding == "Severe" or tag_canSpl == "Severe":
        amen = UNLIKELY_AMENABLE
        affected = True
    elif missSpl.startswith("Gain"):
        affected = True
        if tag_coding == "Moderate" or tag_canSpl == "Moderate":
            amen = POSSIBLY_AMENABLE
        else:
            amen = PROBABLY_AMENABLE
    elif missSpl.startswith("Skipping"):
        amen = POSSIBLY_AMENABLE
        affected = True
    else:
        amen = UNLIKELY_AMENABLE
        affected = False

    if consequence in intergenicL:
        amen = UNLIKELY_AMENABLE

    return amen, affected


def damage2orf(spl_mol_cons):
    return 'ORFdamaged' if any([(cons.split(':')[-1] == 'PTC-gained') or
                                (cons.split(':')[-2:] == ['non3N', 'No_PTC-gained'])
                                for cons in spl_mol_cons.split('/')]) else ''


def coding_only_result(chrom, pos, ref, alt, vepD):
    """Amenability of a variant that carries no usable SpliceAI prediction."""

    dCoding = checkDamageCoding2(chrom, pos, ref, alt, vepD)
    tag_coding = dCoding[0].split(":")[0]

    resD = {k: '' for k in sai_headerL + sai_cons_headerL}
    resD.update({'STATUS_pLOF': dCoding[1],
                 'STATUS_missense': dCoding[2],
                 'STATUS_canon_splice': '',
                 'STATUS_mis_splice': '',
                 'ASO_amenability': UNLIKELY_AMENABLE,
                 'Affected': tag_coding == "Severe"})
    return resD


def amenability_result(chrom, pos, ref, alt, vepD, t, sai_text, masked_saiD):
    """Amenability of a variant-transcript pair.

    t          : transcript record of the refFlat annotation
    sai_text   : '|'-joined SpliceAI/MaxEntScan text of this variant-transcript pair
    masked_saiD: masked SpliceAI result dict of this variant-transcript pair
    """

    dCoding = checkDamageCoding2(chrom, pos, ref, alt, vepD)
    tag_coding = dCoding[0].split(":")[0]

    spl_mol_cons, dmg2orf = '', ''
    try:
        spl_mol_cons = splicing_alteration(chrom, pos, ref, alt, masked_saiD)
        dmg2orf = damage2orf(spl_mol_cons)
    except Exception as e:
        log(f" --- Skipped - molcons: {chrom}-{pos}-{ref}-{alt};{vepD['Feature']} ({e})")

    saiL, msaiL, _tid = flatten_sai_text(sai_text, sep='|')
    saiL, msaiL = sai_delta(saiL), sai_delta(msaiL)

    # damage to the canonical splice sites
    tag_canSpl, dCanSpl = canonical_splice_tag(checkCanonSplice2(chrom, pos, ref, alt, saiL, t))

    # mis-splicing (pseudoexon / exon skipping)
    missSpl = checkMisSplicing2(chrom, pos, ref, alt, saiL, t)

    amen, affected = classify_amenability(tag_coding, tag_canSpl, missSpl, vepD['Consequence'])

    resD = dict(zip(sai_headerL, sai_text.split('|')))
    resD.update({'spliceAI_molecular_consequence': spl_mol_cons,
                 'spliceAI_damage2orf_by_splicing': dmg2orf,
                 'STATUS_pLOF': dCoding[1],
                 'STATUS_missense': dCoding[2],
                 'STATUS_canon_splice': '/'.join(dCanSpl),
                 'STATUS_mis_splice': missSpl,
                 'ASO_amenability': amen,
                 'Affected': affected})
    return resD
