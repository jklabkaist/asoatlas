"""Coding-effect evaluation for ASO amenability (VEP consequence / REVEL)."""

import collections
import re

coding_dmg_vartype = collections.defaultdict(lambda:"pass")
coding_dmg_vartype.update({"stop_gained" : "Severe",
                           "frameshift_variant" : "Severe",
                           "start_lost" : "Severe"})

intragenic_variant_typeL = ['stop_gained','frameshift_variant','start_lost','5_prime_UTR_variant','3_prime_UTR_variant',
                            'splice_acceptor_variant','splice_donor_variant','intron_variant','splice_region_variant',
                            'inframe_insertion','inframe_deletion','missense_variant','protein_altering_variant',
                            'start_retained_variant','stop_retained_variant','synonymous_variant','coding_sequence_variant',
                            'stop_lost','splice_donor_5th_base_variant','splice_donor_region_variant','splice_polypyrimidine_tract_variant']

sift_polyphen_dmg_vartype = {"deleterious": "bad",
                             "deleterious_low_confidence":"bad",
                             "probably_damaging":"bad",
                             "possibly_damaging":"bad",
                             "tolerated":"pass","benign":"pass","-":"pass"}



def isIntergenic(chrom, pos, ref, alt, vepD):
    return sum([vtype in vepD['Consequence'] for vtype in intragenic_variant_typeL]) == 0

def isCodingLOF(chrom, pos, ref, alt, vepD):
    consL = re.split(',|&',vepD['Consequence'])
    sevL  = [cons for cons in consL if coding_dmg_vartype[cons]=="Severe"]
    if len(sevL)>0:
        return f"Severe: {'/'.join(sevL)}"
    else:
        return f"No to little: pLOF variant"

def isMissensePatho(chrom, pos, ref, alt, vepD): #SIFT/PP2
    badSift = any([v=='bad' for k,v in sift_polyphen_dmg_vartype.items() if k==vepD['SIFT'].split('(')[0]])
    badPP   = any([v=='bad' for k,v in sift_polyphen_dmg_vartype.items() if k==vepD['PolyPhen'].split('(')[0]])

    if badSift and badPP:
        return "Severe: SIFT_PP2"
    elif badSift or badPP:
        return "Moderate: SIFT_PP2"
    else:
        return "No to little: SIFT_PP2"

def isMissensePatho2(chrom, pos, ref, alt, vepD, th=0.5): # REVEL
    return "No to little: REVEL" if vepD['REVEL']=='' or float(vepD['REVEL']) <= th else "Severe: REVEL"

def codingseq_indel(chrom, pos, ref, alt, vepD):
    return (len(ref)!=len(alt) and vepD['EXON']!='' and vepD['CDS_position']!='') or ('-' in vepD['Protein_position'])

def checkDamageCoding(chrom, pos, ref, alt, vepD):

    if isIntergenic(chrom, pos, ref, alt, vepD):
        return "intergenic"
    
    cLOF = isCodingLOF(chrom, pos, ref, alt, vepD)
    pMiss = isMissensePatho(chrom,pos, ref, alt, vepD)

    if cLOF or pMiss=="Severe":
        return "Severe"
    elif pMiss=="Moderate":
        return "Moderate"
    else:
        return "No to little"

def checkDamageCoding2(chrom, pos, ref, alt, vepD):

    if isIntergenic(chrom, pos, ref, alt, vepD):
        return ("intergenic", '','')

    cLOF = isCodingLOF(chrom, pos, ref, alt, vepD) 
    pMiss = isMissensePatho2(chrom, pos, ref, alt, vepD)

    if codingseq_indel(chrom, pos, ref, alt, vepD):
        pMiss += '; Severe: coding_sequence_indel'

    if cLOF.split(':')[0]=="Severe" or "Severe" in pMiss:
        return ("Severe", cLOF, pMiss)
    else:
        return ("No to little", cLOF, pMiss)
