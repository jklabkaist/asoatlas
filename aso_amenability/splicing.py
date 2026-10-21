"""Splicing-effect evaluation for ASO amenability (SpliceAI / MaxEntScan / LaBranchoR)."""

import copy
import itertools

import tabix

import jkspliceai
from jklib.genome import locus, loadBlatOutput
from jklib.utils import margin, overlap

from .config import LABRANCHOR_PATH, REFFLAT_MANE_PATH

h = loadBlatOutput(REFFLAT_MANE_PATH)
labr = tabix.open(LABRANCHOR_PATH)

codons_short = {"UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V",
                "UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V",
                "UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
                "UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V",
                "UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A",
                "UCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
                "UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
                "UCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
                "UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D",
                "UAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
                "UAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
                "UAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
                "UGU" : "C", "CGU" : "R", "AGU" : "S", "GGU" : "G",
                "UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
                "UGA" : "*", "CGA" : "R", "AGA" : "R", "GGA" : "G",
                "UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G"}

rc = lambda seq: seq.translate(seq.maketrans('ACGT','TGCA'))[::-1]
def translation(seq, codonD=codons_short):
    res=''
    for i in range(len(seq)//3):
        try:
            res+=codonD[seq[i*3:(i*3)+3].replace('T','U')]
        except:
            continue
    return res

def getRegionType_givenTranscript(loc_,t):

    def region_idx(i, n_region, strand):
        if strand=='+':
            return i+1
        else:
            return n_region-i

    loc = copy.copy(loc_)
    loc.chrSta = loc.chrSta -1
    strand_ = loc.strand if loc.strand else '+'

    regions=[]

    if t['strand']==strand_:
        sense='sense'
    else:
        sense='antisense'

    tstrand = t['strand']
    t_exn_n = len(t['exnList'])
    t_int_n = len(t['intron'])

    flag_default=''
    if t['cdsSta'] == t['cdsEnd']:
        flag_default='lnc_'

        for regionName in ['exnList', 'intron']:
            flag = flag_default
            for i,(s,e) in enumerate(t[regionName]):
                if overlap(loc,locus('%s:%s-%s%s'%(loc.chrom,s,e,t['strand'])))==0:
                    continue

                if regionName == 'exnList':
                    flag += f'exon{region_idx(i, t_exn_n, tstrand)}'
                else:
                    flag += f'intron{region_idx(i, t_int_n, tstrand)}'
                
                marg = margin((loc.chrSta, loc.chrEnd), (s,e))
                regions.append((t['transName'],t['transID'], sense, str([(flag, marg)])))
                flag= flag_default

    else:
        for regionName in ['cdsList','utr5','utr3','intron']:
            flag = flag_default
            for s,e in t[regionName]:
                if overlap(loc, locus('%s:%s-%s%s' % (loc.chrom, s,e,t['strand']))) == 0:
                    continue
                if regionName == 'cdsList':
                    if loc.chrEnd-loc.chrSta==0:
                        offset=0
                        for s_,e_ in t['cdsList']:
                            offset += min(e_, loc.chrSta)-min(s_,loc.chrSta)

                        frame = (offset-1) % 3
                        if t['strand']=='-':
                            frame = 2-frame
                        flag += 'cds_%d_%d' % (offset-1, frame)
                    else:
                        flag += 'cds_m'
                else:
                    flag += regionName

                marg = margin((loc.chrSta, loc.chrEnd), (s,e))
                regions.append((t['transName'],t['transID'], sense, str([(flag, marg)])))
                flag = flag_default

    regions = list(set(regions))

    if len(regions)<2:
        return list(map(lambda x: x[:3] + (eval(x[-1]),), regions ))

    regions.sort()

    regions_ = []
    tmp = regions[0]

    for ann in regions[1:]:
        if tmp[:3]==ann[:3]:
            tmp3=tmp[3]+ '+' + ann[3]
            tmp = tmp[:3] + (tmp3,)
        else:
            regions_.append(tmp)
            tmp = ann

    regions_.append(tmp)
    regions_ = list(map(lambda x: x[:3] + (eval(x[-1]),), regions_))

    return regions_


def iscanonicalsplicesite(loc, t):

    if isinstance(loc, str) and 'ins' in loc:
        return 'Not canonical'

    for region in getRegionType_givenTranscript(loc,t):
        for item in region[3]:
            if 'intron' in item:
                sta, end = item[1]
                if (t['strand'] =='+' and end==-1) or (t['strand']=='-' and sta==-1):
                    return 'acceptor_canonical'
                if (t['strand']=='+' and sta==-1) or (t['strand']=='-' and end==-1):
                    return 'donor_canonical'

    return 'Not canonical'

def labranchor_check(loc,ref,alt,labr):
    def between(t,a,b):
        if a>b:
            a,b = b,a
        return a<t<b

    tabix_query=f'{loc.chrom}:{loc.chrSta-50}-{loc.chrSta+50}+'
    labr_values = [x for x in labr.querys(tabix_query)]

    if not labr_values or len(labr_values)==0:
        return False

    else:
        for labr_value in labr_values:
            #labranchor data: chrom pos(0-based) pos+1 canonical_acceptor(0-based) strand
            #branchpoint = f'{labr_value[0]}:{labr_value[2]}-{labr_value[2]}+'
            #canonical_acc = f'{labr_value[0]}:{int(labr_value[3])+1}-{int(labr_value[3])+1}+'

            #The variant is <=3 nt away from the LaBranchoR-predicted branchpoint
            if overlap(locus(f'{loc.chrom}:{loc.chrSta}-{loc.chrSta+len(ref)-1}+'),locus(f'{labr_value[0]}:{int(labr_value[2])-3}-{int(labr_value[2])+3}+')):
                return True
            #the distance between the LaBranchoR-predicted branchpoint and the canonical splice site is changed by >3 nt by the variant
            if between(loc.chrSta, int(labr_value[2]), int(labr_value[3])+1) and abs(len(alt)-len(ref))>3:
                return True
    return False

def damage_splicing_condition(loss_score, mes_score, labr_status, region_check, delta_th=0.1):

    if (region_check == 'donor_canonical' or region_check=='acceptor_canonical'):
        if loss_score >= delta_th and (mes_score[3] < mes_score[1]) and (mes_score[3] < 2 or mes_score[3] < (0.3 * mes_score[1])):
            return f'Severe: {region_check} loss [spliceAI: {loss_score:.3f} mes: {mes_score[1]} -> {mes_score[3]} branchpoint_dmg: {labr_status}]'
        elif (loss_score >= delta_th and (0.3 * mes_score[1] <= mes_score[3] < mes_score[1] and mes_score[3]>=2)) or labr_status:
            return f'Moderate: {region_check} loss [spliceAI: {loss_score:.3f} mes: {mes_score[1]} -> {mes_score[3]} branchpoint_dmg: {labr_status}]'
        else:
            return 'No to little: damage to canonical splicing'
    else:
        return 'No to little: damage to canoncial splicing' #?


def checkCanonSplice2(chrom, pos, ref, alt, saiL, t, delta_th=0.1):

    cs_loss_status=[]

    loc = locus(f'{chrom}:{pos}-{pos}+')
    labr_status = labranchor_check(loc, ref, alt, labr=labr)

    # AL in canonical Acceptor:
    alL = [(r, float(d), mr, ma) for r,d,mr, ma in saiL[1] if d != "." and float(d)>=delta_th]
    dlL = [(r, float(d), mr, ma) for r,d,mr, ma in saiL[3] if d != "." and float(d)>=delta_th]

    
    for r,d,mr,ma in alL:

        al_pos = int(pos) + int(r)
        al_loc = locus(f'{chrom}:{al_pos}-{al_pos}+')

        al_region_check = ''
        #print(getRegionType_givenTranscript(al_loc,t))

        for al_region in getRegionType_givenTranscript(al_loc,t):
            for item in al_region[3]:
                if 'intron' in item:
                    sta, end = item[1]
                    if t['strand']=='+' and end==-1:
                        al_region_check = 'acceptor_canonical'
                    elif t['strand']=='-' and sta==-1:
                        al_region_check = 'acceptor_canonical'
        
        mr = -999 if mr=='.' or 'N' in mr else float(mr)
        ma = -999 if ma=='.' or 'N' in ma else float(ma)
        mesL = [None, mr, None, ma, None]

        #print(r, d, mesL, labr_status, al_region_check)
        cs_loss_status.append(damage_splicing_condition(d, mesL, labr_status, al_region_check, delta_th=delta_th))

    for r,d,mr,ma in dlL:

        dl_pos = int(pos) + int(r)
        dl_loc = locus(f'{chrom}:{dl_pos}-{dl_pos}+')

        dl_region_check = ''
        #print(getRegionType_givenTranscript(dl_loc,t))

        for dl_region in getRegionType_givenTranscript(dl_loc,t):
            for item in dl_region[3]:
                if 'intron' in item:
                    sta, end = item[1]
                    if t['strand'] == '+' and sta==-1:
                        dl_region_check = 'donor_canonical'
                    elif t['strand'] == '-' and end==-1:
                        dl_region_check = 'donor_canonical'

        mr = -999 if mr=='.' or 'N' in mr else float(mr)
        ma = -999 if ma=='.' or 'N' in ma else float(ma)
        mesL = [None, mr, None, ma, None]

        #print(r, d, mesL, labr_status, dl_region_check)
        #print(damage_splicing_condition(d, mesL, labr_status, dl_region_check, delta_th=delta_th))

        cs_loss_status.append(damage_splicing_condition(d, mesL, labr_status, dl_region_check, delta_th=delta_th))
        
    return cs_loss_status


def checkMisSplicing2(chrom, pos, ref, alt, saiL, t, masked=False, delta_th_gain=0.1, delta_th_skip=0.05):

    #print(saiL)
    agL_, alL_, dgL_, dlL_ = [[ ( r if 'ins' in r else locus(f'{chrom}:{int(pos)+int(r)}-{int(pos)+int(r)}+') , float(d), ma ) for r,d,_,ma in sai if d != "." and float(d)>=min(delta_th_gain,delta_th_skip)] for sai in saiL ]
    #print([agL_, alL_, dgL_, dlL_])
    agL, alL, dgL, dlL     = [[ ( l[0].toString() if isinstance(l[0], locus) else l[0], l[1], iscanonicalsplicesite(l[0],t), float(l[2]) if ('.' != l[2] and 'N' not in l[2]) else -999 ) for l in lL ] for lL in [agL_, alL_, dgL_, dlL_]]

    # gain: AG, DG non-canonical AND score >=delta_th_gain
    status_tmpL = []
    for ag in [ag for ag in agL if ag[2]!='acceptor_canonical' and ag[3]>=2 and ag[1]>=delta_th_gain]:
        status_tmpL.append(f'AG [{ag[0]} {ag[1]:.3f}]')
    for dg in [dg for dg in dgL if dg[2]!='donor_canonical' and dg[3]>=2 and dg[1]>=delta_th_gain]:
        status_tmpL.append(f'DG [{dg[0]} {dg[1]:.3f}]')

    if len(status_tmpL) != 0:
        return 'Gain of mis-splicing:' + '/'.join(status_tmpL)
    
    # skipping/retention: AL,DL canonical AND score >=delta_th_skip AND AG,DG at non-canocnial
    #print(agL, alL, dgL, dlL)
    status_tmpL = []
    for al in [al for al in alL if al[2]=='acceptor_canonical' and al[1]>=delta_th_skip]:
        status_tmpL.append(f'AL [{al[0]} {al[1]:.3f}]')
    for dl in [dl for dl in dlL if dl[2]=='donor_canonical' and dl[1]>=delta_th_skip]:
        status_tmpL.append(f'DL [{dl[0]} {dl[1]:.3f}]')

    if len(status_tmpL) != 0:
        return 'Skipping_retention:' + '/'.join(status_tmpL)

    return 'No miss-splicing'

def isORFchanging(chrom,sta,end,tD):
    #sta, end : abspos (pos, offset)
    return overlap(locus(f"{chrom}:{sta[0]}-{end[0]}+"),locus(f"{tD['chrom']}:{tD['cdsSta']+1}-{tD['cdsEnd']}+")) #0-based > 1-based

#exnL: 1-based [ ( (sta,sta_offset), (end,end_ofset) )]
# ins_n: nth base of insertion or indel sequence (1-based)
#        converted into (pos, n-1) -> offset

def insert_merge_exon(sta, end, exnL):
    # sta, end: (pos, offset)
    # exnL: [ ((sta_pos, sta_offset), (end_pos, end_offset)) ]
    resL = []
    for s,e in exnL:
        if sta > e:
            resL.append((s, e))
        elif s > end:
            resL.append((sta, end))
            sta, end = s, e
        else:
            sta, end = min(sta, s), max(end, e)
    resL.append((sta, end))

    if len(resL)>0:
        mergedL = [resL[0]]
        for now in resL[1:]:
            ns, ne = now
            ps, pe = mergedL[-1]
            if (ns[0]==pe[0]+1 and ns[1]==pe[1]==0) or (pe[0]==sta[0] and ns[1]==pe[1]+1): # merge retention
                mergedL[-1] = (ps,ne)
            else:
                mergedL.append(now)   
        resL = mergedL
    return resL

def truncate_exon(sta, end, exnL):  #1-based
    # sta, end: (pos, offset)
    # exnL: [ ((sta_pos, sta_offset), (end_pos, end_offset)) ]
    resL = []
    for s, e in exnL:
        if e < sta or end < s: # No overlap
            resL.append((s, e))
        else: # Overlap, truncate the existing range
            if s < sta:
                if sta[1]>0:
                    resL.append((s, (sta[0],sta[1]-1)))
                else:
                    resL.append((s, (sta[0]-1,sta[1])))
            if end < e:
                if end[1]>0:
                    resL.append(((end[0],end[1]+1), e))
                else:
                    resL.append(((end[0]+1,end[1]), e))
    return resL

def get_mRNAsequence(pos, ref, alt, exnL, tD, cdsOnly=False): 
    #exnL: 1-based; return spliced mRNA sequence for exonL
    #return: list of exon sequences

    if jkspliceai.sai.MANEseqD[tD['transID']]['dnaseq'] == '':
        jkspliceai.sai.MANEseqD[tD['transID']]['dnaseq'] = locus(f"{tD['chrom']}:{tD['txnSta']+1}-{tD['txnEnd']}+").twoBitFrag() #0-based > 1-based
    
    dnaseq=jkspliceai.sai.MANEseqD[tD['transID']]['dnaseq']
    dnaseq_alt=dnaseq
    offset_shift=0
    chrSta=tD['txnSta']+1 #1-based

    if tD['txnSta']<pos<=tD['txnEnd']:
        mut_offset= pos - tD['txnSta'] #check 0/1 based
        dnaseq_alt = dnaseq[:mut_offset] + alt + dnaseq[mut_offset+len(ref):]
        # make offsetL [(s,e)], regarding shifting by indel variants
        offset_shift = len(alt) - len(ref)

    offsetL = []
    for s,e in exnL:
        s_ = s[0]+s[1]+offset_shift if s[0] > pos else s[0]+s[1]
        e_ = e[0]+e[1]+offset_shift if e[0] > pos else e[0]+e[1]
        offsetL.append((s_-chrSta, e_-chrSta)) # equal or not
    
    if cdsOnly:
        offsetL = []
        cdssta, cdsend = (tD['cdsSta']+1,0), (tD['cdsEnd'],0)
        for s,e in [ (max(s1,cdssta),min(e1,cdsend)) for s1,e1 in exnL if max(s1,cdssta)!=min(e1,cdsend) ]:
            s_ = s[0]+s[1]+offset_shift if s[0] > pos else s[0]+s[1]
            e_ = e[0]+e[1]+offset_shift if e[0] > pos else e[0]+e[1]
            offsetL.append((s_-chrSta, e_-chrSta)) # equal or not

    seqL = [dnaseq_alt[s:e+1].upper() for s,e in offsetL]
    
    if tD['strand']=='-':
        seqL=[rc(seq) for seq in seqL][::-1]

    return seqL

def isNMDsplicing(pos, ref, alt, aa_seq, changed_exnL, tD):
    # 1) stop codon falls in the last exon
    # 2) ~50bp upstream of second to the last exon
    #       if pseudo exon included bw last two exons -> 
    # 3) first 100 cds bases -> escape ! skip this option
    # 4) intronless gene -> escape

    exnL = changed_exnL if tD['strand']=='+' else changed_exnL[::-1]
    if len(tD['exnList'])<=3: #short gene -> no escaping rule
        return True
    abspos_last2up50 = max(exnL[-2][0],  shift_n_left_abspos(exnL[-2][1],pos, ref, alt, n=50)) if tD['strand']=='+' else min(exnL[-2][1],  shift_n_right_abspos(exnL[-2][0], pos, ref, alt, n=50))
    cdssta, cdsend = (tD['cdsSta']+1,0), (tD['cdsEnd'],0)
    changed_cdsL = [(max(s,cdssta),min(e,cdsend)) for s,e in changed_exnL if max(s,cdssta)!=min(e,cdsend)]

    def cpos2abspos(cpos, changed_cdsL, pos, ref, alt, tD):
        tmp_cpos_sum = 0
        for s,e in changed_cdsL: # 1-based
            se_len=len_abspos(s,e,pos, ref, alt)
            if tmp_cpos_sum + se_len >= cpos:
                return shift_n_right_abspos(s, pos, ref, alt, n=(cpos - tmp_cpos_sum - 1)) if tD['strand']=='+' else shift_n_left_abspos(e, pos, ref, alt, n=(cpos - tmp_cpos_sum - 1))
            tmp_cpos_sum += se_len
        raise Exception(f"No matched HGVS c position: {cpos}, {tD['transID']}")
    
    if len(tD['exnList'])==1: #4)
        return False

    try:
        first_ter_ppos = aa_seq.index('*') + 1 # aa_offset of PTC
        first_ter_cpos = (3*first_ter_ppos - 2, 3*first_ter_ppos) # (,) inclusive
        first_ter_abspos = (cpos2abspos(first_ter_cpos[0],changed_cdsL,pos, ref, alt,tD), cpos2abspos(first_ter_cpos[1],changed_cdsL,pos, ref, alt,tD))
    except Exception as e:
        raise e

    if first_ter_ppos == len(aa_seq): # no PTC
        return False
    elif (tD['strand']=='+' and first_ter_abspos[1] >= abspos_last2up50) or (tD['strand']=='-' and first_ter_abspos[0] <= abspos_last2up50): #1, 2)
        return False
    else:
        return True

def relpos2abspos(pos, relpos):
    # relpos (str) or relpos_in_INS (str; 'ins_n')
    relpos = str(relpos)
    return (int(pos), int(relpos.replace('ins_',''))) if 'ins_' in relpos else (int(pos) + int(relpos), 0)

def abspos2relpos(pos, tup):
    pos = int(pos) #casting
    if tup[1]>0:
        if tup[0]==pos:
            return f'ins_{tup[1]}'
        else:
            print(pos, tup)
            raise Exception
    else:
        return f'{tup[0]-pos}'

def shift_right_abspos(tup, pos, ref, alt): # +1
    return shift_n_right_abspos(tup, pos, ref, alt, n=1)
    """
    if tup[1]>0:
        return (tup[0]+1, 0) if tup[1]+1 >= len(alt) else (tup[0],tup[1]+1) # overflow
    else:
        return (tup[0]+1, 0)
    """

def shift_n_right_abspos(tup, pos, ref, alt, n=1):
    if tup[0]==pos:
        return (tup[0] + tup[1]+ n - (len(alt)-len(ref)) , 0) if tup[1]+n >= len(alt) else (tup[0],tup[1]+n)
    elif tup[0]>pos:
        return (tup[0]+n, 0)
    else: #tup[0]<pos
        if tup[0]+n<pos:
            return (tup[0]+n,0)
        else: # tup[0]+n>=pos
            x = tup[0]+n-pos
            return (pos, x) if x<len(alt) else (tup[0]+n-len(alt)+len(ref),0)

def shift_left_abspos(tup, pos, ref, alt):    # -1
    return shift_n_left_abspos(tup, pos, ref, alt, n=1)
    """
    if tup[1]>=1:
        return (tup[0],tup[1]-1)
    else:
        if tup[0]-1==pos and len(ref)<len(alt):
            return (tup[0]-1,len(alt))
        else:
            return (tup[0]-1,0)
    """

def shift_n_left_abspos(tup, pos, ref ,alt, n=1):
    if tup[0]==pos:
        return (tup[0]+tup[1]-n,0) if tup[1]-n<0 else (tup[0],tup[1]-n)
    elif tup[0]<pos:
        return (tup[0]-n,0)
    elif tup[0]>pos:
        if tup[0]-n >= pos + len(ref):
            return (tup[0]-n,0)
        else:
            x = tup[0] - pos + len(alt) - len(ref) - n 
            return (pos, x) if x>=0 else (x + pos, 0)

def predict_splMolCons(chrom, pos, ref, alt, tD, relpos_sta, relpos_end, splicing_type):
    
    # intronic splicing gain - pseudo exon / retention
    # exonic splicing gain/loss - entire exon skipping / truncation
    # exitron cases (intra-exonic splicing gain)
    # stop codon gain in last 2 exons -> ?
    # label if reading frame is shifed

    # relpos_sta: str - relpos
    # relpos_end: str - relpos
    # splicing_type: inclusion
    #                exclusion
    #                retention  -> s, e inclusion
    #                exitron    -> s+1, e-1 exclusion
    def find_mismatch_idx(aa_ref,aa_alt):
        for i, (a1, a2) in enumerate(zip(aa_ref,aa_alt)):
            if a1!=a2:
                return i
        return min(len(aa_ref),len(aa_alt))


    assert splicing_type in ['inclusion','exclusion','exitron']

    sta, end = relpos2abspos(pos, relpos_sta), relpos2abspos(pos, relpos_end)
    if sta>end:
        sta, end = end, sta

    if splicing_type=='exitron': # recursive call for exclusion
        relpos_sta = abspos2relpos(pos, shift_right_abspos(sta, pos, ref, alt)) # sta +1
        relpos_end = abspos2relpos(pos, shift_left_abspos(end, pos, ref, alt)) # end -1
        return predict_splMolCons(chrom, pos, ref, alt, tD, relpos_sta, relpos_end, 'exclusion')

    if not isORFchanging(chrom,sta,end,tD):
        return 'No_damage_to_ORF'

    #[ ( (sta_pos, sta_offset),(end_pos, end_offset) ) ]
    original_cdsL = [((s+1,0),(e,0)) for s,e in tD['cdsList']] #1-based
    original_exnL = [((s+1,0),(e,0)) for s,e in tD['exnList']] #1-based

    if splicing_type=='exclusion':# exclusion
        changed_cdsL = truncate_exon(sta, end, original_cdsL)
        changed_exnL = truncate_exon(sta, end, original_exnL)
    else: # inclusion
        changed_cdsL = insert_merge_exon(sta, end, original_cdsL)
        changed_exnL = insert_merge_exon(sta, end, original_exnL)

    #mRNAL_ref_cutby_exon = get_mRNAsequence(0, '', '', original_cdsL, tD)
    #mRNAL_alt_cutby_exon = get_mRNAsequence(pos, ref, alt, changed_cdsL, tD)
    cdsL_alt_cutby_exon = get_mRNAsequence(pos, ref, alt, changed_exnL, tD, cdsOnly=True)
    
    #aa_ref = translation(''.join(mRNAL_ref_cutby_exon),codonD=codons_short)
    aa_alt = translation(''.join(cdsL_alt_cutby_exon),codonD=codons_short)
    #aa_change_from = (find_mismatch_idx(aa_ref,aa_alt)+1) / len(aa_ref)

    #diff_len=abs(sum([len(x) for x in mRNAL_alt_cutby_exon]) - tD['cdsLen'])
    diff_len = len_abspos(sta, end, pos, ref, alt)
    diff_3n = '3N' if diff_len%3 ==0 else 'non3N'
    diff_cds = f"{diff_len}:{diff_3n}"
    #diff_cds = f"{aa_change_from:.2f}:{diff_len}:{diff_3n}"

    # skip non-NMD
    if '*' not in aa_alt[:-1]:
        return f'{diff_cds}:No_PTC-gained'

    # skip NMD escaping variant
    if not isNMDsplicing(pos, ref, alt, aa_alt, changed_exnL, tD): 
        return f'{diff_cds}:PTC-gained_but_NMD_escaping'

    # NMD variant
    return f'{diff_cds}:PTC-gained'

def both_in_an_exon(acc, don, tD): # acc, don: parsed_relpos (pos, offset)
    if tD['strand']=='+':
        return any([ (s+1,0) <= acc < don <= (e,0) for s,e in tD['exnList']])
    else:
        return any([ (s+1,0) <= don < acc <= (e,0) for s,e in tD['exnList'] ])


def isExitron(acc, don, tD): # acc, don: parsed_relpos (pos, offset)
    # both in an exon but not canonical splice site
    # in order of  canAcc > newDon > newAcc > canDon
    if tD['strand']=='+':
        return any([ (s+1,0) < don < acc < (e,0) for s,e in tD['exnList'] ])
    else:
        return any([ (s+1,0) < acc < don < (e,0) for s,e in tD['exnList'] ])

def len_abspos(tup1, tup2, pos, ref, alt):
    if tup1>tup2:
        tup1, tup2 = tup2, tup1
    
    if tup1[0]<pos<tup2[0]:
        return tup2[0]-tup1[0]+1 + len(alt)-len(ref)
    if tup2[0]<pos or tup1[0]>pos:
        return tup2[0]-tup1[0]+1
    if tup1[0]==pos:
        if tup1[1]>0:
            return tup2[0]-tup1[0]+1 + len(alt)-len(ref) - tup1[1]
        else:
            return tup2[0]+tup2[1] - tup1[0] + 1
    if tup2[0]==pos:
        return tup2[0]+tup2[1] - tup1[0] + 1
    if tup1[0]==tup2[0]:
        return tup2[1]-tup1[1]+1

    raise Exception(f"not proper case:{tup1}, {tup2}, {pos}, {ref}, {alt}")

def both_in_an_intron(acc, don, tD):
    if tD['strand']=='+':
        return any([ (s+1,0) <= acc < don <= (e,0) for s,e in tD['intron'] ])
    else:
        return any([ (s+1,0) <= don < acc <= (e,0) for s,e in tD['intron'] ])

def current_region(pos,seL):
    for s,e in seL:
        if s<pos[0]<=e:
            return (s+1,e) #1-based
    return (0,0)


def splicing_alteration(chrom, pos, ref, alt, sD): #sD: masked splicing dict
    # pick paired start, end site of splicing alteration
    #   pseudoexon inclusion
    #       - acc/don gain >0.1, find paired acc/don: if gain>0.1, elif gain>0 (weak), else original score > 0.2
    #   exon retention/truncation
    #       - acc/don gain > 0.1, 
    #   exon skipping
    #       - canonical acc/don loss, paired loss of the canonical exon splice site

    def format_chrom_pos(chrom,absposA,absposB):
        if absposA > absposB:
            absposA, absposB = absposB, absposA
        formA = f"{absposA[0]}_ins{absposA[1]}" if absposA[1]>0 else f"{absposA[0]}"
        formB = f"{absposB[0]}_ins{absposB[1]}" if absposB[1]>0 else f"{absposB[0]}"
        return f"{chrom}:{formA}-{formB}"
    
    pos=int(pos) #casting

    resL =[]

    strand = sD['annotation']['chr_strand']
    tD = h[sD['annotation']['Ensembl_nuc']][0]

    for key in ['acceptor_gain_top5','acceptor_loss_top5','donor_gain_top5','donor_loss_top5']:
        sD[key]= [e for e in sD[key] if e['delta']>=0.01] # filter out having only MES changes
        for e in sD[key]:
            e['abspos'] = relpos2abspos(pos, e['relpos'])

    # pseudoexon inclusion / exitron
    for acc, don in itertools.product(sD['acceptor_gain_top5'],sD['donor_gain_top5']):
        if both_in_an_intron(acc['abspos'], don['abspos'], tD):
            poscons = format_chrom_pos(chrom, acc['abspos'],don['abspos']) + f":AG_{acc['delta']:.3f}_DG_{don['delta']:.3f}" + ':inclusion'
            molcons = predict_splMolCons(chrom, pos, ref, alt, tD, acc['relpos'], don['relpos'], splicing_type='inclusion')
            resL.append(':'.join([poscons,molcons]))
        if isExitron(acc['abspos'], don['abspos'], tD) and len_abspos(acc['abspos'], don['abspos'], pos, ref, alt)>30:
            if strand=='+':
                don_, acc_ = shift_right_abspos(don['abspos'],pos, ref, alt), shift_left_abspos(acc['abspos'],pos, ref, alt)
            else:
                acc_, don_ = shift_right_abspos(acc['abspos'],pos, ref, alt), shift_left_abspos(don['abspos'],pos, ref, alt)
            poscons = format_chrom_pos(chrom, acc_, don_) + f":AG_{acc['delta']:.3f}_DG_{don['delta']:.3f}" + ':exitron'
            molcons = predict_splMolCons(chrom, pos, ref, alt, tD, acc['relpos'], don['relpos'], splicing_type='exitron')
            resL.append(':'.join([poscons,molcons]))

    # exon skipping 
    for acc, don in itertools.product(sD['acceptor_loss_top5'],sD['donor_loss_top5']):
        if both_in_an_exon(acc['abspos'], don['abspos'], tD):
            poscons = format_chrom_pos(chrom, acc['abspos'],don['abspos']) + f":AL_{acc['delta']:.3f}_DL_{don['delta']:.3f}" + ':exclusion'
            molcons = predict_splMolCons(chrom, pos, ref, alt, tD, acc['relpos'], don['relpos'], splicing_type='exclusion')
            resL.append(':'.join([poscons,molcons]))

    # exon retention/truncation
    alD = {e['abspos'][0]:e for e in sD['acceptor_loss_top5']}
    al_posL = alD.keys()
    for acc_new in sD['acceptor_gain_top5']:
        can_acce = current_region(acc_new['abspos'],tD['exnList'])
        if can_acce!=(0,0): #truncation
            can_acce_pos = can_acce[0] if strand=='+' else can_acce[1]
            delta_format = f":AG_{acc_new['delta']:.3f}"
            if can_acce_pos in al_posL:
                delta_format += f"_AL_{alD[can_acce_pos]['delta']:.3f}"
            
            if strand=='+':
                acc_ = shift_left_abspos(acc_new['abspos'], pos, ref, alt)
            else:
                acc_ = shift_right_abspos(acc_new['abspos'], pos, ref, alt)
            
            poscons = format_chrom_pos(chrom, acc_, (can_acce_pos,0)) + delta_format + ':truncation'
            molcons = predict_splMolCons(chrom, pos, ref, alt, tD, abspos2relpos(pos, acc_), abspos2relpos(pos, (can_acce_pos,0)), splicing_type='exclusion')
            resL.append(':'.join([poscons,molcons]))

        can_acci = current_region(acc_new['abspos'],tD['intron'])
        if can_acci!=(0,0): #retention
            can_acci_pos = can_acci[1] if strand=='+' else can_acci[0]
            #if abs(acc_new['abspos'][0]-can_acci_pos)>500: # threasholding
            if len_abspos(acc_new['abspos'], (can_acci_pos,0), pos, ref, alt)>500:
                continue
            delta_format = f":AG_{acc_new['delta']:.3f}"
            if can_acci_pos in al_posL:
                delta_format += f"_AL_{alD[can_acci_pos]['delta']:.3f}"
            poscons = format_chrom_pos(chrom, acc_new['abspos'],(can_acci_pos,0)) + delta_format + ':retention'
            molcons = predict_splMolCons(chrom, pos, ref, alt, tD, acc_new['relpos'], abspos2relpos(pos, (can_acci_pos,0)), splicing_type='inclusion')
            resL.append(':'.join([poscons,molcons]))

    dlD = {e['abspos'][0]:e for e in sD['donor_loss_top5']}
    dl_posL = dlD.keys()
    for don_new in sD['donor_gain_top5']:
        can_done = current_region(don_new['abspos'],tD['exnList'])
        if can_done!=(0,0): #truncation
            can_done_pos = can_done[1] if strand=='+' else can_done[0]
            delta_format = f":DG_{don_new['delta']:.3f}"
            if can_done_pos in dl_posL:
                delta_format += f"_DL_{dlD[can_done_pos]['delta']:.3f}"
            
            if strand=='+':
                don_ = shift_right_abspos(don_new['abspos'], pos, ref, alt)
            else:
                don_ = shift_left_abspos(don_new['abspos'], pos, ref, alt)
            
            poscons = format_chrom_pos(chrom, don_, (can_done_pos,0)) + delta_format + ':truncation'
            molcons = predict_splMolCons(chrom, pos, ref, alt, tD, abspos2relpos(pos, don_), abspos2relpos(pos, (can_done_pos,0)), splicing_type='exclusion')
            resL.append(':'.join([poscons,molcons]))

        can_doni = current_region(don_new['abspos'],tD['intron'])
        if can_doni!=(0,0): #retention
            can_doni_pos = can_doni[0] if strand=='+' else can_doni[1]
            #if abs(don_new['abspos'][0]-can_doni_pos)>500: # threasholding
            if len_abspos(don_new['abspos'], (can_doni_pos,0), pos, ref, alt)>500:
                continue
            delta_format = f":DG_{don_new['delta']:.3f}"
            if can_doni_pos in dl_posL:
                delta_format += f"_DL_{dlD[can_doni_pos]['delta']:.3f}"
            poscons = format_chrom_pos(chrom, don_new['abspos'],(can_doni_pos,0)) + delta_format + ':retention'
            molcons = predict_splMolCons(chrom, pos, ref, alt, tD, don_new['relpos'], abspos2relpos(pos, (can_doni_pos,0)), splicing_type='inclusion')
            resL.append(':'.join([poscons,molcons]))
            #resL.append([poscons, delta_format, 'retention', molcons])

    resL = sorted(resL, key=lambda res: sum(map(float,[res.split(':')[2].split('_')][1::2])))

    return '/'.join(resL)

