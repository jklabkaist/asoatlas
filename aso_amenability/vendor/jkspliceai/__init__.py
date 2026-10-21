from jkspliceai.utils import *
from jkspliceai.MaxEntScan import *
import numpy as np
from collections import defaultdict
#import multiprocessing as mp
#import tqdm

#sai = SpliceAI()
sai = SpliceAI2()
#sai_tfl = SpliceAI_TFlite()
mes= MaxEntScan2()

rc = lambda seq: seq.translate(seq.maketrans('acgtACGTN','tgcaTGCAN'))[::-1]
outline = lambda xL: '\t'.join(map(str,xL))+'\n'

def trim_ref_alt(ref,alt):
    # trim : "chr11-10001-GA-GAA" > "chr11-10002-A-AA"
    # trim : "chr11-10001-GAAAAT-TAAAAT" > "chr11-10001-G-T"
    # how to use: 
    #       ref_, alt_, offset = trim_ref_alt(ref,alt)
    #       pos = int(pos)+offset
    
    if ref.startswith(alt) or alt.startswith(ref):
        firstN = len(list(zip(ref,alt)))
        ref_, alt_ = ref[firstN-1:], alt[firstN-1:]
        offset = firstN-1
        return ref_, alt_, offset
        
    if len(ref)==len(alt):
        for i in range(len(ref)):
            if ref[i:]==alt[i:]:
                return ref[:i], alt[:i], 0
                
    return ref, alt, 0

def spliceAI(var, max_distance=5000, dev=False, MES=False):

    """
    get spliceAI score for a single variant
    """

    resultD = defaultdict(list)

    chrom, pos, ref, alt = varstr_parse(var)
    tmp = sai.get_sequence(chrom, pos, ref, alt)

    varstrL, annL, x_refL, x_altL, seq_refL, seq_altL = map(list,zip(*tmp))

    strandL = [ann['chr_strand'] for ann in annL]

    y_refL = sai.predict_value_batch(x_refL, strandL)
    y_altL = sai.predict_value_batch(x_altL, strandL)

    recL = [process_by_ref_alt(ref, alt, y_ref, y_alt, top_k=top_k, max_distance=max_distance, verbose=False, dev=dev) for y_ref, y_alt in zip(y_refL, y_altL)]
    
    if MES:
        pass

    for varstr, ann, rec in zip(varstrL, annL, recL):
        recD = {'annotation':ann}
        recD.update(rec)
        resultD[varstr].append(recD)

    return resultD

def spliceAI_onsite(varL, tidL, max_distance=0, batch_size=500, outpath=None):

    batchL = batch((varL, tidL), batch_size=batch_size)
    y_refL = []

    if outpath:
        out = open(outpath, 'w')
        out.write('\t'.join(['#CHROM','POS','TID','Neither','Acceptor','Donor'])+'\n')

    for var_B, tid_B in batchL:

        annB, x_refB = [],[]

        for var, tid in zip(var_B, tid_B):
            chrom, pos, ref, alt = varstr_parse(var)
            tmp = sai.get_sequence(chrom, pos, ref, alt, tid=tid, max_distance=max_distance)

            _, annL_tmp, x_refL_tmp, _, _, _ = map(list,zip(*tmp))

            annB.extend(annL_tmp)
            x_refB.extend(x_refL_tmp)

        strandB = [ann['chr_strand'] for ann in annB]

        if outpath:
            out.write(''.join( [outline( list(varstr_parse(varstr)[:2]) + [tid] + [f"{score:.4f}" for score in list(y_score[0][0])] ) for varstr, tid, y_score in zip(var_B, tid_B, sai.predict_value_batch(x_refB, strandB)) ] ))
            continue

        y_refL.extend(sai.predict_value_batch(x_refB, strandB))

    if outpath:
        out.close()
        return []

    return [(varstr,y_ref) for varstr, y_ref in zip(varL, y_refL)]

def spliceAI_batch(varL, tidL=None, max_distance=5000, top_k=5, batch_size=500, dev=False, verbose=False, MES=False, nproc=32):
    # "chrom-pos-ref-alt"

    # get a list of annotation, ref_seq, alt_ref for each variant

    varstrL, annL, x_refL, x_altL = [],[],[],[]
    seq_refL, seq_altL = [], []
    cntL=[]
    resultD = defaultdict(list)
    if not tidL:
        tidL = [None] * len(varL)

    for var, tid in zip(varL, tidL):
        chrom, pos, ref, alt = varstr_parse(var)
        tmp = sai.get_sequence(chrom, pos, ref, alt, tid=tid, max_distance=max_distance)

        if tmp is None:
            continue

        try:
            varstrL_tmp, annL_tmp, x_refL_tmp, x_altL_tmp, seq_refL_tmp, seq_altL_tmp = map(list,zip(*tmp))
        except Exception as e:
            print("Error!", e)
            print(var, tid, tmp)
            continue

        varstrL.extend(varstrL_tmp), annL.extend(annL_tmp)
        x_refL.extend(x_refL_tmp), x_altL.extend(x_altL_tmp)
        seq_refL.extend(seq_refL_tmp), seq_altL.extend(seq_altL_tmp)
        cntL.append(len(varstrL_tmp))

    #print("normal",list(map(lambda x: x.upper(), seq_refL)), list(map(lambda x: x.upper(), seq_altL)))
    #raise Exception

    strandL = [ann['chr_strand'] for ann in annL]

    y_refL, y_altL = [], []
    batchL = batch((x_refL, x_altL, strandL), batch_size=batch_size)

    for x_refB, x_altB, strandB in batchL:

        y_refL.extend(sai.predict_value_batch(x_refB, strandB))
        y_altL.extend(sai.predict_value_batch(x_altB, strandB))

    recL = [process_by_ref_alt(*varstr_parse(varstr)[-2:], y_ref, y_alt, top_k=top_k, max_distance=max_distance) for varstr, y_ref, y_alt in zip(varstrL, y_refL, y_altL)]    

    if MES: # When true, update MES score for each spliceAI called data
        #with mp.Pool(processes=32) as pool:
        #    recL = pool.starmap(update_MES_score, tqdm.tqdm(zip(recL, seq_refL, seq_altL, strandL), total=len(recL)))
        #multiprocessing.Process(target=lambda args:update_MES_score(*args), args=zip(recL, seq_refL, seq_altL))
        
        #recL = [update_MES_score(rec, seq_ref, seq_alt, strand) for rec, seq_ref, seq_alt, strand in zip(recL, seq_refL, seq_altL, strandL)]

        recL_ = []
        for rec, seq_ref, seq_alt, strand in zip(recL, seq_refL, seq_altL, strandL):
            try:
                recL_.append(update_MES_score(rec, seq_ref, seq_alt, strand))
            except:
                print('rec',rec)
                print('ref',seq_ref)
                print('alt',seq_alt)
                print('strand',strand)
                raise Exception
            
        recL = recL_

    for varstr, ann, rec in zip(varstrL, annL, recL):
        recD = {'annotation':ann}
        recD.update(rec)
        resultD[varstr].append(recD)

    return resultD

def spliceAI_batch_w_mutation(varL, mut, tidL=None, max_distance=5000, top_k=5, batch_size=500, dev=False, verbose=False, MES=False, nproc=32):
    # "chrom-pos-ref-alt"

    # get a list of annotation, ref_seq, alt_ref for each variant

    varstrL, annL, x_refL, x_altL = [],[],[],[]
    seq_refL, seq_altL = [], []
    cntL=[]
    resultD = defaultdict(list)
    if not tidL:
        tidL = [None] * len(varL)

    for var, tid in zip(varL, tidL):
        chrom, pos, ref, alt = varstr_parse(var)
        tmp = sai.get_sequence_w_mutation(chrom, pos, ref, alt, mut, tid=tid, max_distance=max_distance)

        if tmp is None:
            continue

        try:
            varstrL_tmp, annL_tmp, x_refL_tmp, x_altL_tmp, seq_refL_tmp, seq_altL_tmp = map(list,zip(*tmp))
        except Exception as e:
            print("Error!", e)
            print(var, tid, tmp)
            continue

        varstrL.extend(varstrL_tmp), annL.extend(annL_tmp)
        x_refL.extend(x_refL_tmp), x_altL.extend(x_altL_tmp)
        seq_refL.extend(seq_refL_tmp), seq_altL.extend(seq_altL_tmp)
        cntL.append(len(varstrL_tmp))


    #print("w_mutation",list(map(lambda x: x.upper(), seq_refL)), list(map(lambda x: x.upper(), seq_altL)))
    #raise Exception
    strandL = [ann['chr_strand'] for ann in annL]

    y_refL, y_altL = [], []
    batchL = batch((x_refL, x_altL, strandL), batch_size=batch_size)

    for x_refB, x_altB, strandB in batchL:

        y_refL.extend(sai.predict_value_batch(x_refB, strandB))
        y_altL.extend(sai.predict_value_batch(x_altB, strandB))

    recL = [process_by_ref_alt(*varstr_parse(varstr)[-2:], y_ref, y_alt, top_k=top_k, max_distance=max_distance) for varstr, y_ref, y_alt in zip(varstrL, y_refL, y_altL)]    

    if MES: # When true, update MES score for each spliceAI called data
        #with mp.Pool(processes=32) as pool:
        #    recL = pool.starmap(update_MES_score, tqdm.tqdm(zip(recL, seq_refL, seq_altL, strandL), total=len(recL)))
        #multiprocessing.Process(target=lambda args:update_MES_score(*args), args=zip(recL, seq_refL, seq_altL))
        recL = [update_MES_score(rec, seq_ref, seq_alt, strand) for rec, seq_ref, seq_alt, strand in zip(recL, seq_refL, seq_altL, strandL)]

    for varstr, ann, rec in zip(varstrL, annL, recL):
        recD = {'annotation':ann}
        recD.update(rec)
        resultD[varstr].append(recD)

    return resultD

"""
def spliceAI_batch_original_scores(varL, tidL=None, max_distance=5000, top_k=5, batch_size=500, dev=False, verbose=False, MES=False):

    # get a list of annotation, ref_seq, alt_ref for each variant

    varstrL, annL, x_refL, x_altL = [],[],[],[]
    seq_refL, seq_altL = [], []
    cntL=[]
    resultD = defaultdict(list)
    if not tidL:
        tidL = [None] * len(varL)

    for var, tid in zip(varL, tidL):
        chrom, pos, ref, alt = varstr_parse(var)
        tmp = sai.get_sequence(chrom, pos, ref, alt, tid=tid, max_distance=max_distance)

        if tmp is None:
            continue

        try:
            varstrL_tmp, annL_tmp, x_refL_tmp, x_altL_tmp, seq_refL_tmp, seq_altL_tmp = map(list,zip(*tmp))
        except Exception as e:
            print("Error!", e)
            print(var, tid, tmp)
            continue

        varstrL.extend(varstrL_tmp), annL.extend(annL_tmp)
        x_refL.extend(x_refL_tmp), x_altL.extend(x_altL_tmp)
        seq_refL.extend(seq_refL_tmp), seq_altL.extend(seq_altL_tmp)
        cntL.append(len(varstrL_tmp))

    strandL = [ann['chr_strand'] for ann in annL]

    y_refL, y_altL = [], []
    batchL = batch((x_refL, x_altL, strandL), batch_size=batch_size)

    for x_refB, x_altB, strandB in batchL:

        y_refL.extend(sai.predict_value_batch(x_refB, strandB))
        y_altL.extend(sai.predict_value_batch(x_altB, strandB))

    recL = [process_by_ref_alt_original_scores(*varstr_parse(varstr)[-2:], y_ref, y_alt, top_k=top_k, max_distance=max_distance) for varstr, y_ref, y_alt in zip(varstrL, y_refL, y_altL)]    


    if MES: # When true, update MES score for each spliceAI called data
        recL = [update_MES_score(rec, seq_ref, seq_alt, strand) for rec, seq_ref, seq_alt, strand in zip(recL, seq_refL, seq_altL, strandL)]

    for varstr, ann, rec in zip(varstrL, annL, recL):
        recD = {'annotation':ann}
        recD.update(rec)
        resultD[varstr].append(recD)

    return resultD
"""


def spliceAI_batch_saionly(varL, tidL=None, max_distance=5000, top_k=5, batch_size=500):

    # get a list of annotation, ref_seq, alt_ref for each variant

    varstrL, annL, x_refL, x_altL = [],[],[],[]
    seq_refL, seq_altL = [], []
    cntL=[]
    resultD = defaultdict(list)
    if not tidL:
        tidL = [None] * len(varL)

    for var, tid in zip(varL, tidL):
        chrom, pos, ref, alt = varstr_parse(var)
        tmp = sai.get_sequence(chrom, pos, ref, alt, tid=tid, max_distance=max_distance)

        if tmp is None:
            continue

        try:
            varstrL_tmp, annL_tmp, x_refL_tmp, x_altL_tmp, seq_refL_tmp, seq_altL_tmp = map(list,zip(*tmp))
        except Exception as e:
            print("Error!", e)
            print(var, tid, tmp)
            continue

        varstrL.extend(varstrL_tmp), annL.extend(annL_tmp)
        x_refL.extend(x_refL_tmp), x_altL.extend(x_altL_tmp)
        seq_refL.extend(seq_refL_tmp), seq_altL.extend(seq_altL_tmp)
        cntL.append(len(varstrL_tmp))

    strandL = [ann['chr_strand'] for ann in annL]

    y_refL, y_altL = [], []
    batchL = batch((x_refL, x_altL, strandL), batch_size=batch_size)

    for x_refB, x_altB, strandB in batchL:

        y_refL.extend(sai.predict_value_batch(x_refB, strandB))
        y_altL.extend(sai.predict_value_batch(x_altB, strandB))

    recL = [process_by_ref_alt(*varstr_parse(varstr)[-2:], y_ref, y_alt, top_k=top_k, max_distance=max_distance) for varstr, y_ref, y_alt in zip(varstrL, y_refL, y_altL)]    

    for varstr, ann, rec in zip(varstrL, annL, recL):
        recD = {'annotation':ann}
        recD.update(rec)
        resultD[varstr].append(recD)

    return resultD


def spliceAI_batch_saionly_w_mutation(varL, mut, tidL=None, max_distance=5000, top_k=5, batch_size=500):

    # get a list of annotation, ref_seq, alt_ref for each variant

    varstrL, annL, x_refL, x_altL = [],[],[],[]
    seq_refL, seq_altL = [], []
    cntL=[]
    resultD = defaultdict(list)
    if not tidL:
        tidL = [None] * len(varL)

    for var, tid in zip(varL, tidL):
        chrom, pos, ref, alt = varstr_parse(var)
        tmp = sai.get_sequence_w_mutation(chrom, pos, ref, alt, mut, tid=tid)

        if tmp is None:
            continue

        print(tmp[0][4].upper())
        print(tmp[0][5].upper())
        raise Exception

        try:
            varstrL_tmp, annL_tmp, x_refL_tmp, x_altL_tmp, seq_refL_tmp, seq_altL_tmp = map(list,zip(*tmp))
        except Exception as e:
            print("Error!", e)
            print(var, tid, tmp)
            continue

        varstrL.extend(varstrL_tmp), annL.extend(annL_tmp)
        x_refL.extend(x_refL_tmp), x_altL.extend(x_altL_tmp)
        seq_refL.extend(seq_refL_tmp), seq_altL.extend(seq_altL_tmp)
        cntL.append(len(varstrL_tmp))

    strandL = [ann['chr_strand'] for ann in annL]

    y_refL, y_altL = [], []
    batchL = batch((x_refL, x_altL, strandL), batch_size=batch_size)

    for x_refB, x_altB, strandB in batchL:

        y_refL.extend(sai.predict_value_batch(x_refB, strandB))
        y_altL.extend(sai.predict_value_batch(x_altB, strandB))

    recL = [process_by_ref_alt(*varstr_parse(varstr)[-2:], y_ref, y_alt, top_k=top_k, max_distance=max_distance) for varstr, y_ref, y_alt in zip(varstrL, y_refL, y_altL)]    

    for varstr, ann, rec in zip(varstrL, annL, recL):
        recD = {'annotation':ann}
        recD.update(rec)
        resultD[varstr].append(recD)

    return resultD

def spliceAI_batch_saionly_tfl(varL, tidL=None, max_distance=5000, top_k=5, batch_size=500):

    # get a list of annotation, ref_seq, alt_ref for each variant

    varstrL, annL, x_refL, x_altL = [],[],[],[]
    seq_refL, seq_altL = [], []
    cntL=[]
    resultD = defaultdict(list)
    if not tidL:
        tidL = [None] * len(varL)

    for var, tid in zip(varL, tidL):
        chrom, pos, ref, alt = varstr_parse(var)
        tmp = sai_tfl.get_sequence(chrom, pos, ref, alt, tid=tid)

        if tmp is None:
            continue

        try:
            varstrL_tmp, annL_tmp, x_refL_tmp, x_altL_tmp, seq_refL_tmp, seq_altL_tmp = map(list,zip(*tmp))
        except Exception as e:
            print("Error!", e)
            print(var, tid, tmp)
            continue

        varstrL.extend(varstrL_tmp), annL.extend(annL_tmp)
        x_refL.extend(x_refL_tmp), x_altL.extend(x_altL_tmp)
        seq_refL.extend(seq_refL_tmp), seq_altL.extend(seq_altL_tmp)
        cntL.append(len(varstrL_tmp))

    strandL = [ann['chr_strand'] for ann in annL]

    y_refL, y_altL = [], []
    batchL = batch((x_refL, x_altL, strandL), batch_size=batch_size)

    for x_refB, x_altB, strandB in batchL:

        y_refL.extend(sai_tfl.predict_value_batch(x_refB, strandB))
        y_altL.extend(sai_tfl.predict_value_batch(x_altB, strandB))

    recL = [process_by_ref_alt(*varstr_parse(varstr)[-2:], y_ref, y_alt, top_k=top_k, max_distance=max_distance) for varstr, y_ref, y_alt in zip(varstrL, y_refL, y_altL)]    


    for varstr, ann, rec in zip(varstrL, annL, recL):
        recD = {'annotation':ann}
        recD.update(rec)
        resultD[varstr].append(recD)

    return resultD


def update_MES_score(rec, seq_ref, seq_alt, strand):

    def get_messeq(d_from_variant, seq_ref, seq_alt, strand, isAcceptor=True):

        pos = len(seq_ref)//2
        len_delta = len(seq_alt) - len(seq_ref)
        offset = (-20,3) if isAcceptor else (-2, 7)
        if strand=='-':
            offset = (-2, 21) if isAcceptor else (-6,3)

        if isinstance(d_from_variant, str) and 'ins' in d_from_variant:

            d_from_variant = int(d_from_variant.replace('ins_',''))
            idx_alt = pos + d_from_variant 
            mes_seq_alt = seq_alt[idx_alt + offset[0] :idx_alt + offset[1]]

            return None, mes_seq_alt

        idx_ref = pos + d_from_variant
        idx_alt = pos + d_from_variant

        if d_from_variant>0 : # For different length indel, there will be some shift downstream of base position
            idx_alt = idx_alt + len_delta

        mes_seq_ref = seq_ref[idx_ref + offset[0] :idx_ref + offset[1]]
        mes_seq_alt = seq_alt[idx_alt + offset[0] :idx_alt + offset[1]]
        if strand=='-':
            mes_seq_ref, mes_seq_alt = rc(mes_seq_ref), rc(mes_seq_alt)

        return mes_seq_ref, mes_seq_alt
    
    def get_mes(seq):
        try:
            if seq is None:
                return None
            if len(seq)==23:
                return mes.mes(seq, 3)
            elif len(seq)==9:
                return mes.mes(seq, 5)
            else:
                raise Exception('MaxEntScan: Not valid input')
        except:
            print(seq)
            raise Exception

    def topk_mes(k, dL, seq_ref, seq_alt, strand):
        if 'acceptor' in k:
            isAcc, k_ref, k_alt = True, 'MES_acceptor_ref', 'MES_acceptor_alt'
        elif 'donor' in k:
            isAcc, k_ref, k_alt = False, 'MES_donor_ref', 'MES_donor_alt'
        else:
            raise Exception

        rdL = []
        for d in dL:
            mes_seq_ref, mes_seq_alt = get_messeq(d['relpos'], seq_ref, seq_alt, strand, isAcceptor=isAcc)
            rdL.append(d.update({k_ref: get_mes(mes_seq_ref), k_alt:get_mes(mes_seq_alt)}) or d)

        return rdL

    def stats_mes(statd, seq_ref, seq_alt, strand):
        for k, v in statd.items():
            acc_seq_ref, acc_seq_alt = get_messeq(v['relpos'], seq_ref, seq_alt, strand, isAcceptor=True)
            mes_acc_ref, mes_acc_alt = get_mes(acc_seq_ref), get_mes(acc_seq_alt)
            don_seq_ref, don_seq_alt = get_messeq(v['relpos'], seq_ref, seq_alt, strand, isAcceptor=False)
            mes_don_ref, mes_don_alt = get_mes(don_seq_ref), get_mes(don_seq_alt)

            statd[k].update({'mes_acceptor_ref':mes_acc_ref, 'mes_acceptor_alt':mes_acc_alt,
                             'mes_donor_ref': mes_don_ref,   'mes_donor_alt': mes_don_alt})
            
        #return statd

    #for rec, seq_ref, seq_alt in zip(recL, seq_refL, seq_altL):

    # update topk acc/don MES record
    topkL = [k for k in rec.keys() if 'top' in k]
    rec.update({k:topk_mes(k, rec[k], seq_ref, seq_alt, strand) for k in topkL})

    # update stats acc/don MES record
    if 'stats' in rec.keys():
        stats_mes(rec['stats'], seq_ref, seq_alt, strand)

    return rec

def spliceAI_batch_dev(varL, tidL=None, max_distance=5000, top_k=5, batch_size=500, dev=False, verbose=False, MES=False):

    # get a list of annotation, ref_seq, alt_ref for each variant

    varstrL, annL, x_refL, x_altL = [],[],[],[]
    seq_refL, seq_altL = [], []
    cntL=[]
    resultD = defaultdict(list)
    if not tidL:
        tidL = [None] * len(varL)

    for var, tid in zip(varL, tidL):
        chrom, pos, ref, alt = varstr_parse(var)
        tmp = sai.get_sequence(chrom, pos, ref, alt, tid=tid)

        try:
            varstrL_tmp, annL_tmp, x_refL_tmp, x_altL_tmp, seq_refL_tmp, seq_altL_tmp = map(list,zip(*tmp))
        except Exception as e:
            print("Error!", e)
            print(var, tid, tmp)
            continue

        varstrL.extend(varstrL_tmp), annL.extend(annL_tmp)
        x_refL.extend(x_refL_tmp), x_altL.extend(x_altL_tmp)
        seq_refL.extend(seq_refL_tmp), seq_altL.extend(seq_altL_tmp)
        cntL.append(len(varstrL_tmp))

    strandL = [ann['chr_strand'] for ann in annL]

    y_refL, y_altL = [], []
    batchL = batch((x_refL, x_altL, strandL), batch_size=batch_size)

    for x_refB, x_altB, strandB in batchL:

        y_refL.extend(sai.predict_value_batch(x_refB, strandB))
        y_altL.extend(sai.predict_value_batch(x_altB, strandB))

    recL = [process_by_ref_alt(*varstr_parse(varstr)[-2:], y_ref, y_alt, top_k=top_k, max_distance=max_distance, verbose=verbose, dev=dev) for varstr, y_ref, y_alt in zip(varstrL, y_refL, y_altL)]    

    resL = []
    if MES: # When true, update MES score for each spliceAI called data
        for rec, seq_ref, seq_alt, strand in zip(recL, seq_refL, seq_altL, strandL):
            resL.append(update_MES_score_dev(rec, seq_ref, seq_alt, strand))

    for varstr, ann, rec in zip(varstrL, annL, resL):
        recD = {'annotation':ann}
        recD.update(rec)
        resultD[varstr].append(recD)

    return resultD

def update_MES_score_dev(rec, seq_ref, seq_alt, strand):

    def get_messeq(d_from_variant, seq_ref, seq_alt, strand, isAcceptor=True):

        pos = len(seq_ref)//2
        len_delta = len(seq_alt) - len(seq_ref)
        offset = (-20,3) if isAcceptor else (-2, 7)
        if strand=='-':
            offset = (-2, 21) if isAcceptor else (-6,3)

        if isinstance(d_from_variant, str) and 'ins' in d_from_variant:

            d_from_variant = int(d_from_variant.replace('ins_',''))
            idx_alt = pos + d_from_variant 
            mes_seq_alt = seq_alt[idx_alt + offset[0] :idx_alt + offset[1]]

            return None, mes_seq_alt

        idx_ref = pos + d_from_variant
        idx_alt = pos + d_from_variant

        if d_from_variant>0 : # For different length indel, there will be some shift downstream of base position
            idx_alt = idx_alt + len_delta

        mes_seq_ref = seq_ref[idx_ref + offset[0] :idx_ref + offset[1]]
        mes_seq_alt = seq_alt[idx_alt + offset[0] :idx_alt + offset[1]]
        if strand=='-':
            mes_seq_ref, mes_seq_alt = rc(mes_seq_ref), rc(mes_seq_alt)

        print(d_from_variant, strand, mes_seq_ref, mes_seq_alt)
        return mes_seq_ref, mes_seq_alt
    
    def get_mes(seq):
        try: #
            if seq is None:
                return None
            if len(seq)==23:
                return mes.mes(seq, 3)
            elif len(seq)==9:
                print(seq)
                return mes.mes(seq, 5)
            else:
                raise Exception('MaxEntScan: Not valid input')
        except:#
            print (seq)#
            raise Exception#

    def topk_mes(k, dL, seq_ref, seq_alt, strand):
        if 'acceptor' in k:
            isAcc, k_ref, k_alt = True, 'MES_acceptor_ref', 'MES_acceptor_alt'
        elif 'donor' in k:
            isAcc, k_ref, k_alt = False, 'MES_donor_ref', 'MES_donor_alt'
        else:
            raise Exception

        rdL = []
        for d in dL:
            mes_seq_ref, mes_seq_alt = get_messeq(d['relpos'], seq_ref, seq_alt, strand, isAcceptor=isAcc)
            rdL.append(d.update({k_ref: get_mes(mes_seq_ref), k_alt:get_mes(mes_seq_alt)}) or d)

        return rdL

    def stats_mes(statd, seq_ref, seq_alt, strand):
        for k, v in statd.items():
            acc_seq_ref, acc_seq_alt = get_messeq(v['relpos'], seq_ref, seq_alt, strand, isAcceptor=True)
            mes_acc_ref, mes_acc_alt = get_mes(acc_seq_ref), get_mes(acc_seq_alt)
            don_seq_ref, don_seq_alt = get_messeq(v['relpos'], seq_ref, seq_alt, strand, isAcceptor=False)
            mes_don_ref, mes_don_alt = get_mes(don_seq_ref), get_mes(don_seq_alt)

            statd[k].update({'mes_acceptor_ref':mes_acc_ref, 'mes_acceptor_alt':mes_acc_alt,
                             'mes_donor_ref': mes_don_ref,   'mes_donor_alt': mes_don_alt})
            
        #return statd

    #for rec, seq_ref, seq_alt in zip(recL, seq_refL, seq_altL):

    # update topk acc/don MES record
    topkL = [k for k in rec.keys() if 'top' in k]
    rec.update({k:topk_mes(k, rec[k], seq_ref, seq_alt, strand) for k in topkL})

    # update stats acc/don MES record
    if 'stats' in rec.keys():
        stats_mes(rec['stats'], seq_ref, seq_alt, strand)

    return rec
            
            
