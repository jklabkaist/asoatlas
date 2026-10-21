from pkg_resources import resource_filename
from collections import defaultdict
from jklib.genome import locus, loadBlatOutput, REFERENCE
from spliceai.utils import one_hot_encode, load_model, normalise_chrom
import logging
import os, re
import numpy as np
import pandas as pd
import tensorflow as tf
import copy
# codes are adapted from spliceai.utils

def abspos_by_relpos(pos, relpos):
    if isinstance(relpos, str) and 'ins' in relpos:
        return None
    return int(pos) + int(relpos)

def tx_mask(pos, h, s):

    tx = [t for t in h[s['annotation']['Ensembl_nuc']] if normalise_chrom(t['chrom'],s['annotation']['GRCh38_chr'])==s['annotation']['GRCh38_chr'] and t['txnSta']+1==int(s['annotation']['chr_start']) and t['txnEnd']==int(s['annotation']['chr_end']) and t['strand']==s['annotation']['chr_strand']][0]
    acc_boundaryL = [i_s+1 for i_s, i_e in tx['exnList']] #0-based -> 1-based
    don_boundaryL = [i_e for i_s, i_e in tx['exnList']]
    if tx['strand'] == '-':
        acc_boundaryL, don_boundaryL = don_boundaryL, acc_boundaryL

    # skip gain in canonical exon boundaries
    s['acceptor_gain_top5'] = [e for e in s['acceptor_gain_top5'] if abspos_by_relpos(pos, e['relpos']) not in acc_boundaryL]
    s['donor_gain_top5'] = [e for e in s['donor_gain_top5'] if abspos_by_relpos(pos, e['relpos']) not in don_boundaryL]

    # skip loss in non-canonical exon boundaries
    s['acceptor_loss_top5'] = [e for e in s['acceptor_loss_top5'] if abspos_by_relpos(pos, e['relpos']) in acc_boundaryL]
    s['donor_loss_top5'] = [e for e in s['donor_loss_top5'] if abspos_by_relpos(pos, e['relpos']) in don_boundaryL]

    return s

def check_molecular_consequence(s, tD):
    return 

def check_molecular_consequence2(s, tD):
    return

def mask_res_dict(varL, resD, h, tidL=[]):

    tid_match = len(varL) == len(tidL)

    resD_=copy.deepcopy(resD)

    for i, var in enumerate(varL):

        if tid_match:
            resD_[var] = [d for d in resD_[var] if d['annotation']['Ensembl_nuc']==tidL[i]]

        sL = resD_[var]
        pos = int(var.split('-')[1])
        for s in sL:
            s = tx_mask(pos, h, s)

    return resD_

def filter_res_dict(varL, resD, tidL=[], delta_th=0.01, h=None, mask=False):

    resD_ = copy.deepcopy(resD)
    tid_match = len(varL) == len(tidL)

    for i, var in enumerate(varL):
        """
        if tid_match:
            resD[var] = [d for d in resD[var] if d['annotation']['Ensembl_nuc']==tidL[i]]
        """
        sL = resD_[var]
        pos = int(var.split('-')[1])
        for s in sL:

            if h and mask:
                s = tx_mask(pos, h, s)

            s['acceptor_gain_top5'] = [e for e in s['acceptor_gain_top5'] if e['delta'] >= delta_th or e['MES_acceptor_ref']!=e['MES_acceptor_alt']]
            s['acceptor_loss_top5'] = [e for e in s['acceptor_loss_top5'] if e['delta'] >= delta_th or e['MES_acceptor_ref']!=e['MES_acceptor_alt']]
            s['donor_gain_top5'] = [e for e in s['donor_gain_top5'] if e['delta'] >= delta_th or e['MES_donor_ref']!=e['MES_donor_alt']]
            s['donor_loss_top5'] = [e for e in s['donor_loss_top5'] if e['delta'] >= delta_th or e['MES_donor_ref']!=e['MES_donor_alt']]

    return resD_

def filter_res_dict_saionly(varL, resD, tidL=[], delta_th=0.01, h=None, mask=False):

    tid_match = len(varL) == len(tidL)

    for i, var in enumerate(varL):
        """
        if tid_match:
            resD[var] = [d for d in resD[var] if d['annotation']['Ensembl_nuc']==tidL[i]]
        """
        sL = resD[var]
        pos = int(var.split('-')[1])
        for s in sL:

            if h and mask:
                s = tx_mask(pos, h, s)

            s['acceptor_gain_top5'] = [e for e in s['acceptor_gain_top5'] if e['delta'] >= delta_th]
            s['acceptor_loss_top5'] = [e for e in s['acceptor_loss_top5'] if e['delta'] >= delta_th]
            s['donor_gain_top5'] = [e for e in s['donor_gain_top5'] if e['delta'] >= delta_th]
            s['donor_loss_top5'] = [e for e in s['donor_loss_top5'] if e['delta'] >= delta_th]

    return resD

sai_keyL = ['acceptor_gain_top5', 'acceptor_loss_top5', 'donor_gain_top5', 'donor_loss_top5']
def resD_to_list(varL, resD, top_k=5):

    def MES_ref_alt(l):
        refk = [k for k in l.keys() if k in ['MES_acceptor_ref', 'MES_donor_ref']][0]
        altk = [k for k in l.keys() if k in ['MES_acceptor_alt', 'MES_donor_alt']][0]
        mes_ref = l[refk] if l[refk] else '.'
        mes_alt = l[altk] if l[altk] else '.'
        return f"{refk[:-4]}({mes_ref}>{mes_alt})"

    def format_list(lL, top_k):
        if not len(lL):
            return ('.', '0','.')
        return ( '/'.join([str(l['relpos']) for l in lL[:top_k]]), '/'.join(["{:.2f}".format(l['delta']) for l in lL[:top_k]]), '/'.join([MES_ref_alt(l) for l in lL[:top_k]]) )

    def flaten_res(xL):
        dim = len(xL[0])
        ret = []
        for i in range(dim):
            ret += [x[i] for x in xL]
        return ret

    return {var+';'+s['annotation']['Ensembl_nuc']: [s['annotation']['Ensembl_nuc']] + flaten_res([format_list(s[k], top_k=top_k) for k in sai_keyL]) for var in varL for s in resD[var]}

def resD_to_list_refalt(varL, resD, top_k=5):

    def MES_ref_alt(l):
        refk = [k for k in l.keys() if k in ['MES_acceptor_ref', 'MES_donor_ref']][0]
        altk = [k for k in l.keys() if k in ['MES_acceptor_alt', 'MES_donor_alt']][0]
        mes_ref = l[refk] if l[refk] else '.'
        mes_alt = l[altk] if l[altk] else '.'
        return f"{refk[:-4]}({mes_ref}>{mes_alt})"

    def format_list(lL, top_k):
        if not len(lL):
            return ('.', '0','.')
        return ( '/'.join([str(l['relpos']) for l in lL[:top_k]]), '/'.join(["{:.3f}>{:.3f}".format(l['ref'],l['alt']) for l in lL[:top_k]]), '/'.join([MES_ref_alt(l) for l in lL[:top_k]]) )

    def flaten_res(xL):
        dim = len(xL[0])
        ret = []
        for i in range(dim):
            ret += [x[i] for x in xL]
        return ret

    return {var+';'+s['annotation']['Ensembl_nuc']: [s['annotation']['Ensembl_nuc']] + flaten_res([format_list(s[k], top_k=top_k) for k in sai_keyL]) for var in varL for s in resD[var]}

def resD_to_list_saionly(varL, resD, top_k=5):

    def format_list(lL, top_k):
        if not len(lL):
            return ('.', '0')
        return ( '/'.join([str(l['relpos']) for l in lL[:top_k]]), '/'.join(["{:.2f}".format(l['delta']) for l in lL[:top_k]]) )

    def flaten_res(xL):
        dim = len(xL[0])
        ret = []
        for i in range(dim):
            ret += [x[i] for x in xL]
        return ret

    return {var+';'+s['annotation']['Ensembl_nuc']: [s['annotation']['Ensembl_nuc']] + flaten_res([format_list(s[k], top_k=top_k) for k in sai_keyL]) for var in varL for s in resD[var]}


def parse_MES_text(txt): # MES_{}({}>{})
    pattern = r"\((.*?)>(.*?)\)"
    matches = re.findall(pattern, txt)
    return matches[0] if len(matches)==1 else ('.','.')

sai_headerL = ['spliceAI_5k_transcriptID']
for prefix in ['spliceAI_5k_raw_', 'spliceAI_5k_masked_']:
    for suf1 in ['score_', 'relpos_', 'MES_']:
        for suf2 in ['AG','AL','DG','DL']:
            sai_headerL.append(prefix + suf1 + suf2)



h=loadBlatOutput(resource_filename(__name__,'annotations/gencode.v41.annotation.gtf.refFlat.MANE.v1.arranged.txt'), by='transID')
def get_resD_masked(varL, resD, tidL=[], delta_th=0.1, top_k=1, h=h, sep='\t'):
    #resD = filter_res_dict(varL, resD, tidL =tidL, delta_th=delta_th)
    maskD = mask_res_dict(varL, resD, h)
    return {var+';'+s['annotation']['Ensembl_nuc']:s for var in varL for s in maskD[var]}

def outform_vcf_info(varL, resD, tidL=[], delta_th=0.01, top_k=1, h=h, sep='\t'):
    # raw
    resD = filter_res_dict(varL, resD, tidL =tidL, delta_th=delta_th)
    rawD = resD_to_list(varL, resD, top_k=top_k)
    # masked
    resD = mask_res_dict(varL, resD, h)
    maskedD = resD_to_list(varL, resD, top_k=top_k)

    return {vartid: sep.join(rawD[vartid] + maskedD[vartid][1:]) for vartid in rawD.keys()}

def outform_vcf_info_refalt(varL, resD, tidL=[], delta_th=0.1, top_k=1, h=h, sep='\t'):
    # raw
    #resD = filter_res_dict(varL, resD, tidL =tidL, delta_th=delta_th)
    rawD = resD_to_list_refalt(varL, resD, top_k=top_k)
    # masked
    resD = mask_res_dict(varL, resD, h)
    maskedD = resD_to_list_refalt(varL, resD, top_k=top_k)

    return {vartid: sep.join(rawD[vartid] + maskedD[vartid][1:]) for vartid in rawD.keys()}

def outform_vcf_info_saionly(varL, resD, tidL=[], delta_th=0.01, top_k=1, h=h, sep='\t'):
    # raw
    resD = filter_res_dict_saionly(varL, resD, tidL =tidL, delta_th=delta_th)
    rawD = resD_to_list_saionly(varL, resD, top_k=top_k)
    # masked
    resD = mask_res_dict(varL, resD, h)
    maskedD = resD_to_list_saionly(varL, resD, top_k=top_k)

    return {vartid: sep.join(rawD[vartid] + maskedD[vartid][1:]) for vartid in rawD.keys()}

def batch(inL, batch_size=500):
    if isinstance(inL, tuple):
        tmp = [len(e) for e in inL]
        if tmp.count(tmp[0])!=len(tmp):
            raise Exception
        return [tuple([inL[x][i:i+batch_size] for x in range(len(inL))]) for i in range(0, len(inL[0]), batch_size)]
    return [inL[i:i+batch_size] for i in range(0, len(inL), batch_size)]

def divide_list(inL, cntL):
    resL = []
    idx= 0
    for cnt in cntL:
        resL.append(inL[idx:idx+cnt])
        idx += cnt
    return resL

def varstr_parse(var):
    chrom, pos, ref, alt = var.split('-')
    return normalise_chrom(chrom, 'chrN'), int(pos), ref.upper(), alt.upper()

def MANE_parser():
    
    annotations = resource_filename(__name__, 'annotations/MANE.GRCh38.v1.0.summary.txt') #1-based
    # mane grch38 v1.0 https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.0/MANE.GRCh38.v1.0.summary.txt.gz
    # #NCBI_GeneID, Ensembl_Gene, HGNC_ID, symbol, name, RefSeq_nuc, RefSeq_prot, Ensembl_nuc, Ensembl_prot, MANE_status, GRCh38_chr, chr_start, chr_end, chr_strand
    resD = {}
    for line in open(annotations):
        if line.startswith('#NCBI'):
            headerL = line.strip(' \r\n').split('\t')
            continue
        
        tokL = line.strip(' \r\n').split('\t')
        resD[tokL[7]] = {k: v for k, v in zip(headerL[1:], tokL[1:])}

    return resD

def one_hot_encode_N2avg(seq):

    map = np.asarray([[0.25, 0.25, 0.25, 0.25],
                      [1, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])

    seq = seq.upper().replace('A', '\x01').replace('C', '\x02')
    seq = seq.replace('G', '\x03').replace('T', '\x04').replace('N', '\x00')

    return map[np.fromstring(seq, np.int8) % 5]

def one_hot_encode(seq,dtype=int):

    map = np.asarray([[0, 0, 0, 0],
                      [1, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]], dtype=dtype)

    seq = seq.upper().replace('A', '\x01').replace('C', '\x02')
    seq = seq.replace('G', '\x03').replace('T', '\x04').replace('N', '\x00')

    return map[np.fromstring(seq, np.int8) % 5]

def replace_ref_to_alt(dna_sequence,pos,ref,alt,max_distance):
    original_dna_sequence=dna_sequence
    altered_dna_sequence=dna_sequence[:5000+max_distance] + alt + dna_sequence[5000+max_distance+len(ref):]
    return original_dna_sequence,altered_dna_sequence

def padding_nonfixed_arrL(arrL, pad=[0,0,0,0]):
    maxLeninBatch = max([a[0].shape[0] for a in arrL])
    resL = []
    for arr in arrL:
        if arr[0].shape[0]==maxLeninBatch:
            resL.append(arr)
            continue
        
        resL.append(np.concatenate((arr, np.asarray([[pad]*(maxLeninBatch-arr[0].shape[0])])), axis=1))
    
    return resL

class SpliceAI:

    def __init__(self):

        try:
            self.MANE = MANE_parser() # 1-based
            self.chrom_pat = list(self.MANE.values())[0]['GRCh38_chr']

        except IOError as e:
            logging.error('{}'.format(e))
            exit()

        paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
        self.models = [load_model(resource_filename(__name__, x)) for x in paths]

    def get_MANE_annotations(self, chrom, pos, tid=None):

        chrom = normalise_chrom(chrom, self.chrom_pat)
        pos = int(pos)

        if tid:
            #1-based
            return [d for d in self.MANE.values() if d['GRCh38_chr']==chrom and d['Ensembl_nuc']==tid and (int(d['chr_start']) <= pos) and (pos <= int(d['chr_end']))]

        return [d for d in self.MANE.values() if d['GRCh38_chr']==chrom and (int(d['chr_start']) <= pos) and (pos <= int(d['chr_end']))]


    def get_sequence(self, chrom, pos, ref, alt, tid=None, max_distance=5000):

        """
        input: VCF form attributes, i.e., position, ref/alt allele based on postive strand, 1-based, GRCh38 genome
        output: a list of (MANE annotation dictionary, x_ref, x_alt). x_ref/alt are np arrays with one-hot encoded 
        x_ref_, x_alt_ are input DNA (ATCGN) sequence in + strand direction (not RCed)
        """

        pos = int(pos)

        varStr = f"{chrom}-{pos}-{ref}-{alt}"

        anndL = self.get_MANE_annotations(chrom, pos, tid=tid)

        cov = 2*max_distance+1
        wid = 10000 + cov # 10k model

        start_index=pos-5000-max_distance
        end_index = pos+5000+max_distance
        
        seq = locus(f"{chrom}:{str(start_index)}-{str(end_index)}+").twoBitFrag()

        resL = []
        for annD in anndL:

            gene, tid_ens = annD['symbol'], annD['Ensembl_nuc']
            tstart, tend, tstrand = int(annD['chr_start']), int(annD['chr_end']), annD['chr_strand']
            dist_tx_start, dist_tx_end = (tstart - pos), (tend - pos)

            pad_size = [max(wid//2 + dist_tx_start, 0), max(wid//2-dist_tx_end, 0)]
            ref_len = len(ref)
            alt_len = len(alt)

            x_ref_ = 'N'*pad_size[0] + seq[pad_size[0]:wid-pad_size[1]] + 'N'*pad_size[1]
            x_alt_ = x_ref_[:wid//2] + str(alt) + x_ref_[wid//2+ref_len:]

            x_ref = one_hot_encode(x_ref_)[None, :]
            x_alt = one_hot_encode(x_alt_)[None, :]
        
            if tstrand == '-':
                x_ref = x_ref[:, ::-1, ::-1]
                x_alt = x_alt[:, ::-1, ::-1]

            resL.append((varStr, annD, x_ref, x_alt, x_ref_, x_alt_))

        return resL

    def get_sequence_multivar(self,):
        pass

    def predict_value(self, x, strand):
        """
        x: an one_hot encoded sequence array
        strand: '+' or '-'
        """

        y = np.mean([self.models[m].predict(x) for m in range(5)], axis=0)
        if strand=='-':
            y = y[:, ::-1]
        
        return y

    def predict_value_batch(self, xL, strandL, batch_size=32):

        """
        xL: a list of one_hot encoded sequence array
        strandL: a list of strand (+/-)
        """
        resL = []

        lenL = [x[0].shape[0]-10000 for x in xL]

        x_v = np.concatenate(padding_nonfixed_arrL(xL), axis=0)
        y_v = np.mean([self.models[m].predict_on_batch(x_v) for m in range(5)], axis=0)

        for y, l, strand in zip(y_v, lenL, strandL):

            y = np.asarray([y[:l]]) # trimming padded regions
            if strand=='-':
                y = y[:,::-1]

            resL.append(y)
        
        return resL

    def predict_value_batch2(self, x, strandL, batch_size=32):

        """
        xL: a list of one_hot encoded sequence array
        strandL: a list of strand (+/-)
        """
        
        x_vL = []
        for i in range(x.shape[1]-10000):
            x_vL.extend(x[:,i:i+10001,:])


        x_v = np.concatenate(padding_nonfixed_arrL(xL), axis=0)
        y_v = np.mean([self.models[m].predict_on_batch(x_v) for m in range(5)], axis=0)

        for y, l, strand in zip(y_v, lenL, strandL):

            y = np.asarray([y[:l]]) # trimming padded regions
            if strand=='-':
                y = y[:,::-1]

            resL.append(y)
        
        return resL

class SpliceAI2:

    def __init__(self,model_idx=['1','2','3','4','5']):

        try:
            self.MANE = MANE_parser()
            self.chrom_pat = list(self.MANE.values())[0]['GRCh38_chr']

        except IOError as e:
            logging.error('{}'.format(e))
            exit()

        paths = ('models/spliceai{}.h5'.format(x) for x in model_idx)
        self.models = [load_model(resource_filename(__name__, x)) for x in paths]
        self.MANEseqD = {tid:{'tid':tid, 'chrom':tD['GRCh38_chr'], 'start': int(tD['chr_start']), 'end': int(tD['chr_end']), 'strand': tD['chr_strand'], 'dnaseq':''} for tid, tD in self.MANE.items()}
        #MANE.summary: 1-based start,end
		
    def get_MANE_annotations(self, chrom, pos, tid=None):

        chrom = normalise_chrom(chrom, self.chrom_pat)
        pos = int(pos)

        if tid:
            d=self.MANE[tid]
            if d['GRCh38_chr']==chrom and d['Ensembl_nuc']==tid and (int(d['chr_start']) <= pos) and (pos <= int(d['chr_end'])):
                return [d]
            else:
                []
            #return [d for d in self.MANE.values() if d['GRCh38_chr']==chrom and d['Ensembl_nuc']==tid and (int(d['chr_start']) <= pos) and (pos <= int(d['chr_end']))]

        return [d for d in self.MANE.values() if d['GRCh38_chr']==chrom and (int(d['chr_start']) <= pos) and (pos <= int(d['chr_end']))]

    def get_sequence(self,chrom,pos,ref,alt,tid,max_distance=5000):

        pos = int(pos)
        ref_len, alt_len = len(ref), len(alt)

        varStr = f"{chrom}-{pos}-{ref}-{alt}"

        mane_chrom, mane_sta, mane_end, mane_strand = self.MANEseqD[tid]['chrom'], self.MANEseqD[tid]['start'], self.MANEseqD[tid]['end'], self.MANEseqD[tid]['strand']
        
        if self.MANEseqD[tid]['dnaseq'] == '':
            self.MANEseqD[tid]['dnaseq'] = locus(f"{mane_chrom}:{str(mane_sta)}-{str(mane_end)}+").twoBitFrag()
        mane_dnaseq = self.MANEseqD[tid]['dnaseq']

        """
        try:
            annD = self.get_MANE_annotations(chrom, pos, tid=tid)[0]
        except:
            print(f"Exception - out-of-range: {varStr}")
            return None
        """
        try:
            annD = self.MANE[tid]
        except:
            print(f"Exception - Not MANE transcript: {tid}, {varStr}")
            return None

        cov = 2*max_distance+1
        wid = 10000 + cov # 10k model

        start_index=pos-5000-max_distance
        end_index = pos+5000+max_distance

        start_offset = max(start_index-mane_sta, 0)
        up_padding_n = max(mane_sta-start_index,0)
        end_offset = min(end_index-mane_sta+1,mane_end-mane_sta+1)
        down_padding_n = max(end_index-mane_end,0)

        x_ref_ = 'N'*up_padding_n + mane_dnaseq[start_offset:end_offset] + 'N'*down_padding_n
        x_alt_ = x_ref_[:wid//2] + str(alt) + x_ref_[wid//2+ref_len:]

        x_ref = one_hot_encode(x_ref_)[None, :]
        x_alt = one_hot_encode(x_alt_)[None, :]
    
        if mane_strand == '-':
            x_ref = x_ref[:, ::-1, ::-1]
            x_alt = x_alt[:, ::-1, ::-1]
            
        return[(varStr, annD, x_ref, x_alt, x_ref_, x_alt_)]

    def get_sequence_w_mutation(self,chrom,pos,ref,alt,mut,tid,max_distance=5000):

        def apply_mutation(ref_dna, mut): #single snv onnly
            chrom_, pos_, ref_, alt_ = mut.split('-')
            pos_=int(pos_)

            mut_offset = pos_ - self.MANEseqD[tid]['start']
            #print(mut_offset, len(ref_dna), pos_, self.MANEseqD[tid]['start'])

            if not(self.MANEseqD[tid]['start'] <= pos_ <= self.MANEseqD[tid]['end']):
                #print("not in the gene")
                return ref_dna

            if ref_dna[mut_offset].upper() != ref_.upper():
                print (ref_dna[mut_offset-10:mut_offset], ref_dna[mut_offset], ref_dna[mut_offset+1:mut_offset+10])
                raise Exception

            return ref_dna[:mut_offset] + alt_ + ref_dna[mut_offset+len(ref_):]
            

        pos = int(pos)
        ref_len, alt_len = len(ref), len(alt)

        varStr = f"{chrom}-{pos}-{ref}-{alt}"

        mane_chrom, mane_sta, mane_end, mane_strand = self.MANEseqD[tid]['chrom'], self.MANEseqD[tid]['start'], self.MANEseqD[tid]['end'], self.MANEseqD[tid]['strand']
        
        if self.MANEseqD[tid]['dnaseq'] == '':
            self.MANEseqD[tid]['dnaseq'] = locus(f"{mane_chrom}:{str(mane_sta)}-{str(mane_end)}+").twoBitFrag()
        mane_dnaseq = apply_mutation(self.MANEseqD[tid]['dnaseq'], mut).upper()

        """
        try:
            annD = self.get_MANE_annotations(chrom, pos, tid=tid)[0]
        except:
            print(f"Exception - out-of-range: {varStr}")
            return None
        """
        try:
            annD = self.MANE[tid]
        except:
            print(f"Exception - Not MANE transcript: {tid}, {varStr}")
            return None

        cov = 2*max_distance+1
        wid = 10000 + cov # 10k model

        start_index=pos-5000-max_distance
        end_index = pos+5000+max_distance

        start_offset = max(start_index-mane_sta, 0)
        up_padding_n = max(mane_sta-start_index,0)
        end_offset = min(end_index-mane_sta+1,mane_end-mane_sta+1)
        down_padding_n = max(end_index-mane_end,0)

        x_ref_ = 'N'*up_padding_n + mane_dnaseq[start_offset:end_offset] + 'N'*down_padding_n
        x_alt_ = x_ref_[:wid//2] + str(alt) + x_ref_[wid//2+ref_len:]

        x_ref = one_hot_encode(x_ref_)[None, :]
        x_alt = one_hot_encode(x_alt_)[None, :]
    
        if mane_strand == '-':
            x_ref = x_ref[:, ::-1, ::-1]
            x_alt = x_alt[:, ::-1, ::-1]

        #print(x_ref_, x_alt_)
        #raise Exception

        return[(varStr, annD, x_ref, x_alt, x_ref_, x_alt_)]

    def get_sequence_multivar(self,):
        pass

    def predict_value(self, x, strand):
        """
        x: an one_hot encoded sequence array
        strand: '+' or '-'
        """

        y = np.mean([self.models[m].predict_on_batch(x) for m in range(5)], axis=0)
        if strand=='-':
            y = y[:, ::-1]
        
        return y

    def predict_value_batch(self, xL, strandL):
        
        """
        xL: a list of one_hot encoded sequence array
        strandL: a list of strand (+/-)
        """

        resL = []

        lenL = [x[0].shape[0]-10000 for x in xL]

        x_v = np.concatenate(padding_nonfixed_arrL(xL), axis=0)
        y_v = np.mean([self.models[m].predict_on_batch(x_v) for m in range(5)], axis=0)

        for y, l, strand in zip(y_v, lenL, strandL):

            y = np.asarray([y[:l]]) # trimming padded regions
            if strand=='-':
                y = y[:,::-1]

            resL.append(y)
        
        return resL

    def predict_value_batch2(self, x, strandL):

        """
        xL: a list of one_hot encoded sequence array
        strandL: a list of strand (+/-)
        """
        
        x_vL = []
        for i in range(x.shape[1]-10000):
            x_vL.extend(x[:,i:i+10001,:])


        x_v = np.concatenate(padding_nonfixed_arrL(xL), axis=0)
        y_v = np.mean([self.models[m].predict(x_v) for m in range(5)], axis=0)

        for y, l, strand in zip(y_v, lenL, strandL):

            y = np.asarray([y[:l]]) # trimming padded regions
            if strand=='-':
                y = y[:,::-1]

            resL.append(y)
        
        return resL

class TFLite_model:
    def __init__(self,model_path="",gpu=None):
        self.model_path=model_path
        
        self.interpreter = tf.lite.Interpreter(model_path=self.model_path)
        if gpu:
            gpu_delegate = tf.lite.experimental.Delegate('libtensorflowlite_gpu_delegate.so', gpu_delegate_options)
        #self.interpreter.allocate_tensors

        self.input_details = self.get_input_details()
        self.output_details= self.get_output_details()

    def get_input_details(self,):
        return self.interpreter.get_input_details()

    def get_output_details(self,):
        return self.interpreter.get_output_details()

    def reshape_input_layer(self,input_shape):
        # e.g. SNV, cov 10001 (5k model, 5k window) > input_shape=(1,20002,4)
        self.interpreter.resize_tensor_input(self.input_details[0]['index'], input_shape)
        self.interpreter.allocate_tensors()

        self.input_details = self.get_input_details()
        self.output_details= self.get_output_details()

    def predict(self,input_data):
        self.reshape_input_layer(input_data.shape)
        self.interpreter.set_tensor(self.input_details[0]['index'], input_data)
        self.interpreter.invoke()
        return self.interpreter.get_tensor(self.output_details[0]['index'])

class SpliceAI_TFlite:

    def __init__(self,model_idx=['1','2','3','4','5']):

        try:
            self.MANE = MANE_parser()
            self.chrom_pat = list(self.MANE.values())[0]['GRCh38_chr']

        except IOError as e:
            logging.error('{}'.format(e))
            exit()

        paths = ('models/spliceai{}.tflite'.format(x) for x in model_idx)
        self.models = [TFLite_model(resource_filename(__name__, x)) for x in paths]
        self.MANEseqD = {tid:{'tid':tid, 'chrom':tD['GRCh38_chr'], 'start': int(tD['chr_start']), 'end': int(tD['chr_end']), 'strand': tD['chr_strand'], 'dnaseq':''} for tid, tD in self.MANE.items()}
        #MANE.summary: 1-based start,end
		
    def get_MANE_annotations(self, chrom, pos, tid=None):

        chrom = normalise_chrom(chrom, self.chrom_pat)
        pos = int(pos)

        if tid:
            return [d for d in self.MANE.values() if d['GRCh38_chr']==chrom and d['Ensembl_nuc']==tid and (int(d['chr_start']) <= pos) and (pos <= int(d['chr_end']))]

        return [d for d in self.MANE.values() if d['GRCh38_chr']==chrom and (int(d['chr_start']) <= pos) and (pos <= int(d['chr_end']))]

    def get_sequence(self,chrom,pos,ref,alt,tid,max_distance=5000):

        pos = int(pos)
        ref_len, alt_len = len(ref), len(alt)

        varStr = f"{chrom}-{pos}-{ref}-{alt}"

        mane_chrom, mane_sta, mane_end, mane_strand = self.MANEseqD[tid]['chrom'], self.MANEseqD[tid]['start'], self.MANEseqD[tid]['end'], self.MANEseqD[tid]['strand']
        
        if self.MANEseqD[tid]['dnaseq'] == '':
            self.MANEseqD[tid]['dnaseq'] = locus(f"{mane_chrom}:{str(mane_sta)}-{str(mane_end)}+").twoBitFrag()
        mane_dnaseq = self.MANEseqD[tid]['dnaseq']

        annD = self.get_MANE_annotations(chrom, pos, tid=tid)[0]

        cov = 2*max_distance+1
        wid = 10000 + cov # 10k model

        start_index=pos-5000-max_distance
        end_index = pos+5000+max_distance

        start_offset = max(start_index-mane_sta, 0)
        up_padding_n = max(mane_sta-start_index,0)
        end_offset = min(end_index-mane_sta+1,mane_end-mane_sta+1)
        down_padding_n = max(end_index-mane_end,0)

        x_ref_ = 'N'*up_padding_n + mane_dnaseq[start_offset:end_offset] + 'N'*down_padding_n
        x_alt_ = x_ref_[:wid//2] + str(alt) + x_ref_[wid//2+ref_len:]

        x_ref = one_hot_encode(x_ref_, np.float32)[None, :]
        x_alt = one_hot_encode(x_alt_, np.float32)[None, :]
    
        if mane_strand == '-':
            x_ref = x_ref[:, ::-1, ::-1]
            x_alt = x_alt[:, ::-1, ::-1]

        return[(varStr, annD, x_ref, x_alt, x_ref_, x_alt_)]

    def get_sequence_multivar(self,):
        pass

    def predict_value(self, x, strand):
        """
        x: an one_hot encoded sequence array
        strand: '+' or '-'
        """

        y = np.mean([m.predict(x) for m in self.models], axis=0)
        if strand=='-':
            y = y[:, ::-1]
        
        return y

    def predict_value_batch(self, xL, strandL):

        """
        xL: a list of one_hot encoded sequence array
        strandL: a list of strand (+/-)
        """
        resL = []

        lenL = [x[0].shape[0]-10000 for x in xL]

        x_v = np.concatenate(padding_nonfixed_arrL(xL), axis=0)
        y_v = np.mean([m.predict(x_v) for m in self.models], axis=0)

        for y, l, strand in zip(y_v, lenL, strandL):

            y = np.asarray([y[:l]]) # trimming padded regions
            if strand=='-':
                y = y[:,::-1]

            resL.append(y)
        
        return resL

    def predict_value_batch2(self, x, strandL):

        """
        xL: a list of one_hot encoded sequence array
        strandL: a list of strand (+/-)
        """
        
        x_vL = []
        for i in range(x.shape[1]-10000):
            x_vL.extend(x[:,i:i+10001,:])


        x_v = np.concatenate(padding_nonfixed_arrL(xL), axis=0)
        y_v = np.mean([m.predict(x_v) for m in self.models], axis=0)

        for y, l, strand in zip(y_v, lenL, strandL):

            y = np.asarray([y[:l]]) # trimming padded regions
            if strand=='-':
                y = y[:,::-1]

            resL.append(y)
        
        return resL

def process_by_ref_alt(ref, alt, y_ref, y_alt, top_k=5, max_distance=5000, verbose=False, dev=False):

    
    cov = 2*max_distance + 1
    wid = 10000 + cov

    ref_len, alt_len = len(ref), len(alt)

    del_len = max(len(ref) - len(alt), 0)
    ins_len = max(len(alt) - len(ref), 0)

    y_ins = None

    dist_from_variant = np.arange(-1 * max_distance, 0).tolist()
    dist_from_variant.extend(np.arange(0, max_distance + 1).tolist())

    if ref_len == alt_len:

        y_con = np.concatenate([y_ref, y_alt])

    if ref_len > alt_len:
        y_alt = np.concatenate([
            y_alt[:, :cov//2 +alt_len],
            np.zeros((1, del_len, 3)),
            y_alt[:, cov//2 +alt_len:]],
            axis=1)

        y_con = np.concatenate([y_ref, y_alt])

    if ref_len < alt_len:
        #print(y_alt[:,:cov//2 + ref_len].shape, y_alt[:, cov//2 + alt_len:].shape)
        
        y_ins_ref = np.zeros((1, ins_len, 3))
        y_ins_alt = y_alt[:, cov//2 + ref_len : cov//2 + alt_len]
        
        y_alt = np.concatenate([
            y_alt[:,:cov//2 + ref_len],
            y_alt[:, cov//2 + alt_len:]],
            axis=1)

        if dev: print(y_ins_ref, y_ins_alt)
        
        y_con = np.concatenate([y_ref, y_alt])
        y_ins = np.concatenate([y_ins_ref, y_ins_alt])

        y_con = np.concatenate([y_con, y_ins], axis=1)
        #y_ins_ag = y_ins[1, :, 1] - y_ins[0, :, 1]
        #y_ins_dg = y_ins[1, :, 2] - y_ins[0, :, 2]

        dist_from_variant.extend([f'ins_{i}' for i in np.arange(ref_len,alt_len).tolist()])


    acc = y_con[1, :, 1] - y_con[0, :, 1]
    donor = y_con[1, :, 2] - y_con[0, :, 2]

    acc_ref = y_con[0, :, 1].tolist()
    acc_alt = y_con[1, :, 1].tolist()

    don_ref = y_con[0, :, 2].tolist()
    don_alt = y_con[1, :, 2].tolist()

    record = {}

    if top_k and isinstance(top_k, int) and top_k>0:

        idx_ag, idx_al = np_topk_both_indice(acc, top_k=top_k)
        idx_dg, idx_dl = np_topk_both_indice(donor, top_k=top_k)

        ag = [{'relpos':dist_from_variant[idx], 'ref':acc_ref[idx], 'alt':acc_alt[idx], 'delta':acc[idx]} for idx in idx_ag] # [(relpos, ref_score, alt_score, delta_score)]
        al = [{'relpos':dist_from_variant[idx], 'ref':acc_ref[idx], 'alt':acc_alt[idx], 'delta':-1*acc[idx]} for idx in idx_al]
        dg = [{'relpos':dist_from_variant[idx], 'ref':don_ref[idx], 'alt':don_alt[idx], 'delta':donor[idx]} for idx in idx_dg]
        dl = [{'relpos':dist_from_variant[idx], 'ref':don_ref[idx], 'alt':don_alt[idx], 'delta':-1*donor[idx]} for idx in idx_dl]

        record = {f'acceptor_gain_top{top_k}': ag,
        f'acceptor_loss_top{top_k}': al,
        f'donor_gain_top{top_k}': dg,
        f'donor_loss_top{top_k}': dl}

    if verbose:
        stats={}
        for idx in np.arange(0, len(dist_from_variant)):
            stats[dist_from_variant[idx]] = {
                'relpos': dist_from_variant[idx],
                'acceptor_ref': acc_ref[idx],
                'acceptor_alt': acc_alt[idx],
                'acceptor_delta': acc[idx],
                'donor_ref': don_ref[idx],
                'donor_alt': don_alt[idx],
                'donor_delta': donor[idx]
            }
        
        record.update({'stats': stats})

    return record
"""
def process_by_ref_alt_original_scores(ref, alt, y_ref, y_alt, top_k=5, max_distance=5000, original_scores=True, original_score_th=0.1, dev=False):
    
    cov = 2*max_distance + 1
    wid = 10000 + cov

    ref_len, alt_len = len(ref), len(alt)

    del_len = max(len(ref) - len(alt), 0)
    ins_len = max(len(alt) - len(ref), 0)

    y_ins = None

    dist_from_variant = np.arange(-1 * max_distance, 0).tolist()
    dist_from_variant.extend(np.arange(0, max_distance + 1).tolist())

    if ref_len == alt_len:

        y_con = np.concatenate([y_ref, y_alt])

    if ref_len > alt_len:
        y_alt = np.concatenate([
            y_alt[:, :cov//2 +alt_len],
            np.zeros((1, del_len, 3)),
            y_alt[:, cov//2 +alt_len:]],
            axis=1)

        y_con = np.concatenate([y_ref, y_alt])

    if ref_len < alt_len:
        
        y_ins_ref = np.zeros((1, ins_len, 3))
        y_ins_alt = y_alt[:, cov//2 + ref_len : cov//2 + alt_len]
        
        y_alt = np.concatenate([
            y_alt[:,:cov//2 + ref_len],
            y_alt[:, cov//2 + alt_len:]],
            axis=1)

        if dev: print(y_ins_ref, y_ins_alt)
        
        y_con = np.concatenate([y_ref, y_alt])
        y_ins = np.concatenate([y_ins_ref, y_ins_alt])

        y_con = np.concatenate([y_con, y_ins], axis=1)

        dist_from_variant.extend([f'ins_{i}' for i in np.arange(ref_len,alt_len).tolist()])


    acc = y_con[1, :, 1] - y_con[0, :, 1]
    donor = y_con[1, :, 2] - y_con[0, :, 2]

    acc_ref = y_con[0, :, 1].tolist()
    acc_alt = y_con[1, :, 1].tolist()

    don_ref = y_con[0, :, 2].tolist()
    don_alt = y_con[1, :, 2].tolist()

    record = {}

    if top_k and isinstance(top_k, int) and top_k>0:

        idx_ag, idx_al = np_topk_both_indice(acc, top_k=top_k)
        idx_dg, idx_dl = np_topk_both_indice(donor, top_k=top_k)

        ag = [{'relpos':dist_from_variant[idx], 'ref':acc_ref[idx], 'alt':acc_alt[idx], 'delta':acc[idx]} for idx in idx_ag] # [(relpos, ref_score, alt_score, delta_score)]
        al = [{'relpos':dist_from_variant[idx], 'ref':acc_ref[idx], 'alt':acc_alt[idx], 'delta':-1*acc[idx]} for idx in idx_al]
        dg = [{'relpos':dist_from_variant[idx], 'ref':don_ref[idx], 'alt':don_alt[idx], 'delta':donor[idx]} for idx in idx_dg]
        dl = [{'relpos':dist_from_variant[idx], 'ref':don_ref[idx], 'alt':don_alt[idx], 'delta':-1*donor[idx]} for idx in idx_dl]

        record = {f'acceptor_gain_top{top_k}': ag,
        f'acceptor_loss_top{top_k}': al,
        f'donor_gain_top{top_k}': dg,
        f'donor_loss_top{top_k}': dl}

    if original_scores:
        ori_scores=[]
        for idx in np.arange(0, len(dist_from_variant)):
            if any([ e>=original_score_th for e in [acc_ref[idx], acc_alt[idx], don_ref[idx], don_alt[idx]] ]):
                ori_scores.append({
                    'relpos': dist_from_variant[idx],
                    'acc_ref': acc_ref[idx],
                    'acc_alt': acc_alt[idx],
                    'acc_delta': acc[idx],
                    'don_ref': don_ref[idx],
                    'don_alt': don_alt[idx],
                    'don_delta': donor[idx]
                })
        record.update({'original_scores': ori_scores})

    return record
"""
def get_delta_by_pos(target, pos, ref, alt, y_ref, y_alt, max_distance=5000):
    pass

def np_topk_indice(arr, top_k=5):
    return arr.argsort()[-1*top_k:][::-1]

def np_topk_both_indice(arr, top_k=5):
    idx = arr.argsort()
    return idx[-1*top_k:][::-1], idx[:top_k] # top5, btm5

def process_by_ref_alt_multivar(ratupL, y_ref, y_alt, max_distance=5000):
    pass