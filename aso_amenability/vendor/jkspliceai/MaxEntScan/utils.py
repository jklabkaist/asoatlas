import sys
import os
import re
import traceback

from pkg_resources import resource_filename
from jklib import DATA

MES_path = resource_filename(__name__,'')

def makescorematrix_5(mefile='me2x5'):
    me_path = f'{MES_path}/{mefile}'
    matrix = {}
    with open(me_path) as me_file:
        for index, item in enumerate(me_file.readlines()):
            matrix[index]=item.replace('\n', '').replace(' ', '')
    return matrix

def makesequencematrix_5(splice='splicemodels/splice5sequences'):
    splice_path = f'{MES_path}/{splice}'
    matrix = {}
    with open(splice_path) as splice_file:
        for index, item in enumerate(splice_file.readlines()):
            matrix[item.replace('\n','')] = index
    return matrix

def getrest_5(seq):
    seq = list(seq.upper()[:9])
    seq = ''.join(seq)
    seq = seq[:3] + seq[5:9]
    return seq

def scoreconsensus_5(seq):
    seqa = list(seq.upper())
    bgd = {'A':0.27,'C':0.23, 'G':0.23, 'T': 0.27}
    cons1 = {'A' : 0.004, 'C':0.0032,'G':0.9896,'T':0.0032 }
    cons2 = {'A' : 0.0034,'C':0.0039,'G':0.0042,'T':0.9884 }
    addscore = cons1[seqa[3]]*cons2[seqa[4]]/(bgd[seqa[3]]*bgd[seqa[4]])
    return addscore


def makemaxentscores_3(splice_model='splicemodels/'):
    model_path = f'{MES_path}/{splice_model}'

    me_list = ['me2x3acc1','me2x3acc2','me2x3acc3','me2x3acc4','me2x3acc5','me2x3acc6','me2x3acc7','me2x3acc8','me2x3acc9']
    metables = []
    for me in me_list:
        with open(f'{model_path}/{me}') as me_file:
            new_dic={}
            for index, me_line in enumerate(me_file.readlines()):
                new_dic[index] = float(me_line.replace('\n', ''))
            metables.append(new_dic)
    return metables

def getrest_3(seq):
    seq_noconsensus = seq[:18] + seq[20:23]
    return seq_noconsensus

def scoreconsensus_3(seq):
    seq = list(seq.upper())
    bgd = {'A' : 0.27, 'C' : 0.23, 'G' :0.23, 'T': 0.27 }
    cons1 = {'A':  0.9903, 'C' : 0.0032, 'G' : 0.0034, 'T' :0.0030 }
    cons2 = {'A' : 0.0027,'C' :0.0037, 'G' : 0.9905, 'T' : 0.0030 }
    addscore = (cons1[seq[18]] * cons2[seq[19]]) / (bgd[seq[18]] * bgd[seq[19]])
    return addscore

def hashseq(seq):
    seq = seq.upper()
    gene = {'A' : 0, "C" : 1, "T" : 3, "G" : 2 }
    seqa = [ gene[item] for item in list(seq)]
    four = [1,4,16,64,256,1024,4096,16384]

    line_index = 0
    for index in range(len(seq)):
        line_index += (seqa[index] * four[len(seq) - index  -1 ])
    return line_index

def maxentscore_3(seq, metables):
    sc = []
    sc.append(metables[0][hashseq(seq[0:7])])
    sc.append(metables[1][hashseq(seq[7:14])])
    sc.append(metables[2][hashseq(seq[14:21])])
    sc.append(metables[3][hashseq(seq[4:11])])
    sc.append(metables[4][hashseq(seq[11:18])])
    sc.append(metables[5][hashseq(seq[4:7])])
    sc.append(metables[6][hashseq(seq[7:11])])
    sc.append(metables[7][hashseq(seq[11:14])])
    sc.append(metables[8][hashseq(seq[14:18])])
    finalscore = sc[0] * sc[1] * sc[2] * sc[3] * sc[4] / (sc[5] * sc[6] * sc[7] * sc[8])
    return finalscore

