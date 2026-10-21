import os, re
import numpy as np
from enum import Enum

import math
from jkspliceai.MaxEntScan.utils import makescorematrix_5, makesequencematrix_5, getrest_5, scoreconsensus_5, makemaxentscores_3, getrest_3, scoreconsensus_3, maxentscore_3

from jklib.utils.basic import *
from jklib.genome import locus

def float_format(x):
    return float("{0:0.2f}".format(x))

class MaxEntScan2:
    
    def __init__(self):
        self.me2x5 = makescorematrix_5()
        self.seq_matrix = makesequencematrix_5()
        self.metables = makemaxentscores_3()

    def mes(self, seq, score):
        """
        make the score for give maxentscan
        
        Parameters
        ----------
        seq: str or list
            DNA sequence
        score: int or str
            mes score selector

        Returns
        -------
        float or list
            maxentscan score
        """
        if '3' in str(score):
            if len(seq) != 23:
                print (seq)
                raise Exception
            
            try:
                result_dic = {}
                #metables = makemaxentscores_3()
                if type(seq) == list:
                    for s in seq:
                        addscore = scoreconsensus_3(s) * maxentscore_3(getrest_3(s), self.metables)
                        result = math.log(addscore, 2)
                        result = float(f'{result:.2f}')
                        result_dic[seq] = result
                    return result_dic
                elif type(seq) == str:
                    addscore = scoreconsensus_3(seq) * maxentscore_3(getrest_3(seq), self.metables)
                    result = math.log(addscore, 2)
                    return float(f'{result:.2f}')
                else:
                    raise TypeError

            except Exception as e:
                return e


        elif '5' in str(score):
            if len(seq) != 9:
                print (seq)
                raise Exception
                
            try:
                #me2x5 = makescorematrix_5()
                #seq_matrix = makesequencematrix_5()
                if type(seq) == list:
                    result_dic = {}
                    for s in seq:
                        result = math.log(scoreconsensus_5(s) * float(self.me2x5[self.seq_matrix[getrest_5(s)]]),2)
                        result_dic[s] = float(f'{result:.2f}')
                    return result_dic
                elif type(seq) == str:
                    result = math.log(scoreconsensus_5(seq) * float(self.me2x5[self.seq_matrix[getrest_5(seq)]]),2)
                    return float(f'{result:.2f}')
                else:
                    raise TypeError
            except Exception as e:
                return e
