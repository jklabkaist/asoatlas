import os, sys, re

import pickle


def margin(innerRange,outerRange): # 0-base
    """
    Parameters
    ----------
    innerRange: tuple
    outerRange: tuple

    Results
    -------
    
    """
    i_s,i_e = innerRange
    o_s,o_e = outerRange

    return i_s-o_s, o_e-i_e


def overlap(x, y): # based on UCSC coordinate
    """
    default overlap processing

    Parameters
    ----------
    x: locus
    y: locus
    
    Returns
    -------
    e-s or 0: int
    """
    c1,s1,e1 = x.chrom, x.chrSta, x.chrEnd
    c2,s2,e2 = y.chrom, y.chrSta, y.chrEnd

    if c1 != c2 or e2<s1 or e1<s2:
        return 0

    s = max(s1,s2)
    e = min(e1,e2)

    if s <= e:
        return e-s+1
    else:
        return 0


class Lookup:

    def __init__(self, pileupDirL, sId, chrom, loc, ref, alt, flag='T'):

        self.pileupDirL = pileupDirL
        self.sId = sId
        self.chrom = chrom
        self.loc = loc
        self.ref = ref
        self.alt = alt
        self.flag = flag


    def lookupPileup(self):
        """
        Parameters
        ----------
        None

        Results
        -------
        (altCount, refCount): tuple

        """
        inputFileNL = []

        if flag == 'T':
            for pileupDir in self.pileupDirL:
                inputFileNL += os.popen('find %s -name %s_T_*%s.pileup_proc' % (pileupDir, self.sId, self.chrom)).readlines()
        else:
            for pileupDir in pileupDirL:
                inputFileNL += os.popen('find %s -name %s_N_*%s.pileup_proc' % (pileupDir, self.sId, self.chrom)).readlines()
                inputFileNL += os.popen('find %s -name %s_B_*%s.pileup_proc' % (pileupDir, self.sId, self.chrom)).readlines()

        if len(inputFileNL) > 1:
            inputFileNL = [x for x in inputFileNL if not re.match('.*KN.*', x)]

        if len(inputFileNL) == 0:
            return None

        resultL = os.popen('grep -m 1 "^%s:%s," %s' % (self.chrom, self.loc,inputFileNL[0].rstrip()), 'r').readlines()

        if len(resultL)==0:
            return None
        else:
            tL = resultL[0].rstrip().split(',')
            if self.ref != tL[2]:
                sys.exit(1)
            refCount = int(tL[3])
            altCount = tL[4].count(self.alt)
            return (altCount, refCount)


    ## batch version of lookupPileup() without sample id
    ## output dictionary of 'sample':'refCount|altCount'
    def lookupPileup_batch(self, useFlag=True):
        """
        batch version of lookupPileup() without sample id
        output dictionary of 'sample':'refCount|altCount'
        
        Parameters
        ----------
        useFlag: boolean
            default value is True
        
        Returns
        -------
        resultH: dictionary
            key is 'sampN'
            value is ref|alt count or 'NA'
        """
        ## same critera as mutScan
        minCover = 3
        minMutReads = 2
        minFreq = 0.01

        inputFileNL = []

        if useFlag:
            if self.flag == 'T':
                for pileupDir in self.pileupDirL:
                    inputFileNL += os.popen('find %s -name *_T_*%s.pileup_proc' % (pileupDir, self.chrom)).readlines()
            else:
                for pileupDir in self.pileupDirL:
                    inputFileNL += os.popen('find %s -name *_N_*%s.pileup_proc' % (pileupDir, self.chrom)).readlines()
                    inputFileNL += os.popen('find %s -name *_B_*%s.pileup_proc' % (pileupDir, self.chrom)).readlines()
        else:
            for pileupDir in self.pileupDirL:
                inputFileNL += os.popen('find %s -name *%s.pileup_proc' % (self.pileupDir, self.chrom)).readlines()

        if len(inputFileNL) > 1:
            inputFileNL = [x for x in inputFileNL if not re.match('.*KN.*', x)]

        if len(inputFileNL) == 0:
            return None

        resultH = {}
        for inputFile in inputFileNL:
            sampN = inputFile.rstrip().split('/')[-1].split('_')[0]
            resultL = os.popen('grep -m 1 -P "^%s\\t%s\\t" %s' % (self.chrom, self.loc, inputFile.rstrip()), 'r').readlines()

            if len(resultL) > 0:
                tL = resultL[0].rstrip().split('\t')
                tot = int(tL[3])
                altCount = tL[4].upper().count(self.alt)
                refCount = tot - altCount
                resultH[sampN] = '%s|%s' % (refCount, altCount)
            else:
                resultH[sampN] = 'NA'

        return resultH    