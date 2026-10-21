import collections, gzip, re
import twobitreader

from pkg_resources import resource_filename
from jklib.utils.basic import rc
from jklib import REFERENCE


def loadBlatOutput(blatOutputPath,by='transID',blacklist=['NR_106988']):
    """
    Parameters
    ----------
    blatOutputPath: string
    by: string
        user select the target name
        default value is 'transID'
    blacklist: array
        delete the blacklist
        default value is ['NR_106988']

    Results
    -------
    h: collections.defaultdict
        target are 'txnSta', 'txnEnd'
    """

    h = collections.defaultdict(list)

    if blatOutputPath.endswith('.gz'):
        f = gzip.open(blatOutputPath)
    else:
        f = open(blatOutputPath)

    for line in f:

        if line[0] == '#':
            continue

        r = _processBlatLine(line)

        if r['transID'] in blacklist:
            continue

        h[r[by]].append(r)

    from operator import itemgetter, attrgetter

    for k,vL in list(h.items()):
        h[k] = sorted(vL,key=itemgetter('txnSta','txnEnd'))

    return h


def _processBlatLine(line):
    """
    user not use this method
    coordinates: 0-base [,)
    you need to convert coordinate system to use this in jklib functions

    Parameters
    ----------
    line: string

    Results
    -------
    h : dictionary

    """
    tokL = line.rstrip().split('\t')

    h = {}

    h['transName'] = tokL[0]
    h['transID'] = tokL[1]
    h['chrom'] = tokL[2]
    h['chrNum'] = tokL[2][3:]
    h['strand'] = tokL[3]
    h['txnSta'] = int(tokL[4])
    h['txnEnd'] = int(tokL[5])
    h['txnLen'] = h['txnEnd'] - h['txnSta']
    h['cdsSta'] = int(tokL[6])
    h['cdsEnd'] = int(tokL[7])
    h['exnList'] = list(map(lambda x,y: (int(x),int(y)), tokL[9].split(',')[:-1], tokL[10].split(',')[:-1]))
    h['exnLenList'] = [e-s for (s,e) in h['exnList']]
    h['exnLen'] = sum(h['exnLenList'])

    if len(tokL) > 12:
        h['geneName'] = tokL[12]

    h['cdsList'] = []
    frontL, backL = [],[]

    if h['cdsSta'] != h['cdsEnd']:

        for (s,e) in h['exnList']:

            frontL.append((min(s,h['cdsSta']),min(e,h['cdsSta'])))
            h['cdsList'].append((max(s,h['cdsSta']),min(e,h['cdsEnd'])))
            backL.append((max(s, h['cdsEnd']), max(e, h['cdsEnd'])))

        frontL = [x for x in frontL if x[0] < x[1]]
        h['cdsList'] = [x for x in h['cdsList'] if x[0] < x[1]]
        backL = [x for x in backL if x[0] < x[1]]

    if h['strand'] == '+':
        h['utr5'] = frontL
        h['utr3'] = backL
    elif h['strand'] == '-':
        h['utr5'] = backL
        h['utr3'] = frontL
    else:
        raise Exception

    h['utr5Len'] = sum([e-s for (s,e) in h['utr5']])
    h['utr3Len'] = sum([e-s for (s,e) in h['utr3']])
    h['cdsLen'] = sum([e-s for (s,e) in h['cdsList']])

    h['intron'] = []

    for i in range(len(h['exnList'])-1):
        h['intron'].append((h['exnList'][i][1],h['exnList'][i+1][0]))

    return h



class locus: # UCSC type

    def __init__(self,loc,id='',assembly='hg38'):
        """
        '[chr]:[sta]-[end][+-]' / '[chr]:[sta]-[end]' are valid only
        for the latter format, the strand value is None.
        locus coordinate is same to UCSC: 1-base, both closed
        """

        rm = re.match('([^:]+):([0-9,]+)-([0-9,]+)([+-])',loc)
        if rm:
            self.strand = rm.group(4)

        else:
            rm = re.match('([^:]+):([0-9,]+)-([0-9,]+)',loc)

            if rm:
                self.strand = None
            else:
                raise Exception

        self.chrom = rm.group(1)
        if self.chrom[:3] == 'chr':
            self.chrNum = rm.group(1)[3:]
        else:
            self.chrNum = self.chrom
            self.chrom = 'chr'+self.chrom

        self.chrSta = int(rm.group(2))
        self.chrEnd = int(rm.group(3))

        self.id = id
        self.assembly = assembly

    def toString(self):
        """
        Parameters
        ----------

        Returns
        -------
        (chrom, chrSta, chrEnd, strand): string
            return to format
        """
        if self.strand:
            return '%s:%s-%s%s' % (self.chrom,self.chrSta,self.chrEnd,self.strand)
        else:
            return '%s:%s-%s' % (self.chrom,self.chrSta, self.chrEnd)

    def twoBitFrag(self, assembly=None, buffer5p=0, buffer3p=0):
        """
        Paramters
        ---------
        assembly: string
            default value is 'hg38'
        buffer5p: int
            default value is 0
        beffer3p: int
            default value is 0

        Returns
        -------
        seq: string
            DNA sequence
        jkbio.rc(seq): string
            DNA sequence
        """

        return self.twoBitFrag2(assembly, buffer5p, buffer3p)

    def twoBitFrag2(self, assembly=None, buffer5p=0, buffer3p=0):
        """
        Paramters
        ---------
        assembly: string
            if user select the twoBitFrag method, assembly is `hg38`
        buffer5p: int
            default value is 0
        beffer3p: int
            default value is 0

        Returns
        -------
        seq: string
            DNA sequence
        jkbio.rc(seq): string
            DNA sequence
        """
        if self.strand == '+' or not self.strand:
            staPos = self.chrSta - 1 - buffer5p
            endPos = self.chrEnd + buffer3p
        elif self.strand == '-':
            staPos = self.chrSta - 1 - buffer3p
            endPos = self.chrEnd + buffer5p

        if not assembly:
            assembly = self.assembly

        if assembly=='hg38':
            hg382bit = twobitreader.TwoBitFile(resource_filename(__name__, f'datafiles/{assembly}.2bit'))
            seq = str(hg382bit[self.chrom][staPos:endPos]) 
            if self.strand == '+' or not self.strand:
                return seq
            else:
                return rc(seq)

        # if assembly is 'hg19',
        twoBitFilePath = resource_filename(__name__, f'datafiles/{assembly}.2bit')
        fragFile = twobitreader.TwoBitFile(twoBitFilePath)
        seq = str(fragFile[self.chrom][staPos:endPos]).upper()

        if self.strand == '+' or not self.strand:
            return seq
        else:
            return rc(seq)
