import types
import socket
import sys


def mean(l):
    """
    return to mean value of list

    Parameters
    ----------
    l: list

    Returns
    -------
    mean: float
        list mean value

    """
    return float(sum(l)/len(l))

def incHash(h,key,val=1):
    """
    Parameters
    ----------
    h: dictionary
    key: str
    val: int
        default value is 1

    Returns
    -------
    None
        only calculate in dictionary
    """
    try:
        h[key] += val
    except:
        h[key] = val


def addHash(h,key,val):
    """
    append to dictionary

    Parameters
    ----------
    h: dictionary
    key:

    Returns
    -------

    """

    try:
            h[key].append(val)
    except:
            h[key] = [val]


def pushHash(h,key,val):
    """
    add vlaue to dictionary

    Parameters
    ----------
    h: dictionary
    key: str
        dictionary key name
    val: int

    Returns
    -------
    None
        just add value in dictionary
        delete the duplicate value in dictionary
    """
    try:
            h[key].add(val)
    except:
            h[key] = set([val])


def checkType(seq):
    """
    just check type to seq
    bytes or not
    
    Paramters
    ---------
    seq: str
        DNA sequence

    Returns
    -------
    None
        just checking sequence format type
        if sequence is bytes, system close
    """

    if isinstance(seq,bytes):
        print ("error: string type required!")
        print((type(seq)))
        sys.exit(1)


def rev(seq):
    """
    return to change reverse DNA sequence

    Paramters
    ---------
    seq: str
        DNA sequence

    Returns
    -------
    seq_rev: str
        reverse DNA sequence
    
    """
    checkType(seq)

    seq_rev = ""

    for i in range(0,len(seq)):
            seq_rev = seq[i] + seq_rev

    return seq_rev


def compl_v1(seq, seqType):
    """
    translate according to sequence type
    
    Parameters
    ----------
    seq: str
        DNA sequence
    seqType: str
        RNA or not

    Returns
    -------
    seq_ret: str
        translate DNA sequence
    """
    checkType(seq)

    seq_ret = seq.upper()

    if seqType == 'RNA':
            tab = str.maketrans('ATUGCRYSW-','UAACGYRSW-')
    else:
            tab = str.maketrans('ATUGCRYSW-','TAACGYRSW-')
    
    seq_ret = seq_ret.translate(tab)

    return seq_ret


def compl_v2(seq, seqType='DNA'):
    """

    Parameters
    ----------
    seq: str
        DNA or RNA sequence
    seqType: str
        default value is 'DNA'
        selection case is 'DNA' or 'RNA'

    Returns
    -------
    seq_ret: str
        first call array format
        change str format
    """
    wc_DNA = {'N':'N','.':'.','C':'G','G':'C','A':'T','T':'A','*':'*'}
    wc_RNA = {'N':'N','.':'.','C':'G','G':'C','A':'U','U':'A','*':'*'}
    
    checkType(seq)

    seq_ret = []

    for c in seq.upper():
        if seqType=='DNA':
            try:
                tmp = wc_DNA[c]
            except KeyError:
                tmp = c
            seq_ret.append(tmp)
        elif seqType=='RNA':
            try:
                tmp = wc_RNA[c]
            except KeyError:
                tmp = c
            seq_ret.append(tmp)
        else:
            raise Exception

    return ''.join(seq_ret)



def rc(seq,type='DNA'):
    """
    reverse and translate to sequence 

    Parameters
    ----------
    seq: str
        DNA sequence
    type: str
        default value is 'DNA'

    Returns
    -------
    seq: str
        return to using rev and compl sequence value
    """
    checkType(seq)

    seq = rev(seq)
    seq = compl_v1(seq,type)

    return seq

