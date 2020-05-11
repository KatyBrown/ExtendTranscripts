#!/usr/bin/env python3
import numpy as np
import operator

def FastaToDict(infile):
    '''
    Converts a fasta file to a dictionary
    
    Keys are the sequence names (without ">") and values are the
    sequences.

    Parameters
    ----------
    infile: str
        path to fasta file
        
    Returns
    -------
    dict
        dictionary where keys are sequence names and values are sequence or
        sequence lengths
    '''
    D = dict()
    seq = []
    nam = ""
    x = 0
    with open(infile) as input:
        for line in input:
            line = line.strip()
            if len(line) != 0:
                if line[0] == ">":
                    if len(seq) != 0:
                        seq = "".join(seq)
                        D[nam] = seq
                        seq = []
                    nam = line.replace(">", "")

                else:
                    seq.append(line)
            x += 1
    seq = "".join(seq)
    if x != 0:
        D[nam] = seq
    return D

def makeConsensus(F):
    seqs = np.array([list(seq.upper()) for seq in F.values]())
    consensus = []
    for i in range(0,len(seqs[0,:])):
        unique, counts = np.unique(seqs[:,i], return_counts=True)
        unique_ng = unique[unique != "-"]
        counts_ng = counts[unique != "-"]
        count_ng = dict(zip(unique_ng, counts_ng))
        maxChar_ng, maxCount_ng = max(count_ng.items(),
                                      key=operator.itemgetter(1))
        consensus.append(maxChar_ng)
    return ("".join(consensus))
