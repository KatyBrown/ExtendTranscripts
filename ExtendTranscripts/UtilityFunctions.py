#!/usr/bin/env python3
import numpy as np
import re
import itertools
import sys
import configparser


def makepDTesting(conf_file):
    parser = configparser.ConfigParser()
    parser.read(conf_file)
    confdict = {section: dict(
        parser.items(section)) for section in parser.sections()}
    C = dict()
    for c in confdict['alignment']:
        val = confdict['alignment'][c]
        if val == "None":
            C['alignment_%s' % c] = None
        elif val == "stepwise":
            C['alignment_%s' % c] = val
        else:
            C['alignment_%s' % c] = float(val)
    C['print'] = True
    C['stdout_verbosity'] = 3
    return (C)


def FastaToDict(infile, rna=False):
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


def AlignmentArray(seqs):
    seqs = [list(x.upper()) for x in seqs]
    arr = np.array(seqs)
    return (arr)


def IUPAC():
    '''
    Make two dictionaries representing the IUPAC characters for ambiguous
    bases.
    source: https://droog.gs.washington.edu/parc/images/iupac.html

    Returns
    -------
    D1: dict
        Every IUPAC ambiguity character as keys and lists of all the
        nucleotides they represent as values.
    D2: dict
        Every combination of nucleotides (as strings) as keys and their
        ambiguity characters as values.
    '''
    D1 = {'A': ['A', 'M', 'R', 'W', 'V', 'H', 'D', 'N'],
          'G': ['G', 'R', 'S', 'K', 'V', 'D', 'B', 'N'],
          'T': ['T', 'W', 'Y', 'K', 'H', 'D', 'B', 'N'],
          'C': ['C', 'M', 'S', 'Y', 'V', 'H', 'B', 'N'],
          'M': ['A', 'C'],
          'R': ['A', 'G'],
          'W': ['A', 'T'],
          'S': ['C', 'G'],
          'Y': ['C', 'T'],
          'K': ['G', 'T'],
          'V': ['A', 'C', 'G'],
          'H': ['A', 'C', 'T'],
          'D': ['A', 'G', 'T'],
          'B': ['C', 'G', 'T'],
          'N': ['A', 'C', 'T', 'G']}
    D2 = dict()
    for key, val in D1.items():
        if key not in ['A', 'G', 'T', 'C']:
            for comb in itertools.permutations(val, len(val)):
                D2["".join(comb)] = key
    return (D1, D2)


def getRCDict():
    rcdict = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N", "Y": "R",
              "K": "M", "R": "Y", "M": "K", "B": "V", "V": "B", "D": "H",
              "H": "D", "W": "W", "S": "S", "-": "-"}
    return (rcdict)


def reverseComplement(seq):
    '''
    Reverse complements a sequence

    Parameters
    ----------
    seq: str
        single string containing the sequence

    Returns
    -------
    str
        string containing the reverse complement of the sequence
    '''
    rcdict = getRCDict()
    seq = list(seq)[::-1]
    seq = [rcdict[s] for s in seq]
    seq = "".join(seq)
    return (seq)


def reverseComplementAlignment(ali):
    '''
    Reverse complement a numpy array containing a 2D multiple sequence
    alignment.  Dimension 0 is rows in the alignment, dimension 1 is
    columns.
    '''
    rcdict = getRCDict()
    new = []
    for row in ali:
        rev = row[::-1]
        revC = [rcdict[X] for X in rev]
        new.append(revC)
    new = np.array(new)
    return (new)


def readCIGAR(cigar):
    return(re.findall(r'\d+\D*', cigar))


def lengthFromCIGAR(cigar, excludeI=False, excludeD=False, mOnly=False):
    if excludeI:
        return (sum([int(x[:-1]) for x in re.findall(r'\d+[M|D]', cigar)]))
    if excludeD:
        return (sum([int(x[:-1]) for x in re.findall(r'\d+[M|I]', cigar)]))
    if mOnly:
        return (sum([int(x[:-1]) for x in re.findall(r'\d+[M]', cigar)]))
    return(sum([int(x) for x in re.findall(r'\d+', cigar)]))


def logPrint(string, level, pD):
    if 'log' in pD:
        log = pD['log']
        if level <= pD['logging_verbosity']:
            log.info(string)
    if 'print' in pD:
        if level <= pD['stdout_verbosity']:
            sys.stdout.write("%s\n" % (string))
