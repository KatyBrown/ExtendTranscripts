#!/usr/bin/env python3
from skbio.alignment import StripedSmithWaterman
import re
import numpy as np
import math


def SWalign(seq1, seq2, pD):
    '''
    Aligns two sequences using the StripedSmithWaterman algorithm
    as implements in skbio.
    All default parameters for this function are used except for
    gap opening penalty, gap extension penalty, match score and
    mismatch score, which are user defined.
    
    Parameters
    ----------
    seq1: str
        string of sequence 1
    seq2: str
        string of sequence 2
    min_score: int
        minimum alignment score to consider
    min_length: int
        minimum alignment length to consider
    gap_open: int
        gap opening penalty
    gap_extend: int
        gap extension penalty
    match_score: int
        match score
    mismatch_score: int
        mismatch score
    
    Returns
    -------
    dict:
        Dictionary with the following alignment outputs from
        StripedSmithWaterman
        optimal_alignment_score': optimal alignment score under SW
        query_start: start position of the optimal alignment in seq1 
        query_end': end position of the optimal alignment in seq1 
        target_start': start position of the optimal alignment in seq2
        target_end': end position of the optimal alignment in seq2
        cigar: cigar string describing the alignment as at
        http://genome.sph.umich.edu/wiki/SAM
    '''   
    ali = StripedSmithWaterman(seq1, gap_open_penalty=pD['gap_open'],
                               gap_extend_penalty=pD['gap_extend'],
                               match_score=pD['match_score'],
                               mismatch_score=pD['mismatch_score'])(seq2)

    aliD = {'optimal_alignment_score': ali.optimal_alignment_score,
            'query_start': ali.query_begin,
            'query_end': ali.query_end + 1,
            'target_start': ali.target_begin,
            'target_end': ali.target_end_optimal + 1,
            'cigar': ali.cigar}
    return (aliD)


def readCIGAR(cigar):
    return(re.findall(r'\d+\D*', cigar))

def lengthFromCIGAR(cigar, excludeI=False, excludeD=False):
    if excludeI:
        return (sum([int(x[:-1]) for x in re.findall(r'\d+[M|D]', cigar)]))
    if excludeD:
        return (sum([int(x[:-1]) for x in re.findall(r'\d+[M|I]', cigar)]))
    
    return(sum([int(x) for x in re.findall(r'\d+', cigar)]))


def alignFromCIGAR(seq1, seq2, cigar):
    '''
    I = insertion in query
    D = insertion in target
    '''
    if cigar == "%sM" % lengthFromCIGAR(cigar):
        return (seq1, seq2)
    cigar_L = readCIGAR(cigar)
    s1 = []
    s2 = []
    n1 = 0
    n2 = 0
    for cig in cigar_L:
        typ = cig[-1]
        count = int(cig[:-1])
        for i in range(count):
            if typ == "M":
                s1.append(seq1[n1])
                s2.append(seq2[n2])
                n1 += 1
                n2 += 1
            elif typ == "I":
                s1.append(seq1[n1])
                s2.append("-")
                n1 += 1
            elif typ == "D":
                s1.append("-")
                s2.append(seq2[n2])
                n2 += 1
            else:
                print (typ)
    s1 = "".join(s1)
    s2 = "".join(s2)
    if (len(s1) != len(s2)):
        raise RuntimeError ("Aligned sequences are not the same length.")
    return (s1, s2)

def getClips(result, query_seq, target_seq, pD):
    Q_length = len(query_seq)
    T_length = len(target_seq)

    Q_start_nclipped = result['query_start']
    T_start_nclipped = result['target_start']

    Q_end_nclipped = Q_length - result['query_end']
    T_end_nclipped = T_length - result['target_end']
    
    if pD['clip_perc']:
        max_clip_Q = math.ceil(pD['clip_perc'] * Q_length)
        max_clip_T = math.ceil(pD['clip_perc'] * T_length)
    else:
        max_clip_Q, max_clip_T = pD['max_clip'], pD['max_clip']
    return (Q_start_nclipped, T_start_nclipped,
            Q_end_nclipped, T_end_nclipped,
            max_clip_Q, max_clip_T)
    

def findOverlapType(result, query_seq, target_seq, pD, keep_failed=False):
    '''
    Add various additional annotations to the SW alignment output
    to help determine if the sequences overlap and how much.
    '''
    (Q_start_nclipped, T_start_nclipped, Q_end_nclipped, T_end_nclipped,
     max_clip_Q, max_clip_T) = getClips(result, query_seq, target_seq, pD)

    # If the starts of both sequences are clipped there is a mismatched region
    if Q_start_nclipped > max_clip_Q and T_start_nclipped > max_clip_T:
        if not keep_failed:
            return None
        else:
            start_type = 4
    elif T_start_nclipped <= max_clip_T and Q_start_nclipped <= max_clip_Q:
        # Start type 1 - target start and query start not clipped
        start_type = 1    
    
    elif T_start_nclipped <= max_clip_T:
        # Start type 2 - query start is clipped but target start is not >
        # query overhangs at 5' end
        start_type = 2

    elif Q_start_nclipped <= max_clip_Q:
        # Start type 3 - target start is clipped but query start is not >
        # target overhangs at 5' end
        start_type = 3
    else:
        raise RuntimeError ("Unexpected overlap at sequence start.")
        
        
    # If the ends of both sequences are clipped there is a mismatched region
    if Q_end_nclipped > max_clip_Q and T_end_nclipped > max_clip_T:
        if not keep_failed:
            return None
        else:
            end_type = 4
    elif T_end_nclipped <= max_clip_T and Q_end_nclipped <= max_clip_Q:
        # End type 1 - target end and query end not clipped
        end_type = 1    
    
    elif T_end_nclipped <= max_clip_T:
        # End type 2 - query end is clipped but target end is not >
        # query overhangs at 3' end
        end_type = 2

    elif Q_end_nclipped <= max_clip_Q:
        # End type 3 - target end is clipped but query end is not >
        # target overhangs at 3' end
        end_type = 3
    else:
        raise RuntimeError ("Unexpected overlap at sequence end.")
    return (start_type, end_type)


def getAlignmentLocal(result, query_seq, target_seq, pD, keep_failed=False):
    
    Q_subseq = query_seq[result['query_start']:result['query_end']]
    T_subseq = target_seq[result['target_start']:result['target_end']]
    
    Q_subseq, T_subseq = alignFromCIGAR(Q_subseq, T_subseq, result['cigar'])
    
    alignment_length = lengthFromCIGAR(result['cigar'])
    Q_subseq = np.array(list(Q_subseq))
    T_subseq = np.array(list(T_subseq))
    
    ident = sum(Q_subseq == T_subseq) / alignment_length
    if ident >= pD['min_perc_ident'] or keep_failed:
        return(Q_subseq, T_subseq, ident)
    else:
        return (None, None, None)

 
def getAlignmentFull(result, query_seq, target_seq, pD):
    subseqQ = "".join(list(result['alignment'][0]))
    subseqT = "".join(list(result['alignment'][1]))
    
    (Q_start_nclipped, T_start_nclipped, Q_end_nclipped, T_end_nclipped,
     max_clip_Q, max_clip_T) = getClips(result, query_seq, target_seq, pD)
    
    if Q_start_nclipped <= max_clip_Q:
        subseqQ = "".join("-" * Q_start_nclipped) + subseqQ[Q_start_nclipped:]
        Q_start_nclipped = 0
    if T_start_nclipped <= max_clip_T:
        subseqT = "".join("-" * T_start_nclipped) + subseqT[T_start_nclipped:]
        T_start_nclipped = 0

    if Q_end_nclipped <= max_clip_Q:
        subseqQ = "".join("-"*Q_end_nclipped) + subseqQ[Q_end_nclipped:]
        Q_end_nclipped = 0
    if T_end_nclipped <= max_clip_T:
        subseqT = "".join("-" * T_end_nclipped) + subseqT[T_end_nclipped:]
        T_end_nclipped = 0
    
    Qpad_start = "".join("-" * T_start_nclipped)
    Tpad_start = "".join("-" * Q_start_nclipped)
    
    Qpad_end = "".join("-" * T_end_nclipped)
    Tpad_end = "".join("-" * Q_end_nclipped)
    
    fullseqQ = Qpad_start + query_seq[:result['query_start']] + subseqQ + query_seq[result['query_end']:] + Qpad_end
    fullseqT = Tpad_start + target_seq[:result['target_start']] + subseqT + target_seq[result['target_end']:] + Tpad_end
    print (fullseqQ)
    print (fullseqT)