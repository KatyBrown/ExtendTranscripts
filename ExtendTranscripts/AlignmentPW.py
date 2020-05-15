#!/usr/bin/env python3
from skbio.alignment import StripedSmithWaterman
import numpy as np
import math
import UtilityFunctions
import Consensus


def subMatrixIUPAC(match_score, mismatch_score):
    '''
    Generates a substitution matrix for the SWalign algorithm where
    ambiguous bases are scored equally to non-ambiguous bases.
    In the first step of finding a consensus sequence with just two sequences
    there will be many positions where each base is equally likely - this
    prevents them from being penalised.
    
    Parameters
    ----------
    match_score: int
        Score awarded to a match under the SW model
    mismatch_score: int
        Score awarded to a mismatch under the SW model
    
    Returns
    -------
    scoreD: dict
        2D dictionary of every possible combination of nucleotides, including
        IUPAC characters, and the score to assign to this combination.
    '''
    amb = Consensus.IUPAC()[0]
    scoreD = dict()
    for base1, L in amb.items():
        scoreD.setdefault(base1, dict())
        for base2 in amb.keys():
            if base2 in L:
                scoreD[base1][base2] = match_score
            else:
                scoreD[base1][base2] = mismatch_score
    return (scoreD)


def SWalign(seq1, seq2, pD, useSub=False):
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
    if useSub:
        subs = UtilityFunctions.subMatrixIUPAC(pD['match_score'],
                                               pD['mismatch_score'])
        
        ali = StripedSmithWaterman(seq1, gap_open_penalty=pD['gap_open'],
                                   gap_extend_penalty=pD['gap_extend'],
                                   match_score=pD['match_score'],
                                   mismatch_score=pD['mismatch_score'],
                                   substitution_matrix=subs)(seq2)
    else:
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


def alignFromCIGAR(seq1, seq2, cigar):
    '''
    I = insertion in query
    D = insertion in target
    '''
    if cigar == "%sM" % UtilityFunctions.lengthFromCIGAR(cigar):
        return (seq1, seq2)
    cigar_L = UtilityFunctions.readCIGAR(cigar)
    s1 = ""
    s2 = ""
    n1 = 0
    n2 = 0
    for cig in cigar_L:
        typ = cig[-1]
        count = int(cig[:-1])
        if typ == 'M':
            s1 += seq1[n1:(n1 + count)]
            s2 += seq2[n2:(n2 + count)]
            n1 += count
            n2 += count
        elif typ == 'I' or typ == "Y":
            s1 += seq1[n1:(n1 + count)]
            s2 += '-' * count
            n1 += count
        elif typ == 'D' or typ == "X":
            s1 += '-' * count
            s2 += seq2[n2:n2 + count]
            n2 += count
        else:
            print(typ)
    if (len(s1) != len(s2)):
        raise RuntimeError ("Aligned sequences are not the same length.")
    return (s1, s2)

def getClips(result, query_seq, target_seq):
    Q_length = len(query_seq)
    T_length = len(target_seq)

    Q_start_nclipped = result['query_start']
    T_start_nclipped = result['target_start']

    Q_end_nclipped = Q_length - result['query_end']
    T_end_nclipped = T_length - result['target_end']
    
    return (Q_start_nclipped, T_start_nclipped,
            Q_end_nclipped, T_end_nclipped)
    

def findOverlapType(result, query_seq, target_seq, pD, keep_failed=False):
    '''
    Add various additional annotations to the SW alignment output
    to help determine if the sequences overlap and how much.
    '''
    (Q_start_nclipped, T_start_nclipped,
     Q_end_nclipped, T_end_nclipped) = getClips(result, query_seq,
                                                target_seq)
    
    if pD['clip_perc']:
        max_clip_Q = math.ceil(pD['clip_perc'] * len(query_seq))
        max_clip_T = math.ceil(pD['clip_perc'] * len(target_seq))
    else:
        max_clip_Q, max_clip_T = pD['max_clip'], pD['max_clip']

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


def getAlignmentLocal(result, query_seq, target_seq, pD):
    
    Q_subseq = query_seq[result['query_start']:result['query_end']]
    T_subseq = target_seq[result['target_start']:result['target_end']]

    Q_subseq, T_subseq = alignFromCIGAR(Q_subseq, T_subseq, result['cigar'])
    
    alignment_length = UtilityFunctions.lengthFromCIGAR(result['cigar'])
    Q_subseq = np.array(list(Q_subseq))
    T_subseq = np.array(list(T_subseq))
    
    ident = sum(Q_subseq == T_subseq) / alignment_length
    if ident >= pD['min_perc_ident']:
        return(Q_subseq, T_subseq, ident)
    else:
        return (None, None, None)

 
def getAlignmentFull(result, query_seq, target_seq, pD):
    subseqQ = "".join(list(result['alignment'][0]))
    subseqT = "".join(list(result['alignment'][1]))
    
    (Q_start_nclipped, T_start_nclipped,
     Q_end_nclipped, T_end_nclipped) = getClips(result, query_seq, target_seq)

    Qpad_start = "".join("-" * T_start_nclipped)
    Tpad_start = "".join("-" * Q_start_nclipped)
    
    Qpad_end = "".join("-" * T_end_nclipped)
    Tpad_end = "".join("-" * Q_end_nclipped)

    if Q_start_nclipped == 0 or T_start_nclipped == 0:
        start_cigar_start = ""
    else:
        if Q_start_nclipped < T_start_nclipped:
            start_cigar_start = "%iY" % (Q_start_nclipped)
        else:
            start_cigar_start = "%iX" % (T_start_nclipped)

    if T_start_nclipped >= Q_start_nclipped:
        if T_start_nclipped == 0:
            start_cigar = "%s" % start_cigar_start
        else:
            start_cigar = "%s%iD" % (start_cigar_start, T_start_nclipped)
    else:
        start_cigar = "%s%iI" % (start_cigar_start, Q_start_nclipped)
        
    if Q_end_nclipped == 0 or T_end_nclipped == 0:
        end_cigar_end = ""
    else:
        if Q_end_nclipped < T_end_nclipped:
            end_cigar_end = "%iY" % (Q_end_nclipped)
        else:
            end_cigar_end = "%iX" % (T_end_nclipped)
    if T_end_nclipped >= Q_end_nclipped:
        if T_end_nclipped == 0:
            end_cigar = "%s" % end_cigar_end
        else:
            end_cigar = "%iD%s" % (T_end_nclipped, end_cigar_end)
    else:
        end_cigar = "%iI%s" % (Q_end_nclipped, end_cigar_end)
    cigar_updated = "%s%s%s" % (start_cigar, result['cigar'], end_cigar)
    fullseqQ = Qpad_start + query_seq[:result['query_start']] + subseqQ + query_seq[result['query_end']:] + Qpad_end
    fullseqT = Tpad_start + target_seq[:result['target_start']] + subseqT + target_seq[result['target_end']:] + Tpad_end
    result['query_seq_aligned'] = fullseqQ
    result['target_seq_aligned'] = fullseqT
    if len(fullseqT) != len(fullseqQ):
        raise RuntimeError ("Aligned sequences are not the same length")
    assert UtilityFunctions.lengthFromCIGAR(cigar_updated) == len(fullseqQ), "CIGAR length is not consistent with sequence length"
    result['cigar_updated'] = cigar_updated

    return(result)