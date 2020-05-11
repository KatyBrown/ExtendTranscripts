#!/usr/bin/env python3
import AlignmentPW

def runAlignmentPW(fasta_dict, pD, quick=False, keep_failed=False,
                   rename=False):
    '''
    Runs SWalign on every pair of sequences in a Fasta file converted to
    a dictionary.
    
    Parameters
    ----------
    fasta_dict: dict
        dictionary with sequence names as keys and sequences (strings) as
        values

    min_perc_ident: float
        minimum proportion identity to expect within overlapping regions
    max_clip: int
        maximum number of bases to clip from the ends of sequences
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
    D: dict
        Dictionary of dictionaries containing the SWalign results for every
        pair of sequences in the input dictionary.
        Each dictionary key is a sequence and the value is another dictionary
        with the name of another sequence as the key and the alignment
        dictionary for the two sequences as the value.
    
    '''
    D = dict()
    if rename:
        namD = {nam:i for i, nam in enumerate(fasta_dict.keys())}
        rev_namD = {i:nam for i, nam in enumerate(fasta_dict.keys())}
    else:
        namD = None
        rev_namD = None

    nams = fasta_dict.keys()
    seqs = fasta_dict.values()
    Z = sorted(zip(nams, seqs), key=lambda x: len(x[1]))[::-1]
    
    for i, (N1, S1) in enumerate(Z):
        for j, (N2, S2) in enumerate(Z[i+1:]):
            if len(S2) > len(S1):
                query_nam = N2
                query_seq = S2
                target_nam = N1
                target_seq = S1
            else:
                query_nam = N1
                query_seq = S1
                target_nam = N2
                target_seq = S2
    
            result = AlignmentPW.SWalign(query_seq, target_seq,
                                         pD)
            length = AlignmentPW.lengthFromCIGAR(result['cigar'])
            if length >= pD['min_length'] or keep_failed:
                SE = AlignmentPW.findOverlapType(result, query_seq, target_seq,
                                                 pD, keep_failed=keep_failed)
                if SE:
                    alignment = AlignmentPW.getAlignment(result,
                                                         query_seq,
                                                         target_seq, pD,
                                                         keep_failed=keep_failed)
                    alignment_identity = alignment[2]
                    if alignment_identity or keep_failed:
                        result['SE'] = SE
                        result['alignment'] = alignment
                        
                        if rename:
                            D.setdefault(namD[query_nam], dict())
                            D[namD[query_nam]][namD[target_nam]] = result
                        else:
                            D.setdefault(query_nam, dict())
                            D[query_nam][target_nam] = result
                        if quick:
                            break
                        else:  
                            if rename:
                                D.setdefault(namD[target_nam], dict())
                                D[namD[target_nam]][namD[query_nam]] = result
                            else:
                                D.setdefault(target_nam, dict())
                                D[target_nam][query_nam] = result
    if rename:
        return (namD, rev_namD, D)
    else:
        return (D)