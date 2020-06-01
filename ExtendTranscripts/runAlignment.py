#!/usr/bin/env python3
import AlignmentPW
import UtilityFunctions
import copy
import Consensus
import numpy as np


def runAlignment(fasta_dict, pD, alignment_type='pairwise',
                 quick=False, rename=False):
    '''

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
    seqdict = dict()
    
    Z = sorted(zip(nams, seqs), key=lambda x: len(x[1]))[::-1]
    
    if alignment_type == 'pairwise':
        D = runAlignmentPW(Z, namD, fasta_dict, pD,
                           quick=quick,
                           rename=rename)
    elif alignment_type == 'stepwise':
        D = runAlignmentSW(Z,
                           namD,
                           fasta_dict,
                           pD,
                           seqdict,
                           rename=rename)
    if rename:
        return (namD, rev_namD, D)
    else:
        return (D)

def alignmentMeetsCriteria(result, query_seq, target_seq, pD):
    '''
    Check that the alignment meets the criteria:
        * alignment length > pD['min_length']
        * either one sequence overlaps at one or both ends, the whole of one
          sequence is contained within the other or the sequences are the same
          length
        * alignment identity > pD['min_perc_ident']

    Each step is only run if the previous criteria is passed, they are ordered
    by how long they take.
    '''
    # check the length of the aligned region
    length = UtilityFunctions.lengthFromCIGAR(result['cigar'])
    # this is quickest - if it's too short dont' check anything
    # else
    if length >= pD['min_length']:
        wo_indels = UtilityFunctions.lengthFromCIGAR(result['cigar'],
                                                     mOnly=True)
        indels = wo_indels / length
        if indels > (1 - pD['max_indels']):
            # check how the sequences overlap - the aligned region
            # needs to either include the end of a sequence or 
            
            SE = AlignmentPW.findOverlapType(result, query_seq, target_seq,
                                             pD, keep_failed=False)
            if SE:
                alignment = AlignmentPW.getAlignmentLocal(result,
                                                         query_seq,
                                                         target_seq, pD)
                alignment_identity = alignment[2]
                if alignment_identity:
                    return (True, alignment)
                else:
                    # all the False values are one element tuples so that
                    # the alignment can be returned and accessed if the
                    # result is True
                    return (False, )
            else:
                return (False, )
        else:
            return (False, )
    else:
        return (False, )
                
    
def runAlignmentSW(Z, namD, fasta_dict, pD, seqdict, rename=False):
    '''
    '''
    done = set([Z[0][0]])
    Z = list(Z)
    j = 0
    # X holds the current set of unassigned sequences - sequences are removed
    # from X as they are assigned to a group
    X = copy.copy(Z)
    D = dict()
    prevlen = len(X) + 1
    # keep running until there are no sequences left
    while len(done) < len(Z) and len(X) > 1 and len(X) < prevlen:
        print (j, len(X), prevlen)
        prevlen = len(X)
        # the query sequence is the current consensus - start off with the longest
        # sequence in the input
        query_seq = X[0][1]
        query_nam = X[0][0]
        if "consensus" not in query_nam:
            p_ali = UtilityFunctions.AlignmentArray([query_seq])
            iscons = False
        else:
            consensus_n = int(X[0].replace("consensus_", ""))
            p_ali = D[consensus_n]['alignment']
            iscons = True
        matrix, nt_inds = Consensus.makeAlignmentMatrix(p_ali)
        X = X[1:]
        i = 0
        groupdone = list([X[0][0]])
        done = done | set(groupdone)
        level = 0
        consensusD = dict()
        # keep starting at the top again with the new consensus sequence
        while i != len(X):
            # Run through all the remaining sequences not in the consensus
            while True:
                
                if i >= len(X):
                    break
                # the target sequence is the next longest sequence in the input
                target_seq = X[i][1]
                
                result = AlignmentPW.SWalign(query_seq, target_seq,
                                             pD, useSub=False)
                is_match = alignmentMeetsCriteria(result, query_seq,
                                                  target_seq, pD)
                if not is_match[0]:
                    target_seq = UtilityFunctions.reverseComplement(target_seq)
                    result = AlignmentPW.SWalign(query_seq, target_seq,
                                                 pD, useSub=False)
                    is_match = alignmentMeetsCriteria(result, query_seq,
                                                      target_seq, pD)                    
                if is_match[0]:
                    result['alignment'] = is_match[1]
                    result = AlignmentPW.getAlignmentFull(result,
                                                          query_seq,
                                                          target_seq,
                                                          pD)
                    # keep track of which sequences are in this group
                    groupdone.append(X[i][0])
                    # keep track of which sequence are in any group
                    done.add(X[i][0])
                    # this sequence has been assigned - remove it from X
                    X = X[:i] + X[i+1:]
                    ali, matrix = Consensus.expandAlignment(result,
                                                            matrix,
                                                            nt_inds,
                                                            consensusD,
                                                            level)
                    consensus = Consensus.collapseAlignment(matrix, nt_inds)
                    consensusD[level] = dict()
                    consensusD[level]['alignment'] = ali
                    consensusD[level]['consensus'] = consensus
                    if level - 2 in consensusD:
                        del consensusD[level-2]
                    query_seq = consensusD[level]['consensus']
                    level += 1
                    i = 0
                    break
                i += 1
        D[j] = dict()
        if level != 0:
            ali, matrix, consensus, seqdict = AlignmentPW.cleanAlignmentCIAlign(ali,
                                                                                groupdone,
                                                                                consensus,
                                                                                matrix,
                                                                                nt_inds,
                                                                                seqdict)
            D[j]['consensus'] = consensus
            D[j]['alignment'] = ali
            D[j]['names'] = groupdone
            D[j]['log'] = seqdict
            if len(X) != 0 and (iscons is False or
                  np.shape(ali) != np.shape(p_ali)):
                    X.append(("%s_consensus" % j,
                              consensusD[level-1]['consensus']))
        else:
            D[j]['consensus'] = query_seq
            D[j]['alignment'] = None
            D[j]['names'] = None
            D[j]['log'] = None
        j += 1
    if len(X) == 1:
        D[j] = dict()
        D[j]['consensus'] = X[0][1]
        D[j]['alignment'] = None
        D[j]['names'] = None
        D[j]['log'] = None
    return (D)
        
    
def runAlignmentPW(Z, namD, fasta_dict, pD, quick=False,
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
            is_match = alignmentMeetsCriteria(result, query_seq,
                                  target_seq, pD)
            if not is_match[0]:
                target_seq = UtilityFunctions.reverseComplement(target_seq)
                result = AlignmentPW.SWalign(query_seq, target_seq,
                                             pD, useSub=False)
                is_match = alignmentMeetsCriteria(result, query_seq,
                                                  target_seq, pD)
            if is_match[0]:
                result['alignment'] = is_match[1]
                
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
    return (D)