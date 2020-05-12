#!/usr/bin/env python3
import numpy as np
import UtilityFunctions

def matchFromCIGAR(old_alignment_array,
                   new_alignment_array, cigar,
                   cigar_updated):
    '''
    I = insertion in query
    D = insertion in target
    '''
    if UtilityFunctions.lengthFromCIGAR(cigar_updated) != np.shape(new_alignment_array)[1]:
        raise RuntimeError ("Something is wrong with the cigar string")
    if cigar_updated == "%sM" % UtilityFunctions.lengthFromCIGAR(cigar_updated):
        full_alignment_array = np.concatenate((old_alignment_array,
                                               new_alignment_array))

        
    else:
        cigar_L = UtilityFunctions.readCIGAR(cigar_updated)
        if "D" not in cigar_updated and "X" not in cigar_updated:
            full_alignment_array = np.concatenate((old_alignment_array,
                                                  new_alignment_array))
        else:
            n = 0
            old_alignment_array_updated = old_alignment_array
            for cig in cigar_L:
                typ = cig[-1]
                count = int(cig[:-1])
                for i in range(count):
                    if typ == "D" or typ == "X":
                        newA_left = old_alignment_array_updated[:,:n]
                        newA_right = old_alignment_array_updated[:,n:]
                        x = np.array([["-"] * np.shape(old_alignment_array_updated)[0]]).T
                        old_alignment_array_updated = np.concatenate((
                                newA_left, x, newA_right), axis=1)
                    n += 1

            full_alignment_array = np.concatenate((
                    old_alignment_array_updated, new_alignment_array))
    return (full_alignment_array)


def collapseAlignment(A, consensusD, level):
    D1, D2 = UtilityFunctions.IUPAC()
    C = []
    for i in range(np.shape(A)[1]):
        this_column = A[:, i]
        this_column = this_column[this_column != "-"]
        unique, counts = np.unique(this_column, return_counts=True)
        if len(unique) == 0:
            # the column only has one nt
            C.append("-")
        elif len(unique) == 1:
            # the column is all gaps
            C.append(unique[0])
        else:
            # the column has >1 non-gap character
            maxi = unique[counts == max(counts)]
            maxi_string = "".join(list(maxi))
            if maxi_string not in D2:
                # there is one majority nt
                C.append(maxi_string)
            else:
                # there is >1 - use an IUPAC character
                iupac = D2[maxi_string]
                C.append(iupac)
    # collapse the consensus into a string
    C = "".join(C)
    consensusD[level] = dict()
    consensusD[level]['alignment'] = A
    consensusD[level]['consensus'] = C
    return (consensusD)


def expandAlignment(result, consensusD, level):
    query_seq = result['query_seq_aligned']
    target_seq = result['target_seq_aligned']
    cigar = result['cigar']
    cigar_updated = result['cigar_updated']
    if level == 0:
        full_alignment_array = UtilityFunctions.AlignmentArray([query_seq,
                                                                target_seq])
    else:
        old_alignment_array = consensusD[level - 1]['alignment']
        new_alignment_array = UtilityFunctions.AlignmentArray([target_seq])
        full_alignment_array= matchFromCIGAR(old_alignment_array,
                                             new_alignment_array,
                                             cigar,
                                             cigar_updated)
    
    return (full_alignment_array)
    
