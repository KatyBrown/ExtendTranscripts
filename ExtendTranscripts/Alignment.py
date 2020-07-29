#!/usr/bin/env python3
from skbio.alignment import StripedSmithWaterman
import numpy as np
import math
import UtilityFunctions
import Consensus
import sys
import os
from UtilityFunctions import logPrint as lp
# temporary until CIAlign is updated
sys.path.insert(0, "/home/katy/CIAlign_P")
import CIAlign.parsingFunctions
    
    
def writeFastas(result, k, rround, outdir, candidates=False,
                reference=None):
    '''
    Output three fasta files from a results dictionary generated with
    runAlignment.runAlignment.
    The dictionary should include:
        ali - a np.array containing the alignment
        consensus - the consensus sequence

    sequence array (where each row is an array
    consisting of an aligned nucleotide sequecne), a list of names

    Parameters
    ----------
    fasta_dict: dict
        Dictionary where keys are sequence names and values are nucleotide
        sequences
    outfile: str
        path to output file
    '''
    if not os.path.exists("%s/round_%i" % (outdir, rround)):
        os.mkdir("%s/round_%i" % (outdir, rround))
    ali_out = open("%s/round_%i/cluster_%s_alignment.fasta" % (
        outdir, rround, k), "w")
    cons_out = open("%s/round_%i/cluster_%s_consensus.fasta" % (
        outdir, rround, k), "w")
    both_out = open("%s/round_%i/cluster_%s_ali_plus_cons.fasta" % (
        outdir, rround, k), "w")

    ali = result['alignment']
    cons = result['consensus']
    names = result['names']
    for i, nam in enumerate(names):
        ali_out.write(">%s\n%s\n" % (nam, "".join(list(ali[i]))))
        both_out.write(">%s\n%s\n" % (nam, "".join(list(ali[i]))))
    cons_out.write(">*consensus_%i\n%s\n" % (k, cons))
    both_out.write(">*consensus_%i\n%s\n" % (k, cons))
    ali_out.close()
    cons_out.close()
    both_out.close()


def alignmentMeetsCriteria(result, query_seq, target_seq, pD,
                           is_candidate=False):
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
    qlength = result['query_end'] - result['query_start']
    tlength = result['target_end'] - result['target_start']
    # this is quickest - if it's too short dont' check anything
    if pD['alignment_length_min_n']:
        min_length_Q = pD['alignment_length_min_n']
        min_length_T = pD['alignment_length_min_n']
    else:
        min_length_Q = math.ceil(
            pD['alignment_length_min_perc'] * len(query_seq))
        min_length_T = math.ceil(
            pD['alignment_length_min_perc'] * len(target_seq))
    # else
    if qlength >= min_length_T and tlength >= min_length_Q:
        wo_indels = UtilityFunctions.lengthFromCIGAR(result['cigar'],
                                                     mOnly=True)
        indels = wo_indels / length
        if pD['alignment_indels_max_n']:
            max_indels = pD['alignment_indels_max_n']
        else:
            max_indels = math.ceil(pD['alignment_indels_max_perc'] * length)
        if indels > (1 - max_indels):
            # check how the sequences overlap - the aligned region
            # needs to either include the end of a sequence or neither
            # sequence is heavily clipped

            SE = findOverlapType(result, query_seq, target_seq,
                                 pD, keep_failed=False)
            if SE:
                alignment = getAlignmentLocal(result,
                                              query_seq,
                                              target_seq,
                                              pD)
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


def gappedToUngapped(row, removed_pos):
    '''
    Takes an aligned sequence (with gaps) and a set of indices of residues
    which have been removed from this sequence, and returns the indicies of
    the residues within the same sequence but
    unaligned (without gaps).  Indices which correspond to gaps in the input
    sequence are excluded.
    Returns a list of tuples. In each of these tuples, tuple[0] is a list of
    adjacent positions (in the unaligned original sequence) which
    have been removed and tuple[1] is a string containing the part
    of the sequence which has been removed.

    Parameters
    ----------
    row: np.array
        numpy array containing the aligned sequence
    removed_pos: np.array
        numpy array containing the positions which were removed from the
        aligned sequence

    Returns
    -------
    result: list
        list of tuples showing which positions have been removed from the
        unaligned sequence and which residues were at these positions.

    '''
    # Make a True / False array of non-gap positions
    non_gap_pos = row != "-"
    # Make an np.array of row with gaps removed
    row_ng = row[non_gap_pos]
    result = []
    # We are only interested when non-gap residues have been removed
    if sum(non_gap_pos[removed_pos]) != 0:
        # find the indices in the gap-free sequence of the
        # removed residues
        # np.cumsum counts the number of non-gap residues before each position
        non_gap_ind = np.cumsum(non_gap_pos)
        # find how many non-gap residues precede each removed position - this
        # represents it's index in the unaligned sequence
        this_row_ind = non_gap_ind[removed_pos]
        # remove duplicates
        this_row_ind = sorted(np.unique(this_row_ind))
        current_L = [0]
        prev_r = this_row_ind[0]
        i = 1
        for r in this_row_ind[1:]:
            if r - prev_r == 1:
                current_L.append(i)
            else:
                substring = "".join(row_ng[current_L[0]:current_L[-1] + 1])
                result.append((current_L, substring))
                current_L = [i]
            prev_r = r
            i += 1
        substring = "".join(row_ng[current_L[0]:current_L[-1] + 1])
        result.append((current_L, substring))
    return (result)


def updateSeqDict(orig_arr, nams, logD, seqdict):
    '''
    Records which parts of the original sequences have been clipped out
    using CIAlign when refining a consensus sequence.

    Parameters
    ----------
    orig_arr: np.array
        The original array on which CIAlign was run
    nams: list
        A list of sequence names in the same order as orig_arr
    logD: dict
        A dictionary containing the CIAlign log for each function - for
        crop_ends this is another dictionary where keys are sequence names
        and vals are the indicies of residues which have been replaced with
        gaps, for remove_insertions and remove_gaponly these are sets
        of the indicies of columns which have been removed.  These are all
        relative positions - the position in the original alignment (orig_arr),
        rather than the position at this point in processing.
    seqdict: dict
        A dictionary with any processing previously perfomed on any sequences.

    Returns
    -------
    seqdict: dict
        An updated dictionary where keys are sequence IDs and values are
        dictionaries, within these dictionaries keys are CIAlign function
        names and values are lists of tuples.
        In each of these tuples, tuple[0] is a list of adjacent positions
        (in the unaligned original sequence) which have been removed using
        this CIAlign function and tuple[1] is a string containing the part
        of the sequence which has been removed.
    '''
    # These steps need to be performed in the same order they were originally
    # carried out, as the positions in the final sequences will change after
    # each step.
    for function in logD['order']:
        # Check if anything was removed by this function in this alignment
        if len(logD[function]) != 0:
            # log for remove_insertions function
            if function == 'remove_insertions':
                # sort the list of positions which were removed
                removed_pos = sorted(np.array(list(logD[function])))
                for i, row in enumerate(orig_arr):
                    nam = nams[i]
                    # find  the location of the removed residues in the
                    # ungapped sequence and the piece of sequence
                    # which has been removed
                    result = gappedToUngapped(row, removed_pos)
                    if len(result) != 0:
                        # store this result in seqdict
                        seqdict.setdefault(nam, dict())
                        seqdict[nam][function] = result

            elif function == 'crop_ends':
                cropped = logD[function]
                for seqnam in cropped:
                    seqdict.setdefault(seqnam, dict())
                    row_ind = nams.index(seqnam)
                    row = orig_arr[row_ind, :]
                    removed_pos = np.concatenate(cropped[seqnam])
                    result = gappedToUngapped(row, removed_pos)
                    seqdict[seqnam][function] = result
    return(seqdict)


def cleanAlignmentCIAlign(arr,
                          nams,
                          consensus, matrix, nt_inds, seqdict, pD,
                          functions=['remove_insertions',
                                     'crop_ends',
                                     'remove_gaponly']):

    # make a list of integers - one for each postiion
    # this is used to keep track of which positions from the input alignment
    # have been removed - as after columns and rows are removed the indices
    # change
    relativePositions = list(range(0, len(arr[0])))
    # dictionary to keep track of what has been removed
    logD = dict()
    
    # make a copy of the input array
    orig_arr = arr
    # these will always be run in the same order for now but might as
    # well make it possible to change the order in case it's needed later

    if pD['alignment_clip_max_n']:
        clipmax = pD['alignment_clip_max_n'] / len(consensus)
    else:
        clipmax = pD['alignment_clip_max_perc']

    logD['order'] = functions
    for function in functions:
        # check there's no weird functions in there
        assert function in ['remove_insertions',
                            'crop_ends',
                            'remove_gaponly'], (
                                "CIAlign function %s not found" % function)
        if function == "remove_insertions" and pD[
                'remove_insertions_maxlen'] != 0:
            lp("Removing insertions", 3, pD)
            # store the previous set of indices
            p_relative = np.array(relativePositions)
            # run the CIAlign remove insertions function
            # returns a cleaned array and a  list of removed positions
            # only small insertions are removed - there shouldn't be big
            # ones using this alignment method
            R = CIAlign.parsingFunctions.removeInsertions(arr,
                                                          relativePositions,
                                                          logtype='dict',
                                                          min_size=1,
                                                          max_size=30)
            arr, r, relativePositions = R
            # convert the removed positions to an np.array so they can be
            # used as an index
            removed_relative = np.array(list(r))
            # find out which postions in the original alignment these
            # correspond to
            keep_absolute = np.where(
                np.invert(np.in1d(p_relative, removed_relative)))[0]
            matrix = matrix[:, keep_absolute]

        elif function == "crop_ends" and clipmax != 0:
            lp("Cropping ends", 3, pD)

            # store the previous set of indices
            p_relative = np.array(relativePositions)
            # store the input array
            p_arr = arr
            # run the CIAlign crop ends function
            # returns a cleaned array and a dictionary where keys are
            # sequence IDs and vals are tuples - tuple[0] is the positions
            # removed (replaced with "-") at the beginning
            # and tuple[1] is the positions removed from the end
            arr, r = CIAlign.parsingFunctions.cropEnds(arr, nams,
                                                       relativePositions,
                                                       redefine_perc=clipmax,
                                                       mingap=0.001,
                                                       logtype='dict')
            # iterate through the dictionary
            for nam in r:
                cs, ce = r[nam]
                rm_absolute_cs = np.where(np.in1d(p_relative, cs))[0]
                rm_absolute_ce = np.where(np.in1d(p_relative, ce))[0]
                rm_absolute = np.concatenate((rm_absolute_cs, rm_absolute_ce))
                row_index = nams.index(nam)
                this_row = p_arr[row_index, :]
                rm_nucs = this_row[rm_absolute]
                nuc_inds = np.array([nt_inds[x] for x in rm_nucs])
                matrix[nuc_inds, rm_absolute] -= 1

        elif function == "remove_gaponly":
            p_relative = np.array(relativePositions)
            R = CIAlign.parsingFunctions.removeGapOnly(arr,
                                                       relativePositions,
                                                       logtype='dict')
            arr, r, relativePositions = R
            removed_relative = np.array(list(r))
            keep_absolute = np.where(
                np.invert(np.in1d(p_relative, removed_relative)))[0]
            matrix = matrix[:, keep_absolute]
        logD[function] = r
    consensus = Consensus.collapseAlignment(matrix, nt_inds)
    seqdict = updateSeqDict(orig_arr, nams, logD, seqdict)
    return (arr, matrix, consensus, seqdict)


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
    amb = UtilityFunctions.IUPAC()[0]
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
    gap_open = pD['alignment_gap_open']
    gap_extend = pD['alignment_gap_extend']
    match_score = pD['alignment_match_score']
    mismatch_score = pD['alignment_mismatch_score']

    if useSub:
        subs = subMatrixIUPAC(pD['alignment_match_score'],
                              pD['alignment_mismatch_score'])

        ali = StripedSmithWaterman(seq1, gap_open_penalty=gap_open,
                                   gap_extend_penalty=gap_extend,
                                   match_score=match_score,
                                   mismatch_score=mismatch_score,
                                   substitution_matrix=subs)(seq2)
    else:
        ali = StripedSmithWaterman(seq1, gap_open_penalty=gap_open,
                                   gap_extend_penalty=gap_extend,
                                   match_score=match_score,
                                   mismatch_score=mismatch_score)(seq2)

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
        raise RuntimeError("Aligned sequences are not the same length.")
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

    if pD['alignment_clip_max_n']:
        max_clip_Q = pD['alignment_clip_max_n']
        max_clip_T = pD['alignment_clip_max_n']
    else:
        max_clip_Q = math.ceil(pD['alignment_clip_max_perc'] * len(query_seq))
        max_clip_T = math.ceil(pD['alignment_clip_max_perc'] * len(target_seq))

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
        raise RuntimeError("Unexpected overlap at sequence start.")

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
        raise RuntimeError("Unexpected overlap at sequence end.")
    return (start_type, end_type)


def getAlignmentLocal(result, query_seq, target_seq, pD):
    Q_subseq = query_seq[result['query_start']:result['query_end']]
    T_subseq = target_seq[result['target_start']:result['target_end']]

    Q_subseq, T_subseq = alignFromCIGAR(Q_subseq, T_subseq, result['cigar'])

    alignment_length = UtilityFunctions.lengthFromCIGAR(result['cigar'])
    Q_subseq = np.array(list(Q_subseq))
    T_subseq = np.array(list(T_subseq))

    ident = sum(Q_subseq == T_subseq) / alignment_length
    if 'alignment_ident_min_n' in pD:
        if pD['alignment_ident_min_n']:
            min_ident = pD['alignment_ident_min_n'] / alignment_length
        else:
            min_ident = pD['alignment_ident_min_perc']
        if ident >= min_ident:
            return(Q_subseq, T_subseq, ident)
        else:
            return (None, None, None)
    else:
        return (Q_subseq, T_subseq, ident)


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
        # if one end or the other has nothing clipped
        start_cigar_start = ""
    else:
        # if both are clipped
        if Q_start_nclipped < T_start_nclipped:
            start_cigar_start = "%iY" % (Q_start_nclipped)
        else:
            start_cigar_start = "%iX" % (T_start_nclipped)
    # start_cigar = start_cigar_start
    if T_start_nclipped > Q_start_nclipped:

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
    if T_end_nclipped > Q_end_nclipped:
        if T_end_nclipped == 0:
            end_cigar = "%s" % end_cigar_end
        else:
            end_cigar = "%iD%s" % (T_end_nclipped, end_cigar_end)
    else:
        end_cigar = "%iI%s" % (Q_end_nclipped, end_cigar_end)
    cigar_updated = "%s%s%s" % (start_cigar, result['cigar'], end_cigar)

    qstart = query_seq[:result['query_start']]
    qend = query_seq[result['query_end']:]

    tstart = target_seq[:result['target_start']]
    tend = target_seq[result['target_end']:]

    fullseqQ = Qpad_start + qstart + subseqQ + qend + Qpad_end
    fullseqT = Tpad_start + tstart + subseqT + tend + Tpad_end

    result['query_seq_aligned'] = fullseqQ
    result['target_seq_aligned'] = fullseqT
    if len(fullseqT) != len(fullseqQ):
        raise RuntimeError("Aligned sequences are not the same length")
    assert UtilityFunctions.lengthFromCIGAR(cigar_updated) == len(fullseqQ), (
        "CIGAR length is not consistent with sequence length")
    result['cigar_updated'] = cigar_updated

    return(result)


def mergeFastas(rround, nfastas, outdir):
    F = dict()
    for i in range(1, nfastas+1):
        F.update(
            UtilityFunctions.FastaToDict(
                "%s/round_%i/cluster_%i_consensus.fasta" % (
                    outdir, rround, i)))
    nams = sorted(list(F.keys()), key=lambda x: int(x.split("_")[1]))
    out = open("%s/round_%i/consensus_combined.fasta" % (outdir, rround), "w")
    for nam in nams:
        out.write(">%s\n%s\n" % (nam, F[nam]))
    out.close()
