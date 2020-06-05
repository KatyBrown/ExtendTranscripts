#!/usr/bin/env python3
import Alignment
import UtilityFunctions
import copy
import Consensus
from UtilityFunctions import logPrint as lp


def runAlignment(fasta_dict, pD, outdir, alignment_type='pairwise',
                 quick=False):

    nams = fasta_dict.keys()
    seqs = fasta_dict.values()
    seqdict = {nam: dict() for nam in nams}
    lp("Reading input file %s" % pD['infile'], 1, pD)
    lp("%i sequences in input files" % (len(nams)), 2, pD)
    for nam in nams:
        seqdict[nam]['is_rc'] = False

    Z = sorted(zip(nams, seqs), key=lambda x: len(x[1]))[::-1]

    if alignment_type == 'pairwise':
        D = runAlignmentPW(Z, fasta_dict, pD,
                           quick=quick)
    elif alignment_type == 'stepwise':
        D = runAlignmentSW(Z,
                           fasta_dict,
                           pD,
                           seqdict)
    writeFastas(D, outdir)
    return (D)


def runAlignmentSW(Z, fasta_dict, pD, seqdict):
    '''
    '''
    UtilityFunctions.logPrint("Starting step-wise alignment.", 1, pD)

    # Z is a list of tuples arranged as (nam, sequence) representing
    # all the sequences in the input file
    # make a copy of Z to manipulate
    # X holds the current set of unassigned sequences - sequences are removed
    # from X as they are assigned to a group
    X = copy.copy(Z)
    done = set()
    # Three conditions must be True to run again
    # Not all the input sequences are in done
    # There is at least one sequence left in X

    # dictionary to store the results
    D = dict()
    k = 0
    while len(done) < len(Z) and len(X) > 1:
        # We will make this True if something changes
        query_nam, current_consensus = X[0]
        lp("Remaining unaligned sequences %i" % (len(X)), 2, pD)
        lp("Starting new cluster with %s" % (query_nam), 2, pD)

        # if the query is a previous consensus sequence
        if "*consensus" in query_nam:
            query_cons_ind = int(query_nam.split("_")[0])
            current_alignment = D[query_cons_ind]['alignment']
            current_names = D[query_cons_ind]['names']
        else:
            current_alignment = UtilityFunctions.AlignmentArray(
                    [current_consensus])
            current_names = [query_nam]
        # convert the query alignment into a matrix
        current_mat, nt_inds = Consensus.makeAlignmentMatrix(
                current_alignment)
        # remove the query sequence from X
        X = X[1:]
        # keep track of position in X
        j = 0
        # while we are not at the bottom of x
        while j != len(X) and len(X) != 0:
            i = 0
            any_matches_inner = False
            n_new = 0
            # until we find a match and change the query sequence
            while True and i != len(X):
                query_seq = current_consensus
                query_ali = current_alignment
                query_names = current_names
                target_nam, target_seq = X[i]

                if "*consensus" in target_nam:
                    target_cons_ind = int(target_nam.split("_")[0])
                    target_ali = D[target_cons_ind]['alignment']
                    target_names = D[target_cons_ind]['names']
                else:
                    target_ali = UtilityFunctions.AlignmentArray([target_seq])
                    target_names = [target_nam]

                lp("Testing %s" % ", ".join(target_names), 3, pD)
                # Align the query and target consensus sequences
                result = Alignment.SWalign(query_seq, target_seq,
                                           pD, useSub=False)

                # Check the if the consensus sequences are a good match
                is_match = Alignment.alignmentMeetsCriteria(result, query_seq,
                                                            target_seq, pD)
                # if they are not try the reverse complement

                if not is_match[0]:
                    target_seq = UtilityFunctions.reverseComplement(target_seq)
                    result = Alignment.SWalign(query_seq, target_seq,
                                               pD, useSub=False)
                    is_match = Alignment.alignmentMeetsCriteria(result,
                                                                query_seq,
                                                                target_seq,
                                                                pD)
                    ####################
                    UtilityFunctions.reverseComplementAlignment(target_ali)
                    target_ali = UtilityFunctions.AlignmentArray([target_seq])
                    for nam in target_names:
                        seqdict[nam]['is_rc'] = True

                if is_match[0]:
                    lp("Match found.", 3, pD)
                    # We found a match - something has changed
                    any_matches_inner = True
                    n_new += 1
                    # remove the current value from X
                    X = X[:i] + X[i+1:]
                    result['alignment'] = is_match[1]
                    # get the full alignment for the two consensus sequences
                    result = Alignment.getAlignmentFull(result,
                                                        query_seq,
                                                        target_seq,
                                                        pD)
                    current_names = query_names + target_names

                    lp("Expanding current alignment to include %s" % (
                            ", ".join(target_names)), 3, pD)
                    ali, current_mat = Consensus.expandAlignment(result,
                                                                 query_ali,
                                                                 target_ali,
                                                                 current_mat,
                                                                 nt_inds)
                    # make a new sequence based on the new alignment
                    current_consensus = Consensus.collapseAlignment(
                            current_mat, nt_inds)
                    current_alignment = ali
                    i = 0
                    # now you have a match and the consensus is updated,
                    # start at the top again
                    break
                else:
                    lp("No match.", 3, pD)
                # keep going through the other sequences
                i += 1
            j += 1
            done = done | set(current_names)
            if any_matches_inner:
                # if anything has changed, clean up the alignment etc
                lp("Cluster %i updated - %s sequences" % (k,
                                                          len(current_names)),
                   2, pD)
                lp("Cleaning cluster %i with CIAlign" % (k), 3, pD)
                R = Alignment.cleanAlignmentCIAlign(ali,
                                                    current_names,
                                                    query_seq,
                                                    current_mat,
                                                    nt_inds,
                                                    seqdict)
                current_alignment, current_mat, current_consensus, seqdict = R
            else:
                break
        D.setdefault(k, dict())
        lp("Saved cluster %i with %s sequences" % (k, len(current_names)),
           1, pD)
        D[k]['consensus'] = current_consensus
        D[k]['alignment'] = current_alignment
        D[k]['names'] = current_names
        D[k]['log'] = seqdict
        k += 1
    if len(X) == 1:
        D.setdefault(k, dict())
        lp("Saved cluster %i with %s sequences" % (k, 1),
           1, pD)
        D[k]['consensus'] = X[0][1]
        D[k]['alignment'] = UtilityFunctions.AlignmentArray([X[0][1]])
        D[k]['names'] = [X[0][0]]
        D[k]['log'] = dict()
    return (D)


def runAlignmentPW(Z, fasta_dict, pD, quick=False):
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

            result = Alignment.SWalign(query_seq, target_seq,
                                       pD)
            is_match = Alignment.alignmentMeetsCriteria(result, query_seq,
                                                        target_seq, pD)
            if not is_match[0]:
                target_seq = UtilityFunctions.reverseComplement(target_seq)
                result = Alignment.SWalign(query_seq, target_seq,
                                           pD, useSub=False)
                is_match = Alignment.alignmentMeetsCriteria(result, query_seq,
                                                            target_seq, pD)
            if is_match[0]:
                result['alignment'] = is_match[1]
                D.setdefault(query_nam, dict())
                D[query_nam][target_nam] = result

                if quick:
                    break
                else:
                    D.setdefault(target_nam, dict())
                    D[target_nam][query_nam] = result
    return (D)


def writeFastas(results, outdir):
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
    for i in results.keys():
        ali_out = open("%s/cluster_%s_alignment.fasta" % (outdir, i), "w")
        cons_out = open("%s/cluster_%s_consensus.fasta" % (outdir, i), "w")
        both_out = open("%s/cluster_%s_ali_plus_cons.fasta" % (outdir, i), "w")

        ali = results[i]['alignment']
        cons = results[i]['consensus']
        names = results[i]['names']

        for i, nam in enumerate(names):
            ali_out.write(">%s\n%s\n" % (nam, "".join(list(ali[i]))))
            both_out.write(">%s\n%s\n" % (nam, "".join(list(ali[i]))))
        cons_out.write(">consensus_%i\n%s\n" % (i, cons))
        both_out.write(">consensus_%i\n%s\n" % (i, cons))
        ali_out.close()
        cons_out.close()
        both_out.close()
