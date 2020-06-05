#!/usr/bin/env python3
import Alignment
import UtilityFunctions
import copy
import Consensus
from UtilityFunctions import logPrint as lp
import os


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
        current = dict()
        current['name'], current['consensus'] = X[0]
        current['names'] = [current['name']]
        lp("Remaining unaligned sequences %i" % (len(X)), 2, pD)

        # convert the query alignment into a matrix

        current['alignment'] = UtilityFunctions.AlignmentArray(
            [current['consensus']])

        current['matrix'], current['nt_inds'] = Consensus.makeAlignmentMatrix(
                                                current['alignment'])

        current['seqdict'] = seqdict
        # remove the query sequence from X
        X = X[1:]
        consensusD = dict()
        # Build a cluster based on the current query sequence
        X, C = Alignment.buildCluster(X, current, consensusD, pD, k)

        D.setdefault(k, dict())
        lp("Saved cluster %i with %s sequences" % (k, len(C['current_names'])),
           1, pD)
        D[k]['consensus'] = C['current_consensus']
        D[k]['alignment'] = C['current_alignment']
        D[k]['names'] = C['current_names']
        D[k]['log'] = C['seqdict']
        writeFastas(D[k], k, 1, pD['outdir'])
        k += 1
    if len(X) == 1:
        D.setdefault(k, dict())
        lp("Saved cluster %i with %s sequences" % (k, 1),
           1, pD)
        D[k]['consensus'] = X[0][1]
        D[k]['alignment'] = UtilityFunctions.AlignmentArray([X[0][1]])
        D[k]['names'] = [X[0][0]]
        D[k]['log'] = dict()
        writeFastas(D[k], k, 1, pD['outdir'])
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


def writeFastas(result, i, rround, outdir):
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
        outdir, rround, i), "w")
    cons_out = open("%s/round_%i/cluster_%s_consensus.fasta" % (
        outdir, rround, i), "w")
    both_out = open("%s/round_%i/cluster_%s_ali_plus_cons.fasta" % (
        outdir, rround, i), "w")

    ali = result['alignment']
    cons = result['consensus']
    names = result['names']

    for i, nam in enumerate(names):
        ali_out.write(">%s\n%s\n" % (nam, "".join(list(ali[i]))))
        both_out.write(">%s\n%s\n" % (nam, "".join(list(ali[i]))))
    cons_out.write(">consensus_%i\n%s\n" % (i, cons))
    both_out.write(">consensus_%i\n%s\n" % (i, cons))
    ali_out.close()
    cons_out.close()
    both_out.close()
