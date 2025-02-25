#!/usr/bin/env python3
import Alignment
import UtilityFunctions
import copy
import Consensus
from UtilityFunctions import logPrint as lp
import os
import shutil
import AlignmentSW
import pyfaidx
import Candidates


def runAlignment(fasta_dict, pD, outdir, alignment_type='pairwise',
                 min_length=30,
                 quick=False, candidates=False, reference_dict=False,
                 log=None):

    nams_orig = sorted(list(fasta_dict.keys()))
    nams = []
    for nam in nams_orig:
        if len(fasta_dict[nam][:].seq) > min_length:
            nams.append(nam)
    seqdict = {nam: dict() for nam in nams}
    lp("Reading input file %s" % pD['alignment_infile'], 1, pD)
    lp("%i sequences in input file" % (len(nams)), 2, pD)
    for nam in nams:
        seqdict[nam]['is_rc'] = False

    Z = list(enumerate(nams))

    if alignment_type == 'pairwise':
        D = runAlignmentPW(Z, fasta_dict, pD,
                           quick=quick)
    elif alignment_type == 'stepwise':
        lp("Running clustering and alignment round 1", 1, pD)
        if not candidates:
            D = AlignmentSW.runClusters(Z, fasta_dict, pD, seqdict, 1,
                                        candidates, reference_dict,
                                        log=log)
            Alignment.mergeFastas(1, len(D), outdir)
        else:
            D = Candidates.runCandidates(Z, fasta_dict,
                                         seqdict, reference_dict,
                                         pD, outdir, rround=1)
        lp("Identified %i clusters in round 1" % len(D), 2, pD)
        
        prev_len = len(D)
        i = 2
        current_len = 0
        if candidates:
            D = Candidates.splitAlignment(D, outdir)
            p = "%s/round_1/consensus_combined.fasta" % outdir
            if os.path.exists(p):
                new_cands = "%s/round_1/consensus_combined.fasta" % outdir
                D = Candidates.runCandidates(Z, fasta_dict,
                                             seqdict, new_cands,
                                             pD, outdir, rround=2,
                                             currentD=D)
            i += 1

        currentD = copy.copy(D)
        K = copy.copy(D)
        print (K)
        while current_len < prev_len:
            lp("Running clustering and alignment round %i" % i, 1, pD)
            p = "%s/round_%i/consensus_combined.fasta" % (outdir, i - 1)
            if os.path.exists(p):
                fasta_dict = pyfaidx.Fasta(p)
            else:
                fasta_dict = dict()
            cons_n = []
            cons_s = []
            for key in currentD:
                cons_n.append('*consensus_%i' % key)
                cons_s.append(currentD[key]['consensus'])
                seqdict['*consensus_%i' % key] = dict()
            Z = list(enumerate(cons_n))
            prev_len = len(currentD)

            currentD = AlignmentSW.runClusters(Z, fasta_dict, pD, seqdict, i,
                                               cons=True, currentD=currentD,
                                               log=log)
            current_len = len(currentD)
            if current_len < prev_len and current_len != 0:
                lp("Identified %i clusters in round %i" % (len(currentD), i),
                   2, pD)
                Alignment.mergeFastas(i, len(currentD), outdir)
                K = copy.copy(currentD)
            else:
                lp("Identified no additional clusters in round %i" % (
                    i), 2, pD)
                if os.path.exists("%s/round_%i" % (pD['outdir'], i)):
                    shutil.rmtree("%s/round_%i" % (pD['outdir'], i))

                currentD = K
            i += 1
    lp("Finished aligning and clustering %s" % (pD['alignment_infile']), 1, pD)
    lp("Identified %i clusters from %i fragments" % (len(currentD),
                                                     len(nams)), 1, pD)
    if os.path.exists("%s/final_clusters" % (pD['outdir'])):
        shutil.rmtree("%s/final_clusters" % (pD['outdir']))
    if os.path.exists("%s/round_%i" % (pD['outdir'], i-2)):
        shutil.copytree("%s/round_%i" % (pD['outdir'], i-2),
                        "%s/final_clusters" % pD['outdir'])

    if candidates and len(currentD) != 0:
        Candidates.alignWithRef(currentD, reference_dict, pD, outdir)
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
