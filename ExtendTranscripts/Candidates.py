#!/usr/bin/env python3
import Alignment
import AlignmentSW
import UtilityFunctions
import Consensus
import os
import pyfaidx
import Alignment
import copy
import shutil
import numpy as np
from UtilityFunctions import logPrint as lp


def runCandidates(Z, fasta_dict, seqdict, candfile, pD, outdir, rround,
                  currentD=None):
    '''
    Allows the user to specify a set of reference sequences - in the first
    round of clustering, only contigs which align to these references
    (meeting the same minimum criteria as for normal clustering) are
    selected.  In the second round, the fragments identified in the query
    file which matched the reference sequence are used to identify
    further fragments.  From this point clustering of consensus
    sequences continues as normal.
    
    Parameters
    ----------
    Z: list
        List of two item tuples where the first element is an integer and
        the second the sequence name for all sequences in the main input file
    fasta_dict: pyfaidx.Fasta
        pyfaidx indexed Fasta object containing the main input file of contigs
    seqdict: dict
        A dictionary where keys are sequence IDs and values are empty
        dictionaries - these are used to store CIAlign logs for sequences
        later but it is not run at this stage when candidate sequences are
        used
    candfile: str
        Path to a file containing the reference sequences to match the contigs
        to in the first round of clustering
    pD: dict
        Dictionary containing the initial parameters set by the user
    outdir: str
        Path to directory in which to save all output files
    rround: int
        Round number - which round of clustering is this - used to
        determine where to save the output
    currentD: dict
        Dictionary containing the results of previous rounds of clustering
        used in this case to expand consensus sequences from round 1
    
    Returns
    -------
    D: dict
        Updated version of currentD containing the results of this
        round of clustering
    '''

    X = copy.copy(Z)
    candidates = UtilityFunctions.FastaToDict(candfile)
    D = dict()
    k = 0
    # iterate through the reference sequences
    for c_nam, c_seq in candidates.items():
        current = dict()
        # store candidate name and sequence in the results dictionary
        current['name'] = c_nam
        current['consensus'] = c_seq
        if rround == 1:
            # in the first round, none of the input sequences are
            # consensus sequences
            current['alignment'] = UtilityFunctions.AlignmentArray([c_seq])
            current['seqdict'] = seqdict
            current['names'] = [c_nam]
        elif rround == 2:
            # in the second round, expand the consensus sequences based
            # on the output of the previous round
            consn = int(c_nam.replace("*consensus_", ""))
            current['alignment'] = currentD[consn]['alignment']
            current['seqdict'] = currentD[consn]['seqdict']
            current['names'] = currentD[consn]['names']

        current['matrix'], current['nt_inds'] = Consensus.makeAlignmentMatrix(
                                                    current['alignment'])

        k += 1
        X, C = AlignmentSW.buildCluster(Z, fasta_dict, current, pD, k,
                                        candidate=True, skip=True)

        D.setdefault(k, dict())
        D[k]['consensus'] = C['current_consensus']
        D[k]['alignment'] = C['current_alignment']
        D[k]['names'] = C['current_names']
        D[k]['seqdict'] = C['seqdict']
        Alignment.writeFastas(D[k], k, rround, outdir, candidates=True,
                              reference=candidates)
        if rround == 2:
            Alignment.mergeFastas(2, len(D), outdir)
    return(D)


def splitAlignment(currentD, outdir, minlen=40):
    m = 1
    newD = dict()
    comb_out = open("%s/round_1/consensus_combined.fasta" % outdir, "w")
    
    for k in currentD:

        ali = currentD[k]['alignment']
        nams = np.array(currentD[k]['names'][1:])
        cons = currentD[k]['consensus']
        ali = ali[1:, :]

        not_gap_only = np.where(np.sum(ali == "-", 0) != np.shape(ali)[0])[0]
        cs = np.where(np.diff(not_gap_only) != 1)[0] + 1
        if len(cs) != 0:
            sub_arrs = [((s[0], s[-1])) for s in np.split(not_gap_only, cs)]
            i = 0
            for sub_arr in sub_arrs:
                if (sub_arr[1] - sub_arr[0]) > minlen:
                    stem = "%s/round_1/cluster_%i_split_%i" % (outdir, k, i+1)
                    ali_out = open("%s_alignment.fasta" % (stem), "w")
                    both_out = open("%s_ali_plus_cons.fasta" % (stem), "w")
                    cons_out = open("%s_consensus.fasta" % (stem), "w")
                    sub_ali = ali[:, sub_arr[0]:sub_arr[1]]
                    sub_ali_ng = np.where(
                        np.sum(sub_ali == "-", 1) != np.shape(sub_ali)[1])[0]
                    sub_ali = sub_ali[sub_ali_ng, :]
                    sub_names = nams[sub_ali_ng]
                    newD.setdefault(m, dict())
                    newD[m]['alignment'] = sub_ali
                    newD[m]['names'] = list(sub_names)
                    newD[m]['log'] = None
                    newD[m]['seqdict'] = currentD[1]['seqdict']
                    matrix, nt_inds = Consensus.makeAlignmentMatrix(sub_ali)
                    cons = Consensus.collapseAlignment(matrix, nt_inds)
                    newD[m]['consensus'] = cons
                    for j, nam in enumerate(list(sub_names)):
                        seq = "".join(list(sub_ali[j]))
                        ali_out.write(">%s\n%s\n" % (nam, seq))
                        both_out.write(">%s\n%s\n" % (nam, seq))
                    cons_out.write(">*consensus_%s\n%s\n" % (m, cons))
                    both_out.write(">*consensus_%s\n%s\n" % (m, cons))
                    comb_out.write(">*consensus_%s\n%s\n" % (m, cons))
                    ali_out.close()
                    both_out.close()
                    cons_out.close()
                    i += 1
                    m += 1
    comb_out.close()
    return (newD)


def alignWithRef(currentD, reference_dict, pD, outdir):
    for i in currentD:
        query_seq = currentD[i]['consensus']
        query_ali = currentD[i]['alignment']
        query_names = currentD[i]['names']
        q_matrix, nt_inds = Consensus.makeAlignmentMatrix(query_ali)
        candidates = UtilityFunctions.FastaToDict(reference_dict)
        for target_nam, target_seq in candidates.items():
            qm = copy.copy(q_matrix)
            target_ali = UtilityFunctions.AlignmentArray([target_seq])
            # Align the query and target consensus sequences
            result = Alignment.SWalign(query_seq, target_seq,
                                       pD, useSub=True)

            # Check the if the consensus sequences are a good match
            is_match = Alignment.alignmentMeetsCriteria(result, query_seq,
                                                        target_seq, pD)

            if is_match[0]:
                result['alignment'] = is_match[1]
                # get the full alignment for the two consensus sequences
                result = Alignment.getAlignmentFull(result,
                                                    query_seq,
                                                    target_seq,
                                                    pD)

                ali, matrix = Consensus.expandAlignment(result,
                                                        query_ali,
                                                        target_ali,
                                                        qm,
                                                        nt_inds)
                cons = Consensus.collapseAlignment(matrix, nt_inds)
                names = query_names + [target_nam]
                result2 = Alignment.SWalign(query_seq, cons, pD, useSub=True)
                is_match_2 = Alignment.alignmentMeetsCriteria(result2,
                                                              query_seq,
                                                              cons,
                                                              pD)
                result2['alignment'] = is_match_2[1]
                result2 = Alignment.getAlignmentFull(result2,
                                                     query_seq,
                                                     cons,
                                                     pD)
                a2 = UtilityFunctions.AlignmentArray([query_seq])
                ali, matrix = Consensus.expandAlignment(result2,
                                                        a2,
                                                        ali,
                                                        qm,
                                                        nt_inds)
                names = ["*consensus_%s" % i] + names
                path = "%s/final_clusters/consensus_%s_to_%s_ali.fasta" % (
                    outdir, i, target_nam.split(" ")[0])
                out = open(path, "w")
                for j, nam in enumerate(names):
                    out.write(">%s\n%s\n" % (nam,
                                             "".join(list(ali[j]))))
                out.close()
