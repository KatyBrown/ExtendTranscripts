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

    X = copy.copy(Z)
    candidates = UtilityFunctions.FastaToDict(candfile)
    D = dict()
    k = 0
    for c_nam, c_seq in candidates.items():
        current = dict()
        current['name'] = c_nam
        current['consensus'] = c_seq
        if rround == 1:
            current['alignment'] = UtilityFunctions.AlignmentArray([c_seq])
            current['seqdict'] = seqdict
            current['names'] = [c_nam]
        elif rround == 2:
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


def splitAlignment(currentD, outdir, minlen=50):
    ali = currentD[1]['alignment']
    nams = np.array(currentD[1]['names'][1:])
    cons = currentD[1]['consensus']
    newD = dict()
    ali = ali[1:, :]
    not_gap_only = np.where(np.sum(ali == "-", 0) != np.shape(ali)[0])[0]
    cs = np.where(np.diff(not_gap_only) != 1)[0] + 1
    sub_arrs = [((s[0], s[-1])) for s in np.split(not_gap_only, cs)]
    comb_out = open("%s/round_1/consensus_combined.fasta" % outdir, "w")
    for i, sub_arr in enumerate(sub_arrs):
        if (sub_arr[1] - sub_arr[0]) > minlen:
            stem = "%s/round_1/cluster_1_split_%i" % (outdir, i+i)
            ali_out = open("%s_alignment.fasta" % (stem), "w")
            both_out = open("%s_ali_plus_cons.fasta" % (stem), "w")
            cons_out = open("%s_consensus.fasta" % (stem), "w")
            sub_ali = ali[:, sub_arr[0]:sub_arr[1]]
            sub_ali_ng = np.where(
                np.sum(sub_ali == "-", 1) != np.shape(sub_ali)[1])[0]
            sub_ali = sub_ali[sub_ali_ng, :]
            sub_names = nams[sub_ali_ng]
            newD.setdefault(i+1, dict())
            newD[i+1]['alignment'] = sub_ali
            newD[i+1]['names'] = list(sub_names)
            newD[i+1]['log'] = None
            newD[i+1]['seqdict'] = currentD[1]['seqdict']
            matrix, nt_inds = Consensus.makeAlignmentMatrix(sub_ali)
            cons = Consensus.collapseAlignment(matrix, nt_inds)
            newD[i+1]['consensus'] = cons
            for j, nam in enumerate(list(sub_names)):
                seq = "".join(list(sub_ali[j]))
                ali_out.write(">%s\n%s\n" % (nam, seq))
                both_out.write(">%s\n%s\n" % (nam, seq))
            comb_out.write(">*consensus_%s\n%s\n" % (i+1, cons))
            both_out.write(">*consensus_%s\n%s\n" % (i+1, cons))
            ali_out.close()
            both_out.close()
            cons_out.close()
    comb_out.close()
    return (newD)
