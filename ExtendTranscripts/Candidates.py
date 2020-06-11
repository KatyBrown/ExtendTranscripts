#!/usr/bin/env python3
import Alignment
import AlignmentSW
import UtilityFunctions
import Consensus
import os
import pyfaidx
os.chdir("/home/katy/ETTest")

def runCandidates(infile, candfile, pD, outdir):
    Z = []
    fasta_dict = pyfaidx.Fasta(infile)
    nams = list(fasta_dict.keys())
    seqdict = {nam: dict() for nam in nams}
    for nam in nams:
        seqdict[nam]['is_rc'] = False
    Z = list(enumerate(nams))
    candidates = UtilityFunctions.FastaToDict(candfile)
    D = dict()
    k = 1
    for c_nam, c_seq in candidates.items():
        current = dict()
        current['name'] = c_nam
        current['consensus'] = c_seq
        current['alignment'] = UtilityFunctions.AlignmentArray([c_seq])
        current['matrix'], current['nt_inds'] = Consensus.makeAlignmentMatrix(
                                                    current['alignment'])
        current['seqdict'] = seqdict
        current['names'] = [c_nam]
        k += 1
        X, C = AlignmentSW.buildClusterLM(Z, fasta_dict, current, pD, k)
        D.setdefault(k, dict())

        D[k]['consensus'] = C['current_consensus']
        D[k]['alignment'] = C['current_alignment']
        D[k]['names'] = C['current_names']
        D[k]['log'] = C['seqdict']
        Alignment.writeFastas(D, k, 1, outdir)
    return(D)

pD = UtilityFunctions.makepDTesting("/home/katy/ETTest/config.ini")

