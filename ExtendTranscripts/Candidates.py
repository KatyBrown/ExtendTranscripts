#!/usr/bin/env python3
import Alignment
import UtilityFunctions


def runCandidates(infile, candfile, pD):
    candidates = UtilityFunctions.FastaToDict(candfile)
    seq = []
    with open(infile) as inf:
        for line in inf:
            if line.startswith(">"):
                seq = "".join(seq)
                for c_nam, c_seq in candidates.items():
                    result = Alignment.SWalign(seq, c_seq,
                                               pD, useSub=True)
                    is_match = Alignment.alignmentMeetsCriteria(result, seq,
                                                                c_seq, pD)
                    print (is_match)
                seq = []
                nam = line[1:]
            else:
                seq.append(line)

