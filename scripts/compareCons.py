#!/usr/bin/env python3
import numpy as np
import configargparse
import sys
sys.path.append("/home/katy/ExtendTranscripts/ExtendTranscripts")
import UtilityFunctions
import Alignment


def main():
    parser = configargparse.ArgumentParser(
        description='''Compare two sequences and calculate the % identity''')
    parser.add("--consensus", dest="consensus",
               type=str,
               help="path to consensus sequence fasta file")
    parser.add("--original", dest="original",
               type=str,
               help="path to original input file")
    parser.add("--orig_name", dest="orig_name",
               type=str,
               help="original sequence name")
    parser.add("--gap_open", dest="alignment_gap_open",
               type=int,
               help="gap opening penalty", default=5)
    parser.add("--gap_extend", dest="alignment_gap_extend",
               type=int,
               help="gap extension penalty", default=2)
    parser.add("--match_score", dest="alignment_match_score",
               type=int,
               help="score for each match in SW algorithm", default=5)
    parser.add("--mismatch_score", dest="alignment_mismatch_score",
               type=int,
               help="penalty for each mismatch in SW algorithm", default=-4)

    args = parser.parse_args()

    consD = UtilityFunctions.FastaToDict(args.consensus)
    origD = UtilityFunctions.FastaToDict(args.original)
    cons_seq_f = list(consD.values())[0]
    cons_seq_f = cons_seq_f.upper()

    cons_seq_r = UtilityFunctions.reverseComplement(cons_seq_f)
    orig_seq = origD[args.orig_name]
    orig_seq = orig_seq.upper()
    pD = {'alignment_gap_open': args.alignment_gap_open,
          'alignment_gap_extend': args.alignment_gap_extend,
          'alignment_match_score': args.alignment_match_score,
          'alignment_mismatch_score': args.alignment_mismatch_score}
    ali_f = Alignment.SWalign(cons_seq_f, orig_seq, pD)
    ali_r = Alignment.SWalign(cons_seq_r, orig_seq, pD)

    if ali_f['optimal_alignment_score'] > ali_r['optimal_alignment_score']:
        ali = ali_f
        cons_seq = cons_seq_f
    else:
        ali = ali_r
        cons_seq = cons_seq_r
    q, t, ident, c = Alignment.getAlignmentLocal(ali, cons_seq, orig_seq, pD,
                                                 getcount=True)
    cons_length = len(cons_seq)
    orig_length = len(orig_seq)
    cons_perc_aligned = (ali['query_end'] - ali['query_start']) / cons_length
    orig_perc_aligned = (ali['target_end'] - ali['target_start']) / orig_length

    print (c, ident, cons_perc_aligned, orig_perc_aligned)


if __name__ == "__main__":
    main()