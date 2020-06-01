#!/usr/bin/env python3

import os
import logging
import configargparse
import runAlignment
import UtilityFunctions

def intOpt(val):
    if val.capitalize() == "None":
        return (None)
    else:
        return (int(val))

def floatOpt(val):
    if val.capitalize() == "None":
        return (None)
    else:
        return(float(val))

def main():
    parser = configargparse.ArgumentParser(
                description='''Extend contigs''',
                add_help=False)
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    
    required.add("-i", "--infile", dest='infile',
                 type=str,
                 required=True,
                 help='Path to input alignment file in FASTA format')
    
    required.add("-a", "--alignment_type", dest="alignment_type",
                 type=str,
                 choices=['stepwise', 'pairwise'],
                 required=True,
                 help="Alignment type to run - stepwise or pairwise")

    optional.add("-o", "--outdir", dest='outdir',
                 type=str, 
                 help="Path to output directory. Default: %(default)s",
                 default="ET")

    optional.add("-c", "--conf_file", dest='conf_file', type=str,
                 default=None,
                 help='Path to config file. Default: %(default)s',
                 is_config_file=True)

    optional.add("--gap_open", dest="alignment_gap_open",
                 type=int,
                 help="gap opening penalty", default=5)
    optional.add("--gap_extend", dest="alignment_gap_extend",
                 type=int,
                 help="gap extension penalty", default=2)
    optional.add("--match_score", dest="alignment_match_score",
                 type=int,
                 help="score for each match in SW algorithm", default=5)
    optional.add("--mismatch_score", dest="alignment_mismatch_score",
                 type=int,
                 help="penalty for each mismatch in SW algorithm", default=-4)
    
    
    optional.add("--clip_max_perc", dest="alignment_clip_max_perc",
                 type=floatOpt,
                 help="maximum proportion of nucleotides in relation to total sequence length to clip when aligning with SW, overridden if clip_max_n is specified",
                 default=0.05)

    optional.add("--clip_max_n", dest="alignment_clip_max_n",
                 type=intOpt,
                 help="maximum number of nucleotides to clip when aligning with SW, overrides clip_max_perc if specified",
                 default=None)

    optional.add("--ident_min_perc", dest="alignment_ident_min_perc",
                 type=floatOpt,
                 help="minimum % identical nucleotides in overlapping region when aligning with SW, overridden if ident_min_n is specified",
                 default=0.9)

    optional.add("--ident_min_n", dest="alignment_ident_min_n",
                 type=intOpt,
                 help="maximum number of nucleotides to clip when aligning with SW, overrides ident_min_perc if specified",
                 default=None)

    optional.add("--length_min_perc", dest="alignment_length_min_perc",
                 type=floatOpt,
                 help="minimum percentage of sequence in overlapping region when aligning with SW, overridden length_min_n is specified",
                 default=0.9)

    optional.add("--length_min_n", dest="alignment_length_min_n",
                 type=intOpt,
                 help="minimum number of nucleotides in overlapping region when aligning with SW, overrides length_min_perc if specified",
                 default=None)
    
    optional.add("--indels_max_perc", dest="alignment_indels_max_perc",
                 type=floatOpt,
                 help="minimum percentage of indel positions (gaps in one sequence) in overlapping region when aligning with SW, overridden if indels_max_n is specified",
                 default=0.15)

    optional.add("--indels_max_n", dest="alignment_indels_max_n",
                 type=intOpt,
                 help="minimum number of indel positions (gaps in one sequence) in overlapping region when aligning with SW, overrides indels_max_perc if specified",
                 default=None)
        
    args = parser.parse_args()
    
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    logfile = "%s/ExtendTranscripts_log.txt" % args.outdir
    handler = logging.FileHandler(logfile)
    handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    # add the handlers to the logger
    log.addHandler(handler)

    log.info("\n\nInitial parameters\n------------------\n\n%s\n\n" % str(parser.format_values()))    

    pD = vars(args)
    
    fasta_dict = UtilityFunctions.FastaToDict(args.infile)
    D = runAlignment.runAlignment(fasta_dict, pD,
                                  alignment_type=args.alignment_type,
                                  quick=False, rename=True)
    print (D)
                                  

if __name__ == "__main__":
    main()