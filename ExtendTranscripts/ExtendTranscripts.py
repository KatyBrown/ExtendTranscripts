#!/usr/bin/env python3

import os
import logging
import configargparse
import runAlignment
import UtilityFunctions
from UtilityFunctions import logPrint as lp


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
                description='''Extend contigs''', add_help=True)
    required = parser.add_argument_group('Required Arguments')
    general = parser.add_argument_group('General Options')
    alignment = parser.add_argument_group("Alignment Options")

    required.add("-i", "--infile", dest='infile',
                 type=str,
                 required=True,
                 help="""Path to input alignment file in FASTA format. \
                 Type: %(type)s.""")

    general.add("-o", "--outdir", dest='outdir',
                type=str,
                help="""Path to output directory. \
                Default: %(default)s. Type: %(type)s.""",
                default="ET")

    general.add("-c", "--conf_file", dest='conf_file', type=str,
                default=None,
                help="""Path to config file. \
                Default: %(default)s. Type: %(type)s.""",
                is_config_file=True)

    general.add("-vl", "--logging_verbosity", dest="logging_verbosity",
                type=int,
                help="""Logging verbosity level: 0: silent, 1: minimal,
                2: default 3: detailed. Default %(default)s. Type %(type)s""",
                default=2)

    general.add("-vs", "--stdout_verbosity", dest="stdout_verbosity",
                type=int,
                help="""Verbosity level to STDOUT: 0: silent, 1: minimal,
                2: default 3: detailed. Default %(default)s. Type %(type)s""",
                default=2)

    alignment.add("--type", dest="alignment_type",
                  type=str,
                  choices=['stepwise', 'pairwise'],
                  help="""Alignment type to run - stepwise or pairwise. \
                  Default: %(default)s. Type: %(type)s.""")

    alignment.add("--gap_open", dest="alignment_gap_open",
                  type=int,
                  help="""Gap opening penalty. \
                  Default: %(default)s. Type: %(type)s.""",
                  default=5)

    alignment.add("--gap_extend", dest="alignment_gap_extend",
                  type=int,
                  help="""Gap extension penalty. \
                  Default: %(default)s. Type: %(type)s.""",
                  default=2)

    alignment.add("--match_score", dest="alignment_match_score",
                  type=int,
                  help="""Score for each match in SW algorithm. \
                  Default: %(default)s. Type: %(type)s.""", default=5)

    alignment.add("--mismatch_score", dest="alignment_mismatch_score",
                  type=int,
                  help="""Penalty for each mismatch in SW algorithm.
                  Default: %(default)s. Type: %(type)s.""",
                  default=-4)

    alignment.add("--clip_max_perc", dest="alignment_clip_max_perc",
                  type=floatOpt,
                  help="""Maximum proportion of nucleotides in relation to \
                  total sequence length to clip when aligning with SW, \
                  overridden if clip_max_n is specified.
                  Default: %(default)s. Type: float""",
                  default=0.05)

    alignment.add("--clip_max_n", dest="alignment_clip_max_n",
                  type=intOpt,
                  help="""Maximum number of nucleotides to clip when aligning \
                  with SW, overrides clip_max_perc if specified.
                  Default: %(default)s. Type: int""",
                  default=None)

    alignment.add("--ident_min_perc", dest="alignment_ident_min_perc",
                  type=floatOpt,
                  help="""Minimum percentage identical nucleotides in \
                  overlapping \
                  region when aligning with SW, overridden if ident_min_n \
                  is specified. Default: %(default)s. Type: float""",
                  default=0.9)

    alignment.add("--ident_min_n", dest="alignment_ident_min_n",
                  type=intOpt,
                  help="""Maximum number of nucleotides to clip when aligning \
                  with SW, overrides ident_min_perc if specified.
                  Default: %(default)s. Type: int""",
                  default=None)

    alignment.add("--length_min_perc", dest="alignment_length_min_perc",
                  type=floatOpt,
                  help="""Minimum percentage of sequence in overlapping \
                  region when aligning with SW, overridden length_min_n is \
                  specified. Default: %(default)s. Type: float""",
                  default=0.9)

    alignment.add("--length_min_n", dest="alignment_length_min_n",
                  type=intOpt,
                  help="""Minimum number of nucleotides in overlapping region \
                  when aligning with SW, overrides length_min_perc if \
                  specified.  Default: %(default)s.  Type: int""",
                  default=None)

    alignment.add("--indels_max_perc", dest="alignment_indels_max_perc",
                  type=floatOpt,
                  help="""Minimum percentage of indel positions (gaps in one \
                  sequence) in overlapping region when aligning with SW, \
                  overridden if indels_max_n is specified.
                  Default: %(default)s. Type: float""",
                  default=0.15)

    alignment.add("--indels_max_n", dest="alignment_indels_max_n",
                  type=intOpt,
                  help="""Minimum number of indel positions (gaps in one\
                  sequence) in overlapping region when aligning with SW, \
                  overrides indels_max_perc if specified.
                  Default: %(default)s. Type: int""",
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
    formatter = logging.Formatter("""%(asctime)s: %(message)s""")
    handler.setFormatter(formatter)

    # add the handlers to the logger
    log.addHandler(handler)

    pD = vars(args)
    pD['log'] = log

    lp("Initial Parameters %s" % str(parser.format_values()), level=1, pD=pD)

    fasta_dict = UtilityFunctions.FastaToDict(args.infile)
    runAlignment.runAlignment(fasta_dict, pD, args.outdir,
                              alignment_type=args.alignment_type,
                              quick=False)


if __name__ == "__main__":
    main()
