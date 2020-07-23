#!/usr/bin/env python3

import os
import logging
import configargparse
import sys
import pyfaidx
sys.path.append("/home/katy/ExtendTranscripts/ExtendTranscripts")
import runAlignment
import UtilityFunctions
from UtilityFunctions import logPrint as lp
import Bam

class ExtendAction(configargparse.Action):
    '''
    Adds the "extend" action to argparse - required for Python <3.8
    '''
    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest) or []
        items.extend(values)
        setattr(namespace, self.dest, items)


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
    parser.register('action', 'extend', ExtendAction)
    general = parser.add_argument_group('General Options')
    alignment = parser.add_argument_group("Alignment Options")
    breakpoints = parser.add_argument_group("Breakpoint Options")


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

    general.add("-a", "--alignment", dest="alignment", action="store_true",
                help="""Specify this argument to run the alignment algorithm \
                to combine contigs.  Default %(default)s. Type %(type)s""",
                default=False)

    general.add("-b", "--breakpoints", dest="breakpoints",
                action="store_true",
                help="""Specify this argument to run the algorithm to check \
                for possible breakpoints in a fasta file of contigs based on \
                read coverage. Default %(default)s. Type %(type)s""",
                default=False)

    alignment.add("-i", "--infile", dest='alignment_infile',
                  type=str,
                  help="""Path to input alignment file in FASTA format. \
                  Type: %(type)s.""", default=None)

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

    alignment.add("--remove_insertions_max_n",
                  dest="alignment_remove_insertions_max_n",
                  type=int,
                  help="""Remove insertions which are not in the majority \
                  of sequences up to a maximum of this size from the \
                  alignments as they are being built""",
                  default=50)

    breakpoints.add("--contig_fasta_file", dest="breakpoint_contigs",
                    type=str, help="""Path to fasta file containing contigs \
                    to test for possible breakpoints based on read coverage""",
                    default=None)

    breakpoints.add("--bam_files", dest="breakpoint_bams",
                    action="extend", nargs="+", type=str, help="""List of bam \
                    files to use to test for possible breakpoints in the \
                    contig. All bam files in this list will be combined. \
                    Default %(default)s.  Type: %(type)s.""",
                    default=[])

    breakpoints.add("--coverage_tool", dest="breakpoint_coverage_tool",
                    type=str, help="""Tool to use to calculate coverage of \
                    contig positions - can be bedtools coverage ("bedtools") \
                    or samtools mpileup ("samtools") Default %(default)s.  \
                    Type: %(type)s.""",
                    default="samtools")

    breakpoints.add("--figdpi", dest="breakpoint_figdpi",
                    type=int, help="""DPI for figures showing coverage \
                    of contigs by reads. Default %(default)s.  \
                    Type: %(type)s.""",
                    default=300)

    breakpoints.add("--start_interval", dest="breakpoint_start_interval",
                    type=int, help="""Create a zoomed coverage plot for this \
                    many positions at the beginning of the sequence. \
                    Default %(default)s.  \
                    Type: %(type)s.""",
                    default=100)

    breakpoints.add("--end_interval", dest="breakpoint_end_interval",
                    type=int, help="""Create a zoomed coverage plot for this \
                    many positions at the end of the sequence. \
                    Default %(default)s.  \
                    Type: %(type)s.""",
                    default=100)

    breakpoints.add("--additional_interval",
                    dest="breakpoint_additional_interval",
                    action="extend", nargs="+", type=int, help="""Additional
                    intervals to create zoomed coverage plots for. These
                    are specified as --additional_interval start_pos end_pos \
                    the same argument can be reused as many times as needed. \
                    Default %(default)s.  Type: %(type)s.""",
                    default=[])

    breakpoints.add("--min_coverage", dest="breakpoint_min_coverage",
                    type=int, help="""Minimum number of reads in an interval \
                    before coverage is plotted. Default %(default)s. \
                    Type: %(type)s""", default=1)

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
    pD['print'] = True

    lp("Initial Parameters %s" % str(parser.format_values()), level=1, pD=pD)

    if args.alignment:
        if args.alignment_infile is None:
            raise RuntimeError(
                "To align contigs a fasta file of contigs must be specified")
        fasta_dict = pyfaidx.Fasta(args.alignment_infile)
        runAlignment.runAlignment(fasta_dict, pD, args.outdir,
                                  alignment_type=args.alignment_type,
                                  quick=False)

    if args.breakpoints:
        if args.breakpoint_contigs is None:
            raise RuntimeError(
                "To find breakpoints a fasta file of contigs must be \
specified")
        if len(args.breakpoint_bams) == 0:
            raise RuntimeError(
                "To find breakpoints bam files containing reads must be \
specified")
        fasta_dict = pyfaidx.Fasta(args.breakpoint_contigs)
        Bam.runAll(fasta_dict,
                   sorted(args.breakpoint_bams),
                   args.outdir,
                   coverage_tool=args.breakpoint_coverage_tool,
                   figdpi=args.breakpoint_figdpi,
                   start_interval=args.breakpoint_start_interval,
                   end_interval=args.breakpoint_end_interval,
                   add_intervals=args.breakpoint_additional_interval,
                   min_coverage=args.breakpoint_min_coverage)


if __name__ == "__main__":
    main()
