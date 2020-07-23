#!/usr/bin/env python3
import pybedtools
import pysam
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import UtilityFunctions
import os
import warnings
import pyfaidx
import math
import numpy as np
import random
import copy


def calculateCoverage(contig_file,
                      contig,
                      bam_file,
                      bamnam,
                      bedtool,
                      bednam,
                      outdir,
                      typ='samtools'):

    out_tab = "%s/%s_coverage_%s_%s.tsv" % (outdir, contig,
                                            bamnam, typ)
    # index the bam file (if needed)
    if not os.path.exists("%s.bai" % bam_file):
        pysam.index(bam_file)

    if typ == 'bedtools':
        with warnings.catch_warnings():
            # this raises a warning in Python 3.8 that is not relevant
            # ignore it
            # https://github.com/benoitc/gunicorn/issues/2091
            warnings.filterwarnings("ignore")
            bam = pybedtools.BedTool(bam_file)
            # d=True - calculate depth of coverage for every position
            # in the bam file
            cov = bedtool.coverage(bam, d=True)
        cov_df = cov.to_dataframe(names=['chromosome', 'start',
                                         'end', 'pos', 'coverage'])

        cov_df['bam'] = bamnam
    elif typ == 'samtools':
        # -A - keep anomalous read pairs
        # -B - don't calculate per base alignment quality
        # -C 0 - don't adjust based on mapping quality
        # -d 0 - don't limit coverage depth
        # -Q 0 - minimum quality score of 0
        # -aa - output every position even if coverage is 0
        cov = pysam.mpileup(bam_file, "-f", contig_file,
                            "-A", "-B", "-C", "0",
                            "-d", "0", "-Q", "0", "-aa")

        # convert the coverage table into a dataframe
        cov = [x.split("\t") for x in cov.split("\n")]
        
        cov_df = pd.DataFrame(cov, columns=['chromosome',
                                            'pos',
                                            'reference_base',
                                            'coverage',
                                            'base_quality',
                                            'alignment_quality'])
        cov_df['bam'] = bamnam
        
        # the final position is always NA for some reason
        cov_df['coverage'] = cov_df['coverage'].fillna(0)
        cov_df = cov_df[cov_df['pos'].notnull()]
        
        # convert to int
        cov_df['coverage'] = cov_df['coverage'].astype(int)
        cov_df['pos'] = cov_df['pos'].astype(int)

    # save the coverage table
    cov_df.to_csv(out_tab, sep="\t", index=None)
    return (cov_df)


def bamToDict(bam, contig, contig_seq):
    # make a forward and reverse string
    # of the contig sequence
    ref_fwd = contig_seq
    ref_rev = UtilityFunctions.reverseComplement(contig_seq)

    # read the bam file
    samfile = pysam.AlignmentFile(bam)

    # is the library paired end?
    # changes to True later if it is
    paired = False
    D = dict()
    # number of reads not identical to reference
    mm = 0
    # check this bam is actually mapped to this contig
    if contig in samfile.references:
        for line in samfile.fetch(contig):
            read_ID = line.qname.split("/")[0]

            D.setdefault(read_ID, dict())

            # each read ID is a dictionary in D with the keys
            # read1, read2 and unpaired
            D[read_ID].setdefault('read1', dict())
            D[read_ID].setdefault('read2', dict())
            D[read_ID].setdefault('unpaired', dict())

            # call the leftmost read read1 and the rightmost
            # read2 for convenience
            if line.is_read1 and line.is_proper_pair:
                if line.is_reverse:
                    thisdict = 'read2'
                else:
                    thisdict = 'read1'
                paired = True
            elif line.is_read2 and line.is_proper_pair:
                if line.is_reverse:
                    thisdict = 'read1'
                else:
                    thisdict = 'read2'
            else:
                thisdict = 'unpaired'

            # each read1, read2 or unpaired dict has
            # start - start position
            # end - end position
            # strand - + or -
            # if the sequence is different to the reference
            # they also have "seq"
            D[read_ID][thisdict]['start'] = line.pos
            D[read_ID][thisdict]['end'] = line.aend

            # assign strand and get appropriate reference
            if not line.is_reverse:
                D[read_ID][thisdict]['strand'] = "+"
                ref = ref_fwd
            else:
                D[read_ID][thisdict]['strand'] = "-"
                ref = ref_rev
            read_ref = ref[line.pos:line.aend]
            # if the read is not identical to these
            # positions on the reference
            if line.seq != read_ref:
                this_read = np.array(list(line.seq))
                ref_read = np.array(list(read_ref))
                # positions of differences
                dpos = np.where(this_read != ref_read)[0]
                # nts of differences
                dnt = this_read[dpos]
                D[read_ID][thisdict]['seq'] = (dpos, dnt)
                mm += 1
        return(D, paired, mm)
    else:
        raise RuntimeWarning("""The bam file %s is not mapped \
                             to contig %s""" % bam, contig)
        return (None, None, None)


def getIntervals(start_interval, end_interval,
                 add_intervals,
                 contig_length):
    '''
    Compute the intervals over which to visualise coverage etc.
    add_intervals are additional intervals requested by the user,
    initially processed as a list of ints alternating between
    start and end position.
    '''
    # whole genome
    intervals = [((0, contig_length))]
    # first n nts based on start_interval
    intervals.append((0, start_interval))
    # last n nts based on end_interval
    intervals.append((contig_length - end_interval,
                      contig_length))

    # additional user specified intervals
    for i in range(0, len(add_intervals), 2):
        intervals.append((add_intervals[i], add_intervals[i+1]))

    return (intervals)


def runAll(fasta_dict,
           bam_files,
           outdir,
           coverage_tool='samtools',
           figdpi=300,
           start_interval=100,
           end_interval=100,
           add_intervals=[],
           min_coverage=1):

    # run each contig / segment indivdually
    for nam in fasta_dict.keys():
        # get the contig name
        contig = nam.split(" ")[0].split(".")[0]
        
        contig_file = "%s/%s.fasta" % (outdir, contig)
        # make a fasta file for just this contig
        out = open(contig_file, "w")
        out.write(">%s\n%s\n" % (nam, fasta_dict[nam][0:].seq))
        out.close()

        # index it and load into a dict
        contig_dict = pyfaidx.Fasta(contig_file)

        # check how long it is
        contig_length = len(contig_dict[contig])

        # make a bed file of contig lengths for this contig
        bednam = "%s/%s.bed" % (outdir, contig)
        out = open(bednam, "w")
        out.write("%s\t0\t%s\n" % contig, contig_length)
        out.close()

        # store a pybedtools.BedTool object for the
        # full contig
        bedtool = pybedtools.BedTool(bednam)

        # get the sequence
        contig_seq = contig_dict[contig][:].seq

        # get the intervals to plot for this contig
        intervals = getIntervals(start_interval,
                                 end_interval,
                                 add_intervals,
                                 contig_length)

        # read each bam file for this contig
        for bam in bam_files:
            bamnam = os.path.splitext(os.path.basename(bam))[0]
            bamD, paired, mm = bamToDict(bam, contig, contig_dict)
            if bamD is not None:
                coverage_tab = calculateCoverage(contig_file,
                                                 contig,
                                                 bam_file,
                                                 bamnam,
                                                 bedtool,
                                                 bednam,
                                                 outdir,
                                                 coverage_tool)


'''
        plotCov(coverage_tab,
                bam_files,
                outdir,
                stem,
                contig=nam,
                fasta_dict=fasta_dict,
                coverage_tool=coverage_tool,
                figdpi=figdpi,
                start_interval=start_interval,
                end_interval=end_interval,
                intervals=intervals,
                minimum_coverage=min_coverage)
'''

