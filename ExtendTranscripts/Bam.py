#!/usr/bin/env python3
import pybedtools
import pysam
import pandas as pd
import UtilityFunctions
import os
import warnings
import pyfaidx
import math
import numpy as np
import copy
import BamPlots


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


def calcAltCoverage(bam_dict,
                    contig_length,
                    buffer_prop,
                    rl):
    '''
    Calculate coverage at each position by terminal and
    non-terminal positions in reads, insert region with paired
    end reads.
    '''
    # calculate buffer size based on read length
    buffer = math.ceil(rl * buffer_prop)

    # keep track of the coverage of each type at each
    # postion by adding to these vectors at the
    # appropriate index

    # coverage in left terminal of reads (0:buffer)
    starts = np.zeros((1, contig_length))
    # coverage in right terminal of reads (-buffer:end)
    ends = np.zeros((1, contig_length))

    # coverage in bodies of reads (buffer:-buffer)
    bodies = np.zeros((1, contig_length))

    # coverage in inserts - between read1 and read2 for paired
    # end reads
    inserts = np.zeros((1, contig_length))

    # record the co-ordinates of the read
    for read_ID in bam_dict:
        if bam_dict[read_ID]['has_pair']:
            read1 = bam_dict[read_ID]['read1']
            read2 = bam_dict[read_ID]['read2']
            # make a range covering all the positions in the read
            read1_positions = np.arange(read1['start'], read1['end'])
            read2_positions = np.arange(read2['start'], read2['end'])

            # add one to all positions within the buffers
            starts[0, read1_positions[:buffer]] += 1
            ends[0, read2_positions[-buffer:]] += 1

            # add one to all positions within the read bodies
            bodies[0, read1_positions[buffer:]] += 1
            bodies[0, read2_positions[:-buffer]] += 1

            # add one to all positions between read1 and read2
            inte = np.arange(read1_positions[-1], read2_positions[0])
            inserts[0, inte] += 1

        else:
            # make a range covering all the positions in the read
            read = [bam_dict[read_ID]['unpaired']]
            read_positions = np.arange(read['start'], read['end'])

            # add one to all the positions within the buffers
            starts[0, read_positions[:buffer]] += 1
            ends[0, read_positions[-buffer:]] += 1

            # add one to all the positions within read bodies
            bodies[0, read_positions[buffer:]] += 1

    return (starts, ends, bodies, inserts)


def getReadArray(bam_dict, paired):
    if paired:
        arr = np.empty([4, len(bam_dict)])
    else:
        arr = np.empty([2, len(bam_dict)])
    nams = np.empty(len(bam_dict), dtype='object')
    strands = np.empty(len(bam_dict), dtype='object')
    for i, (nam, read) in enumerate(bam_dict.items()):
        if paired:
            span = np.array([read['read1']['start'],
                             read['read1']['end'],
                             read['read2']['start'],
                             read['read2']['end']])
            strands[i] = read['read1']['strand']
        else:
            span = np.array([read['unpaired']['start'],
                             read['unpaired']['end']])
            strands[i] = read['unpaired']['strand']
        arr[:, i] = span
        nams[i] = nam
        
    return (nams, arr, strands)


def bamToDict(bam, contig, contig_seq):
    # read the bam file
    samfile = pysam.AlignmentFile(bam)

    # is the library paired end?
    # changes to True later if it is
    paired = False
    D = dict()
    # number of reads not identical to reference
    mm = 0
    # read lengths
    rlens = []
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
                if not line.is_reverse:
                    thisdict = 'read1'
                    strand = "+"
                else:
                    thisdict = 'read2'
                    strand = "-"
                paired = True
                D[read_ID]['has_pair'] = True
            elif line.is_read2 and line.is_proper_pair:
                if line.is_reverse:
                    thisdict = 'read2'
                    strand = "+"
                else:
                    thisdict = 'read1'
                    strand = "-"
                D[read_ID]['has_pair'] = True
            else:
                thisdict = 'unpaired'
                D[read_ID]['has_pair'] = False
            # each read1, read2 or unpaired dict has
            # start - start position
            # end - end position
            # strand - + or -
            # if the sequence is different to the reference
            # they also have "seq"
            D[read_ID][thisdict]['start'] = line.pos
            D[read_ID][thisdict]['end'] = line.aend
            D[read_ID][thisdict]['strand'] = strand

            read_ref = contig_seq[line.pos:line.aend]

            if line.seq != read_ref:
                this_read = np.array(list(line.seq))
                ref_read = np.array(list(read_ref))
                # positions of differences
                dpos = np.where(this_read != ref_read)[0]
                # nts of differences
                dnt = this_read[dpos]
                D[read_ID][thisdict]['seq'] = (dpos, dnt)
                mm += 1
            rlens.append(line.rlen)
        return(D, samfile, paired, mm, np.median(rlens))
    else:
        raise RuntimeWarning("""The bam file %s is not mapped \
                             to contig %s""" % bam, contig)
        return (None, None, None, None, None)


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
           min_coverage=1,
           buffer_prop=0.1,
           coverage_lim=1000):

    # run each contig / segment indivdually
    for nam in fasta_dict.keys():
        # get the contig name
        contig = nam.split(" ")[0]

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
        out.write("%s\t0\t%s\n" % (contig, contig_length))
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
        for bam_file in bam_files:
            bamnam = os.path.splitext(os.path.basename(bam_file))[0]
            bamD, samfile, paired, mm, rl = bamToDict(bam_file, contig,
                                                      contig_seq)
            if bamD is not None and len(bamD) > min_coverage:
                coverage_tab = calculateCoverage(contig_file,
                                                 contig,
                                                 bam_file,
                                                 bamnam,
                                                 bedtool,
                                                 bednam,
                                                 outdir,
                                                 coverage_tool)
                altCov = calcAltCoverage(bamD,
                                         contig_length,
                                         buffer_prop,
                                         rl)
                readArr = getReadArray(bamD, paired)
                BamPlots.plotAll(coverage_tab,
                                 altCov,
                                 bam_file,
                                 bamnam,
                                 bamD,
                                 readArr,
                                 paired,
                                 mm,
                                 rl,
                                 outdir,
                                 contig,
                                 contig_dict,
                                 figdpi,
                                 intervals,
                                 coverage_lim)


