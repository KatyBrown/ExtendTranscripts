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
import Bam

def plotbasicCov(subplot, coverage_tab, colour, contig, bam):
    positions = coverage_tab['pos']
    coverage = coverage_tab['coverage']
    if max(coverage_tab['coverage'] != 0):

        subplot.set_xlim(min(positions),
                         max(positions))
        subplot.plot(positions, coverage, color=colour, lw=10)

        subplot.set_ylim(0,
                         max(coverage) * 1.1)
        subplot.text(
            min(positions) + (0.05 * (max(positions) - min(positions))),
            max(coverage) * 1.05,
            "%s, %s" % (contig, bam))
        subplot.set_xlabel("Position in Genome")
        subplot.set_ylabel("Number of Reads")
    else:
        subplot.text(0, 0.5, "No reads in this interval for %s, %s" % (contig,
                                                                       bam))




def plotCov(coverage_tab, bams,
            outdir, stem, contig,
            fasta_dict,
            coverage_tool,
            figdpi, start_interval, end_interval,
            intervals=[], minimum_coverage=1):
    contig_length = max(coverage_tab['pos'])
    intervals.insert(0, ((contig_length - end_interval, contig_length)))
    intervals.insert(0, ((0, start_interval)))
    intervals.insert(0, ((0, contig_length)))
    # Assign one colour to each bam file
    colours = [matplotlib.colors.rgb2hex(x)
               for x in plt.cm.Dark2(np.linspace(0, 1, 8))]
    # add cycles of colours in the case of >20 bam files
    nrounds = math.ceil(len(bams) / 10) - 1
    orig_colours = copy.copy(colours)
    for i in range(nrounds):
        colours += orig_colours

    for bam in bams:
        B, paired, mm = bamToDict(bam, contig, fasta_dict)
        print (len(B))
    for interval in intervals:
        subtab = coverage_tab[(coverage_tab['pos'] >= interval[0]) &
                              (coverage_tab['pos'] <= interval[1])]
        goodbams = []
        for bam in bams:
            subtab2 = subtab[subtab['bam'] == bam]
            if max(subtab2['coverage']) >= minimum_coverage:
                goodbams.append(bam)
        nbams = len(goodbams)
        f = plt.figure(figsize=(25, 5*nbams),
                       dpi=figdpi)
        f = matplotlib.gridspec.GridSpec(nbams, 1)

        # no whitespace between subplots
        f.update(wspace=0, hspace=0)

        for i, b in enumerate(goodbams):
            subplot = plt.subplot(f[i, 0])
            plotbasicCov(subplot, subtab[subtab['bam'] == b],
                     colour=colours[i], contig=stem, bam=b)
            subplot.xaxis.set_visible(False)
        # Put the x axis back for the bottom plot
        subplot.xaxis.set_visible(True)
        
        plt.savefig("%s/%s_coverage_%s_%s_%s.png" % (
            outdir, stem, coverage_tool, interval[0], interval[1]),
            dpi=figdpi, bbox_inches='tight')