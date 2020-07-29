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


def plotAll(coverage_tab,
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
            coverage_lim):

    for interval in intervals:
        subtab = coverage_tab[(coverage_tab['pos'] >= interval[0]) &
                              (coverage_tab['pos'] <= interval[1])]
        if max(subtab['coverage'] < coverage_lim):
            nfunctions = 3
        else:
            nfunctions = 2
        f = plt.figure(figsize=(25, 5*nfunctions),
                       dpi=figdpi)
        f = matplotlib.gridspec.GridSpec(nfunctions, 1)
        subplot = plt.subplot(f[0, 0])
        plotbasicCov(subplot, subtab,
                     colour='crimson', contig=contig, bam=bamnam)
        subplot = plt.subplot(f[1, 0])
        plotAltCov(subplot, altCov, interval, paired, rl)

        if max(subtab['coverage']) < coverage_lim:
            subplot = plt.subplot(f[2, 0])
            subReadArr = getInRange(readArr[0], readArr[1],
                                    readArr[2],
                                    interval[0], interval[1])
            plotReads(subplot, subReadArr,
                          interval, paired, rl)
        plt.savefig("%s/%s_coverage_%s_%s_%s.png" % (
                    outdir, contig, bamnam, interval[0], interval[1]),
                    dpi=figdpi, bbox_inches='tight')
        plt.close()


def getInRange(nams, arr, strands, start_pos, end_pos):
    which = np.sum(
        (arr[0:4, :] >= start_pos) & (arr[0:4, :] <= end_pos), 0) != 0
    whichnams = nams[which]
    whichstrand = strands[which]
    whicharr = arr[:, which]
    return(whichnams, whicharr, whichstrand)


def plotbasicCov(subplot, coverage_tab, colour, contig, bam):
    positions = coverage_tab['pos']
    coverage = coverage_tab['coverage']
    if max(coverage_tab['coverage'] != 0):

        subplot.set_xlim(min(positions),
                         max(positions))
        subplot.plot(positions, coverage, color=colour, lw=3)

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

def plotAltCov(subplot, altCov, interval, paired, rl):
    xrange = np.arange(interval[0], interval[1])
    starts, ends, bodies, inserts = altCov
    
    starts_sub = (starts + 1)[0:,interval[0]: interval[1]][0]
    ends_sub = (ends + 1)[0:, interval[0]: interval[1]][0]
    bodies_sub = (bodies + 1)[0:, interval[0]: interval[1]][0]
    if paired:
        inserts_sub = (inserts + 1)[0:, interval[0]: interval[1]][0]
    

    subplot.scatter(xrange, starts_sub / bodies_sub, color='#D9514E',
                    marker='.', s=50)
    subplot.scatter(xrange, ends_sub / bodies_sub, color='#2DA8D8',
                    marker='.', s=50)
    
    if paired:
        subplot.scatter(xrange, inserts_sub / bodies_sub, color='#2A2B2D',
                        marker=".", s=50)
    
    subplot.vlines(rl, 0, subplot.get_ylim()[1], ls='dotted')
    subplot.vlines(interval[1] - rl, 0, subplot.get_ylim()[1], ls='dotted')
    subplot.set_xlim(interval[0], interval[1])

def plotReads(subplot, readArr, interval, paired, rl):
    subplot.set_xlim(interval)
    interval_len = interval[1] - interval[0]
    
    nams, arr, strands = readArr
    ys = np.arange(0, np.shape(arr)[1])
    
    plus_arr = arr[:,strands == "+"]
    plus_ys = ys[strands == "+"]
    minus_arr = arr[:, strands == "-"]
    minus_ys = ys[strands == "-"]

    lw = 100 / np.shape(arr)[1]

    subplot.plot([plus_arr[0,:], plus_arr[1,:]],
                 [plus_ys, plus_ys], color='red', lw=lw)

    if paired:
        subplot.plot([plus_arr[1,:], plus_arr[2,:]],
             [plus_ys, plus_ys], color='grey',
             ls='dotted', lw=lw)
        subplot.plot([plus_arr[2,:], plus_arr[3,:]],
                     [plus_ys, plus_ys],
                     color='blue', lw=lw)

    subplot.plot([minus_arr[0,:], minus_arr[1,:]],
                 [minus_ys, minus_ys], color='green', lw=lw)
    
    if paired:
        subplot.plot([minus_arr[1,:], minus_arr[2,:]],
                     [minus_ys, minus_ys], color='grey',
                     ls='dotted', lw=lw)
        subplot.plot([minus_arr[2,:], minus_arr[3,:]],
                     [minus_ys, minus_ys], color='orange', lw=lw)

