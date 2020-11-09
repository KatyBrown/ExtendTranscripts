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
import re


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
            contig_seq,
            figdpi,
            intervals,
            coverage_lim,
            variants_g,
            common,
            minperc,
            covD,
            refD,
            altD,
            orfD):

    for interval in intervals:
        subtab = coverage_tab[(coverage_tab['pos'] >= interval[0]) &
                              (coverage_tab['pos'] <= interval[1])]
        nfunctions = 8
        subplots = dict()
        max_cov = max(subtab['coverage'])
        f = plt.figure(figsize=(15, 3*nfunctions),
                       dpi=figdpi)
        f = matplotlib.gridspec.GridSpec(nfunctions, 1)

        subplot_cov1 = plt.subplot(f[len(subplots), 0])
        plotbasicCov(subplot_cov1, interval, subtab,
                     colour='#DC143C', contig=contig, bam=bamnam)
        subplots['Raw Coverage'] = subplot_cov1

        subplot_cov2 = plt.subplot(f[len(subplots), 0])
        plotbasicCov(subplot_cov2, interval, subtab,
                     colour='#1f8aa7', contig=contig, bam=bamnam,
                     log=True)
        subplots['Log Coverage'] = subplot_cov2

        subplot_altcov = plt.subplot(f[len(subplots), 0])
        plotAltCov(subplot_altcov, altCov, interval, paired, rl)
        subplots['Read Position Ratios'] = subplot_altcov

        subplot_variants1 = plt.subplot(f[len(subplots), 0])
        plotVariants(subplot_variants1, variants_g, refD, altD, interval)
        subplots['All Variants > %.1f%% Frequency'
                 % (minperc * 100)] = subplot_variants1

        subplot_variants2 = plt.subplot(f[len(subplots), 0])
        plotVariants(subplot_variants2, variants_g, refD, altD, interval,
                     common_only=True, common=common)
        subplots['All Variants > 50.0% Frequency'] = subplot_variants2

        subplot_reads1 = plt.subplot(f[len(subplots), 0])
        subplot_reads2 = plt.subplot(f[len(subplots)+1, 0])

        subReadArr = getInRange(readArr[0], readArr[1],
                                readArr[2],
                                interval[0], interval[1],
                                rl)
        plotReads(subplot_reads1, subplot_reads2, subReadArr,
                  interval, paired, rl, coverage_lim, max_cov)
        subplots['Read Distribution + Strand'] = subplot_reads1
        subplots['Read Distribution - Strand'] = subplot_reads2

       # subplot_seq = plt.subplot(f[len(subplots), 0])
       # plotSeq(subplot_seq, contig_seq, interval, orfD)
       # subplots['Sequence'] = subplot_seq

        for i, subplot_nam in enumerate(subplots):
            subplot = subplots[subplot_nam]

            xstart = subplot.get_xlim()[0]
            xend = subplot.get_xlim()[1]
            tick_names, tick_positions, line_positions = getTicks(interval,
                                                                  xstart,
                                                                  xend)
            subplot.vlines(line_positions,
                           subplot.get_ylim()[0], subplot.get_ylim()[1],
                           lw=0.5, ls='dotted')
            if i+1 != len(subplots) and i != 0:
                # middle plots
                subplot.set_xticks(tick_positions)
                subplot.set_xticklabels([])

                subplot.xaxis.set_ticks_position('bottom')
            elif i == 0:
                # top plot
                subplot.xaxis.set_ticks_position('top')
                subplot.set_xticks(tick_positions)
                subplot.set_xticklabels(tick_names, fontsize=14)

            else:
                # bottom plot
                subplot.set_xticks(tick_positions)
                subplot.set_xticklabels(tick_names, fontsize=14)
                subplot.yaxis.set_ticks_position('left')
                subplot.xaxis.set_ticks_position('bottom')

            subplot.set_title(subplot_nam, loc='left', fontsize=14)
            subplot.yaxis.set_ticks_position('left')
            subplot.spines['right'].set_visible(False)
            subplot.spines['top'].set_visible(False)
        plt.tight_layout(pad=2)
        plt.savefig("%s/%s_coverage_%s_%s_%s.png" % (
                    outdir, contig, bamnam, interval[0], interval[1]),
                    dpi=figdpi, bbox_inches='tight')
        plt.close()


def getTicks(interval, start, end):
    interval_len = interval[1] - interval[0]

    om_len = math.floor(np.log10(interval_len))
    tick_int = int('5' + '0' * (om_len - 1))

    all_ticks = np.arange(0, math.ceil(interval[1]), tick_int)

    tick_names = all_ticks[all_ticks >= interval[0]]

    if interval[0] == start and interval[1] == end:
        line_positions = tick_names
    else:
        line_positions = tick_names - interval[0]

    tick_names = sorted(list(set(
        [interval[0]] + list(tick_names) + [interval[1]])))
    tick_positions = sorted(list(set(
        [start] + list(line_positions) + [end])))

    return (tick_names, tick_positions, line_positions)


def getInRange(nams, arr, strands, start_pos, end_pos, rl):
    which1 = np.sum(
        (arr[0:4, :] >= start_pos) & (
            arr[0:4, :] <= end_pos), 0) != 0
    which2 = (arr[0] <= start_pos) & (arr[-1] >= end_pos)
    which = which1 | which2
    whichnams = nams[which]
    whichstrand = strands[which]
    whicharr = arr[:, which]
    return(whichnams, whicharr, whichstrand)


def plotbasicCov(subplot, interval, coverage_tab, colour, contig, bam,
                 log=False):
    subplot.set_xlim(interval[0], interval[1])
    positions = coverage_tab['pos']
    coverage = coverage_tab['coverage']

    if max(coverage_tab['coverage'] != 0):
        subplot.plot(positions, coverage, color=colour, lw=3)
        if log:
            subplot.set_yscale('log')
        else:
            subplot.set_ylim(0, max(coverage) * 1.1)
        subplot.set_ylabel("Number of Reads")


def plotAltCov(subplot, altCov, interval, paired, rl):
    xrange = np.arange(interval[0], interval[1])
    starts, ends, bodies, inserts = altCov
    starts_sub = (starts + 1)[0:, interval[0]: interval[1]][0]
    ends_sub = (ends + 1)[0:, interval[0]: interval[1]][0]
    bodies_sub = (bodies + 1)[0:, interval[0]: interval[1]][0]
    # allows user to give a range which overhangs the end of the genome

    xrange = xrange[: len(starts_sub)]
    if paired:
        inserts_sub = (inserts + 1)[0:, interval[0]: interval[1]][0]

    subplot.scatter(xrange, starts_sub / bodies_sub, color='#D9514E',
                    marker='.', s=500)
    subplot.scatter(xrange, ends_sub / bodies_sub, color='#2DA8D8',
                    marker='.', s=500)

    if paired:
        subplot.scatter(xrange, inserts_sub / bodies_sub, color='#2A2B2D',
                        marker=".", s=500)

    subplot.set_xlim(interval[0], interval[1])


def plotReads(subplot1, subplot2,
              readArr, interval, paired, rl, coverage_lim, max_cov):
    colourD = {'+': dict(), '-': dict()}
    colourD['+']['read1'] = '#be0535'
    colourD['+']['read2'] = '#fcb738'
    colourD['+']['overlap'] = '#fc8d59'
    colourD['+']['gap'] = '#d3d2d4'

    colourD['-']['read1'] = '#3f09f5'
    colourD['-']['read2'] = '#1fa727'
    colourD['-']['overlap'] = '#05beaa'
    colourD['-']['gap'] = '#d3d2d4'

    subplots = [subplot1, subplot2]
    nams, arr, strands = readArr
    # ys = np.arange(0, np.shape(arr)[1])

    for i, (strand, subplot) in enumerate(zip(['+', '-'], subplots)):
        subplot.set_xlim(interval)
        if np.shape(arr)[1] != 0:
            strand_arr = arr[:, strands == strand]
            strand_ys = np.arange(np.shape(strand_arr)[1])

            if paired:
                strand_overlaps = (strand_arr[2, :] - strand_arr[1, :]) < 0
                strand_overlap_ys = strand_ys[strand_overlaps]
                strand_gaps = (strand_arr[2, :] - strand_arr[1, :]) >= 0
                strand_gap_ys = strand_ys[strand_gaps]
        else:
            continue

        if max_cov < coverage_lim:
            lw = 0.5

            subplot.plot([strand_arr[0, :], strand_arr[1, :]],
                         [strand_ys, strand_ys],
                         color=colourD[strand]['read1'], lw=lw)

            if paired:
                subplot.plot([strand_arr[2, :], strand_arr[3, :]],
                             [strand_ys, strand_ys],
                             color=colourD[strand]['read2'],
                             lw=lw)

                if len(np.shape(strand_overlaps)) == 2:

                    subplot.plot([strand_arr[1, :strand_overlaps],
                                  strand_arr[2, :strand_overlaps]],
                                 [strand_overlap_ys, strand_overlap_ys],
                                 colour=colourD[strand]['overlap'],
                                 lw=lw)

                if len(np.shape(strand_gaps)) == 2:
                    subplot.plot([strand_arr[1, :strand_gaps],
                                  strand_arr[2, :strand_gaps]],
                                 [strand_gap_ys, strand_gap_ys],
                                 color=colourD[strand]['gap'],
                                 lw=lw)

        else:
            breaks = np.array([strand_arr[0, x+1] < strand_arr[1, x]
                               for x in np.arange(
                                       np.shape(strand_arr)[1] - 1)])
            breaks = np.append(breaks, True)
            subplot.fill_betweenx(strand_ys,
                                  strand_arr[0, :], strand_arr[1, :],
                                  color=colourD[strand]['read1'],
                                  where=breaks)

            if paired:
                subplot.fill_betweenx(strand_ys,
                                      strand_arr[2, :], strand_arr[3, :],
                                      color=colourD[strand]['read2'])

                subplot.fill_betweenx(strand_gap_ys,
                                      strand_arr[1, strand_gaps],
                                      strand_arr[2, strand_gaps],
                                      color=colourD[strand]['gap'])

                subplot.fill_betweenx(strand_overlap_ys,
                                      strand_arr[1, strand_overlaps],
                                      strand_arr[2, strand_overlaps],
                                      color=colourD[strand]['overlap'])


def plotVariants(subplot, variants_g, refD, altD, interval,
                 common_only=False, common=None):
    colours = {'A': '#1ed30f',
               'G': '#f4d931',
               'T': '#f43131',
               'C': '#315af4'}
    nucs = ['A', 'C', 'T', 'G']
    subplot.set_xlim(0, interval[1] - interval[0])
    subplot.set_ylim(-0.1, 1.1)
    subplot.set_yticks([0, 1])

    for nuc in nucs:
        v = variants_g[nuc][interval[0]:interval[1]]
        if not common_only:
            w = np.where(np.sum(v, 1) != 0)[0]
        else:
            common_v = common[interval[0]:interval[1]]
            w = np.where((np.sum(v, 1) != 0) & common_v)[0]
        xs = np.vstack([w, w]).T
        ys = v[w, :]

        subplot.plot(xs.T, ys.T, color=colours[nuc], lw=3)
        for p in w:
            q = p + interval[0]
            if altD[q] == nuc:
                if v[p, 1] - v[p, 0] < 0.5:
                    subplot.text(p, 1.05, refD[q],
                                 ha='center', va='center',
                                 color=colours[refD[q]])
                else:
                    subplot.text(p, 1.05, altD[q],
                                 ha='center', va='center',
                                 color=colours[altD[q]])

                subplot.text(p, -0.05, refD[q],
                             ha='center', va='center',
                             color=colours[refD[q]])


def plotSeq(subplot, contig_seq, interval, orfD):
    frame_subs = {'S1':0, 'S2':1, 'S3':2, 'R1': 0, 'R2': 2, 'R3': 1}
    contig_subseq = contig_seq[interval[0]:interval[1]]
    colours = {'A': '#1ed30f',
               'G': '#f4d931',
               'T': '#f43131',
               'C': '#315af4'}
    fontsize = 1000 / len(contig_subseq)
    interval_len = interval[1] - interval[0]
    #subplot.set_xlim(interval[0], interval[1])

    if interval_len < 150:
        for i, char in enumerate(contig_subseq):
            if char not in colours:
                colour = 'lightgrey'
            else:
                colour = colours[char]
            subplot.text(i+0.5+interval[0],
                         1, char, fontsize=fontsize,
                         fontdict={'family': 'monospace', 'weight': 'bold'},
                         color=colour, ha='center', va='center')
    interval_start, interval_end = interval
    k = 2
    if interval_len < 450:
        for orf_nam in orfD:
            orf = orfD[orf_nam]
            orf_start = orf['span'][0]
            orf_end = orf['span'][1]
            orf_seq = orf['seq']
            orf_frame = orf['frame']
            print (orf_frame)
            orf_rf = frame_subs[orf_frame]
            if "R" in orf_frame:
                orf_seq = orf_seq[::-1]
                colour = 'blue'
            else:
                colour = 'red'
            orf_start_plot = max(orf_start, interval_start)
            orf_end_plot = min(orf_end, interval_end)
    
            if orf_start_plot > orf_start:
                clip_start = math.floor((orf_start_plot - orf_start) / 3)
            else:
                clip_start = 0
            if orf_end_plot < orf_end:
                clip_end = math.ceil((orf_end - orf_end_plot) / 3)
            else:
                clip_end = 0
        is_overlapping = (interval_end >= orf_start) and (orf_end >= interval_start)
        if is_overlapping:
            if clip_end == 0:
                orf_seq_plot = orf_seq[clip_start:]
            else:
                orf_seq_plot = orf_seq[clip_start:-clip_end]
            xp = orf_start_plot
            for char in orf_seq_plot:
                subplot.text(xp + orf_rf + 0.5, k+0.5, char, color=colour,
                             ha='center', va='center', fontsize=fontsize,
                             fontdict={'family': 'monospace',
                                       'weight': 'bold'})
                xp += 3
            
            if "R" in orf_frame:
                ticks = np.arange(orf_end_plot-orf_rf, orf_start_plot-orf_rf, -3)
                subplot.plot([orf_start_plot, orf_end_plot-orf_rf], [k, k], color=colour)
                subplot.plot([ticks, ticks], [k+2, k-2],  color=colour)

            k += 1
    subplot.set_ylim(0, k+1)

