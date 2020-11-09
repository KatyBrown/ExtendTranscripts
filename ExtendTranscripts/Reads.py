#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import random
import pyfaidx
import os
import numpy as np
import re
import networkx
import UtilityFunctions
import copy
import pandas as pd


def match(read, ref, overlap_length, nam, direc):
    i = 0
    M = [False]
    if direc == 'rev':
        read = UtilityFunctions.reverseComplement(read)
    for i in range(len(read), overlap_length, -1):
        ref_start = ref[:i]
        ref_end = ref[-i:]

        subread_f = read[:i]
        subread_r = read[-i:]
        if subread_f == ref_end:
            i += 1
            M = [True, nam, 'end', i-1, read]
            break
        elif subread_r == ref_start:
            i += 1
            M = [True, nam, 'start', i-1, read]
            break
        if M[0]:
            return (M)
    return (M)


def endKmers(seq, kmer_length):
    left = seq[:kmer_length]
    right = seq[-kmer_length:]
    return (left, right,
            UtilityFunctions.reverseComplement(right),
            UtilityFunctions.reverseComplement(left))


def readRef(reference_file):
    F = pyfaidx.Fasta(reference_file)
    return (F)


def iterateFastq(fastq_file, ref, kmer_length):
    newR = ref
    pR = ""
    k = 0
    allresults = pd.DataFrame()
    while pR != newR:
        pR = newR
        min_kmers = endKmers(newR, kmer_length)
        results = []
        ref_rc = UtilityFunctions.reverseComplement(newR)
        with open(fastq_file) as infile:
            for i, line in enumerate(infile):
                if i % 4 == 0:
                    nam = line.strip()
                if i % 4 == 1:
                    read = line.strip()
                    if min_kmers[0] in read or min_kmers[1] in read:
                        M = match(read, newR, kmer_length, nam, 'fwd')
                        if M[0]:
                            results.append(M)
                            continue
                    if min_kmers[2] in read or min_kmers[3] in read:
                        M = match(read, ref_rc, kmer_length, nam, 'rev')
                        if M[0]:
                            results.append(M)
                            continue
        results_t = pd.DataFrame(results)
        results_t = results_t.drop(0, 1)
        results_t.columns = ['read', 'site', 'length', 'seq']
        results_t['k'] = k
        allresults = allresults.append(results_t)
        starts = results_t[results_t['site'] == 'start']
        ends = results_t[results_t['site'] == 'end']
        starts['start_overhang'] = starts['seq'].str.split(
            min_kmers[0]).str.get(0)
        ends['end_overhang'] = ends['seq'].str.split(min_kmers[1]).str.get(1)
        starts = starts[starts['start_overhang'].str.len() > 0]
        ends = ends[ends['end_overhang'].str.len() > 0]

        for i in (0, 1):
            if i == 0:
                R = sorted(list(starts['start_overhang']),
                           key=lambda x: len(x))
            else:
                R = sorted(list(ends['end_overhang']),
                           key=lambda x: len(x))
            if len(R) != 0:
                lens = np.array([len(r) for r in R])
                seqs = [np.array(list(r)) for r in R]
                buffers = [np.array(["-"] * x) for x in max(lens) - lens]
                if i == 0:
                    M = np.array([np.append(b, s) for s, b in zip(
                        seqs, buffers)])
                else:
                    M = np.array([np.append(s, b) for s, b in zip(
                        seqs, buffers)])
                X = np.sum(M != "-", 0)
                nts = np.array(['A', 'G', 'C', 'T'])
                counts = np.array([np.sum(M == n, 0) for n in nts])
                seq = "".join(nts[np.argmax(counts, 0)])
                supp = np.max(counts, 0)
                freq = supp / X
                lims = np.where((freq < 0.9) | (supp < 100))[0]
                if len(lims) != 0:
                    if i == 0:
                        lim = lims[-1]
                        newR = "%s%s" % (seq[-lim:], newR)
                    else:
                        lim = lims[0]
                        newR = "%s%s" % (newR, seq[:lim])        
        k += 1
    return (newR, lens, seqs, buffers, X, M, allresults)
            
def run(reference_file, fastq_file, kmer_length=8):
    F = readRef(reference_file)
    for nam in F.keys():
        seq = F[nam][:].seq.upper()
        newR, lens, seqs, buffers, X, M, allresults = iterateFastq(fastq_file, seq, kmer_length)
    return (newR, lens, seqs, buffers, X, M, allresults)
