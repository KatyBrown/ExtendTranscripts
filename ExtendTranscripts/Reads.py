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
    min_kmers = endKmers(ref, kmer_length)
    results = []
    ref_rc = UtilityFunctions.reverseComplement(ref)
    with open(fastq_file) as infile:
        for i, line in enumerate(infile):
            if i % 4 == 0:
                nam = line.strip()

            if i % 4 == 1:
                read = line.strip()
                if min_kmers[0] in read or min_kmers[1] in read:
                    M = match(read, ref, kmer_length, nam, 'fwd')
                    if M[0]:
                        results.append(M)
                        continue
                if min_kmers[2] in read or min_kmers[3] in read:
                    M = match(read, ref_rc, kmer_length, nam, 'rev')
                    if M[0]:
                        results.append(M)
                        continue
    return (results)


def run(reference_file, fastq_file, kmer_length=8):
    F = readRef(reference_file)
    for nam in F.keys():
        seq = F[nam][:].seq
        results = iterateFastq(fastq_file, seq, kmer_length)
    return (results)

os.chdir("/home/katy/extend")
results = run("seq.fasta", "SRR4040080_final.fastq", 8)
