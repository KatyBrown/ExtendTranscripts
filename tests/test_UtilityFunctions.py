#!/usr/bin/env python3
import sys
import pytest
import numpy as np
sys.path.append("/home/katy/ExtendTranscripts/ExtendTranscripts")
import UtilityFunctions
import unittest
import pickle


def test_FastaToDict():
    T = UtilityFunctions.FastaToDict("tests/test_fasta.fasta")
    assert T == {"Seq_1": "ACGT",
                 "Seq_2": "AGCT",
                 "Seq_3": "CTAG"}


def test_AlignmentArray():
    F = UtilityFunctions.FastaToDict("tests/test_fasta.fasta")
    ali = UtilityFunctions.AlignmentArray(F.values())
    comp = [['A', 'C', 'G', 'T'], ['A', 'G', 'C', 'T'], ['C', 'T', 'A', 'G']]
    comp = np.array(comp)
    assert np.array_equal(comp, ali)


def test_IUPAC():
    P1, P2 = pickle.load(open("tests/IUPAC_test", "rb"))
    D1, D2 = UtilityFunctions.IUPAC()
    assert P1 == D1 and P2 == D2


def test_getRCDict():
    P1 = pickle.load(open("tests/rc_test", "rb"))
    D1 = UtilityFunctions.getRCDict()
    assert P1 == D1


def test_reverseComplement():
    rseq = "ACTGNYKRMBVDHWS-"
    assert UtilityFunctions.reverseComplement("-SWDHBVKYMRNCAGT") == rseq


if __name__ == '__main__':
    pytest.main()
