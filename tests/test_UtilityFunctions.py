#!/usr/bin/env python3
import pytest
import numpy as np
import sys
sys.path.append("ExtendTranscripts")
import UtilityFunctions
import pickle


def test_FastaToDict():
    T = UtilityFunctions.FastaToDict("tests/test_fasta.fasta")
    assert T == {"Seq_1": "ACGTT-",
                 "Seq_2": "AACT",
                 "Seq_3": "CTAG"}


def test_AlignmentArray():
    F = UtilityFunctions.FastaToDict("tests/test_fasta.fasta")
    ali = UtilityFunctions.AlignmentArray(F.values())
    comp = [['A', 'C', 'G', 'T', 'T', "-"],
            ['A', 'A', 'C', 'T'],
            ['C', 'T', 'A', 'G']]
    comp = np.array(comp)
    assert np.array_equal(comp, ali)


def test_IUPAC():
    P1, P2 = pickle.load(open("tests/test_data/IUPAC_test", "rb"))
    D1, D2 = UtilityFunctions.IUPAC()
    assert P1 == D1 and P2 == D2


def test_getRCDict():
    P1 = pickle.load(open("tests/test_data/rc_test", "rb"))
    D1 = UtilityFunctions.getRCDict()
    assert P1 == D1


def test_reverseComplement():
    rseq = "ACTGNYKRMBVDHWS-"
    assert UtilityFunctions.reverseComplement("-SWDHBVKYMRNCAGT") == rseq


def test_reverseComplementAlignment():
    F = UtilityFunctions.FastaToDict("tests/test_data/test_fasta.fasta")
    ali = UtilityFunctions.AlignmentArray(F.values())
    rcomp = [['-', 'A', 'A', 'C', 'G', 'T'],
             ['A', 'G', 'T', 'T'],
             ['C', 'T', 'A', 'G']]
    rcomp = np.array(rcomp)
    assert np.array_equal(
            UtilityFunctions.reverseComplementAlignment(ali), rcomp)


def test_readCIGAR():
    cigar1 = '6M22X2Y12D6I'
    assert UtilityFunctions.readCIGAR(cigar1) == ['6M',
                                                  '22X',
                                                  '2Y',
                                                  '12D',
                                                  '6I']


def test_lengthFromCIGAR():
    cigar1 = '6M22X2Y12D6I'
    assert UtilityFunctions.lengthFromCIGAR(cigar1) == 48
    assert UtilityFunctions.lengthFromCIGAR(cigar1, excludeI=True) == 18
    assert UtilityFunctions.lengthFromCIGAR(cigar1, excludeD=True) == 12
    assert UtilityFunctions.lengthFromCIGAR(cigar1,
                                            excludeI=True,
                                            excludeD=True) == 6
    assert UtilityFunctions.lengthFromCIGAR(cigar1, mOnly=True) == 6


if __name__ == '__main__':
    pytest.main()
