#!/usr/bin/env python3
import numpy as np
import UtilityFunctions


def makeAlignmentMatrix(alignment):
    '''
    Make a matrix which will be used to record the details of the current
    alignment.  Each row is a nucleotide (including IUPAC nts) and each
    column is a position in the alignment.  The cells are the number of
    each nt in this column. This can then be used to make the consensus
    sequence.
    '''
    # import the IUPAC nts
    D1, D2 = UtilityFunctions.IUPAC()
    nts = list(D1.keys()) + ["-"]
    # make a dictionary to show which row in the matrix corresponds with
    # each nt
    nt_inds = dict((n, i) for i, n in enumerate(nts))
    # create an empty matrix with one row, plus one column for each nt in
    # the sequence
    sequence = alignment[0, :]
    alignment_tab = np.zeros((len(nts), len(sequence)))
    # fill in the matrix with the nts in the first sequence
    for i in range(np.shape(alignment)[0]):
        sequence = alignment[i, :]
        for col, char in enumerate(sequence):
            nt_ind = nt_inds[char]
            alignment_tab[nt_ind, col] += 1
    return (alignment_tab, nt_inds)


def getPosCigar(cigar, count_typs=['D', 'X']):
    '''
    Returns the positions at which insertions occur according to the cigar
    string and the number of insertions, used to insert gaps into the target
    sequence relative to the query sequence.
    '''
    # convert the cigar string into a list
    cigar_L = UtilityFunctions.readCIGAR(cigar)
    in_pos = []
    # n keeps track of the current position in the sequence
    n = 0
    for cig in cigar_L:
        typ = cig[-1]
        count = int(cig[:-1])
        if typ not in count_typs:
            n += count
        else:
            # append the position and the number of inserted nts
            in_pos.append((n, count))
            n += count
    return (in_pos)


def updateOneArray(array, cigar, direc='query'):
    if direc == 'query':
        indels = getPosCigar(cigar, count_typs=['D', 'X'])
    elif direc == "target":
        indels = getPosCigar(cigar, count_typs=['I', 'Y'])
    for pos, count in indels:
        dashes = np.chararray((np.shape(array)[0], count))
        dashes[:] = "-"
        newA_left = array[:, :pos]
        newA_right = array[:, pos:]
        array = np.concatenate((newA_left, dashes, newA_right), axis=1)
    return (array)


def updateArrays(old_alignment_array,
                 new_alignment_array,
                 matrix,
                 nt_inds,
                 cigar_updated):
    '''
    Updates the two arrays representing the alignment.
    Parameters
    ----------
    old_alignment_array: np.array
        Numpy array representing the previous MSA
    new_alignment_array: np.array
        Numpy array representing the new target sequence

    aligned with the consensus of the old MSA

    I = insertion in query
    D = insertion in target
    '''

    matrix_updated = matrix
    old_alignment_array_updated = old_alignment_array
    new_alignment_array_updated = new_alignment_array

    if cigar_updated == "%sM" % UtilityFunctions.lengthFromCIGAR(
            cigar_updated):
        full_alignment_array = np.concatenate((old_alignment_array,
                                               new_alignment_array))
    else:
        old_alignment_array_updated = updateOneArray(
                old_alignment_array_updated,
                cigar_updated,
                direc='query')
        new_alignment_array_updated = updateOneArray(
                                      new_alignment_array_updated,
                                      cigar_updated,
                                      direc='target')
    assert UtilityFunctions.lengthFromCIGAR(
            cigar_updated) == np.shape(
                    new_alignment_array_updated)[1], (
                            "Something is wrong with the target cigar string")
    assert UtilityFunctions.lengthFromCIGAR(
            cigar_updated) == np.shape(
                    old_alignment_array_updated)[1], (
                            "Something is wrong with the query cigar string")
    full_alignment_array = np.concatenate((old_alignment_array_updated,
                                          new_alignment_array_updated))
    matrix_updated = makeAlignmentMatrix(full_alignment_array)[0]
    return (full_alignment_array, matrix_updated)


def collapseAlignment(matrix, nt_inds):
    D1, D2 = UtilityFunctions.IUPAC()
    Z = sorted(zip(nt_inds.keys(), nt_inds.values()), key=lambda x: x[1])
    nt_keys = np.array([z[0] for z in Z])[:-1]
    nt_vals = np.array([z[1] for z in Z])[:-1]
    C = []
    for column in matrix.T:
        column = column[nt_vals]
        maxi = column == max(column)
        if sum(maxi) == 0 or sum(column) == 0:
            C.append("-")
        elif sum(maxi) == 1:
            C.append(nt_keys[maxi][0])
        else:
            maxi_nt = nt_keys[maxi]
            maxi_string = "".join(list(maxi_nt))

            if maxi_string not in D2:
                if maxi_string in ['A', 'C', 'G', 'T']:
                    C.append(maxi_string)
                else:
                    C.append("N")
            else:
                # there is >1 - use an IUPAC character
                iupac = D2[maxi_string]
                C.append(iupac)
    cons = "".join(C)
    return (cons)


def expandAlignment(result, query_alignment, target_alignment,
                    matrix, nt_inds):
    cigar_updated = result['cigar_updated']
    full_alignment_array, matrix = updateArrays(query_alignment,
                                                target_alignment,
                                                matrix,
                                                nt_inds,
                                                cigar_updated)

    return (full_alignment_array, matrix)
