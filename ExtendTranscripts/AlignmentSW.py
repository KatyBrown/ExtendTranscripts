import UtilityFunctions
import Consensus
from UtilityFunctions import logPrint as lp
import copy
import Alignment


def runClusters(Z, fasta_dict, pD, seqdict, rround, cons=False,
                currentD=None, candidate=False, reference_dict=None):
    X = copy.copy(Z)
    if not cons:
        Znams = [z[0] for z in Z]
    else:
        Znams = []
        for z in Z:
            n = int(z[1].replace("*consensus_", ""))
            Znams += currentD[n]['names']
    done = set()
    k = 1
    D = dict()
    # keep building clusters until 1 or 0 sequences are left
    while len(done) < len(Znams) and len(X) > 1:
        current = dict()
        current['name'] = X[0][1]
        current['consensus'] = fasta_dict[X[0][1]][0:].seq
        lp("Remaining unaligned sequences %i" % (len(X)), 2, pD)

        # convert the query alignment into a matrix
        if cons:
            cons_n_q = int(current['name'].replace("*consensus_", ""))
            current['alignment'] = currentD[cons_n_q]['alignment']
            current['names'] = currentD[cons_n_q]['names']
        else:
            current['alignment'] = UtilityFunctions.AlignmentArray(
                [current['consensus']])
            current['names'] = [current['name']]

        current['matrix'], current['nt_inds'] = Consensus.makeAlignmentMatrix(
                                                current['alignment'])

        current['seqdict'] = seqdict
        # remove the query sequence from X
        X = X[1:]
        # Build a cluster based on the current query sequence
        X, C = buildCluster(X, fasta_dict, current, pD, k, cons, currentD,
                            candidate=candidate)

        D.setdefault(k, dict())

        lp("Saved cluster %i with %s sequences" % (k, len(C['current_names'])),
           1, pD)

        D[k]['consensus'] = C['current_consensus']
        D[k]['alignment'] = C['current_alignment']
        D[k]['names'] = C['current_names']
        D[k]['log'] = C['seqdict']
        Alignment.writeFastas(D[k], k, rround, pD['outdir'])
        done = done | set(C['current_names'])
        k += 1

    if len(X) == 1:
        D.setdefault(k, dict())
        D[k]['consensus'] = current['consensus']
        nam = X[0][1]
        if not cons:
            D[k]['alignment'] = UtilityFunctions.AlignmentArray(
                [current['consensus']])
            D[k]['names'] = [X[0][1]]
        else:
            cons_n_q = int(nam.replace("*consensus_", ""))
            D[k]['alignment'] = currentD[cons_n_q]['alignment']
            D[k]['names'] = currentD[cons_n_q]['names']
        D[k]['log'] = dict()
        lp("Saved cluster %i with %s sequences" % (k, len(D[k]['names'])),
           1, pD)
        Alignment.writeFastas(D[k], k, rround, pD['outdir'])
    return (D)



def buildCluster(X, fasta_dict, current, pD, k, cons=False, currentD=None,
                 candidate=False, skip=False):
    '''
    Build a cluster based on the current query sequence.
    Adapted for large input files - don't store everything in memory
    '''
    lp("Starting new cluster with %s" % (current['name']), 2, pD)
    # X is a list of names and integers
    j = 0
    # these are updated throughout to represent the current alignment
    matrix = current['matrix']
    nt_inds = current['nt_inds']
    seqdict = current['seqdict']

    current_names = current['names']
    current_consensus = current['consensus']
    current_alignment = current['alignment']

    # Each time a match is found and the the cluster is updated, start again
    # with the new cluster as
    # a query.
    i = 0
    # run a first pass without going back to the top every time a new sequence
    # is added - just switch to the consensus until you get to the end
    first_pass_done = False
    while j != len(X) and len(X) != 0:
        if first_pass_done:
            i = 0
        any_matches_inner = False
        n_new = 0
        # Look through all sequences which are not yet clustered until a match
        # to the current query is found (or you get to the end)
        while True and i != len(X):
            # update the query sequence
            query_seq = current_consensus
            query_ali = current_alignment
            query_names = current_names
            target_nam = X[i][1]
            target_seq = fasta_dict[target_nam][0:].seq
            if cons:
                cons_n_t = int(target_nam.replace("*consensus_", ""))
                target_ali = currentD[cons_n_t]['alignment']
                target_names = currentD[cons_n_t]['names']
            else:
                target_ali = UtilityFunctions.AlignmentArray([target_seq])
                target_names = [target_nam]
            lp("Testing %s" % ", ".join(target_names), 3, pD)
            # Align the query and target consensus sequences
            result = Alignment.SWalign(query_seq, target_seq,
                                       pD, useSub=True)

            # Check the if the consensus sequences are a good match
            is_match = Alignment.alignmentMeetsCriteria(result, query_seq,
                                                        target_seq, pD)
            # if they are not try the reverse complement

            if skip:
                skipnames = query_names
            else:
                skipnames = []
            if not is_match[0]:
                target_seq = UtilityFunctions.reverseComplement(target_seq)
                result = Alignment.SWalign(query_seq, target_seq,
                                           pD, useSub=False)
                is_match = Alignment.alignmentMeetsCriteria(result,
                                                            query_seq,
                                                            target_seq,
                                                            pD, candidate)
                target_ali = UtilityFunctions.reverseComplementAlignment(
                    target_ali)
                for nam in target_names:
                    seqdict[nam]['is_rc'] = True

            if is_match[0] and target_nam not in skipnames:
                lp("Match found.", 2, pD)
                # We found a match - something has changed
                any_matches_inner = True
                n_new += 1
                # remove the current value from X
                X = X[:i] + X[i+1:]
                result['alignment'] = is_match[1]
                # get the full alignment for the two consensus sequences
                result = Alignment.getAlignmentFull(result,
                                                    query_seq,
                                                    target_seq,
                                                    pD)
                current_names = query_names + target_names

                lp("Expanding current alignment to include %s" % (
                        ", ".join(target_names)), 3, pD)
                ali, matrix = Consensus.expandAlignment(result,
                                                        query_ali,
                                                        target_ali,
                                                        matrix,
                                                        nt_inds)
                # make a new sequence based on the new alignment
                current_consensus = Consensus.collapseAlignment(
                                    matrix, nt_inds)
                current_alignment = ali
                if first_pass_done:
                    i = 0
                # now you have a match and the consensus is updated,
                # start at the top again
                break
            else:
                lp("No match.", 3, pD)
            # keep going through the other sequences
            i += 1
        j += 1

        if any_matches_inner:
            # if anything has changed, clean up the alignment etc
            lp("Cluster %i updated - %s sequences" % (k,
                                                      len(current_names)),
               2, pD)
            lp("Cleaning cluster %i with CIAlign" % (k), 3, pD)
            if not candidate:
                funcs = ['remove_insertions', 'crop_ends', 'remove_gaponly']
            else:
                funcs = ['remove_gaponly']
            R = Alignment.cleanAlignmentCIAlign(current_alignment,
                                                current_names,
                                                query_seq,
                                                matrix,
                                                nt_inds,
                                                seqdict, pD, functions=funcs)

            current_alignment, matrix, current_consensus, seqdict = R
        elif not first_pass_done:
            first_pass_done = True
            i = 0
        else:
            break

    C = dict()
    C['current_alignment'] = current_alignment
    C['current_consensus'] = current_consensus
    C['current_names'] = current_names
    C['seqdict'] = seqdict
    C['matrix'] = matrix
    C['nt_inds'] = nt_inds
    return (X, C)
