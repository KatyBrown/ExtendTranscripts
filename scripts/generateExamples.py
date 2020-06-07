#!/usr/bin/env python3
import numpy as np
import configargparse
import sys
import math

sys.path.append("../ExtendTranscripts")
import UtilityFunctions


def getNucFreqs(sequence):
    counts = np.unique(list(sequence), return_counts=True)
    freqs = dict(zip(counts[0], counts[1] / len(sequence)))
    return(freqs)


def randomBase(freqs):
    return np.random.choice(list(freqs.keys()), p=list(freqs.values()))


def randomPos(fragment, size):
    return (np.random.choice(range(len(fragment) - size)))


def minmax(minx, maxx):
    if int(minx) == minx and int(maxx) == maxx:
        return(np.random.choice(range(minx, maxx)))
    else:
        return(np.random.uniform(minx, maxx))


def addSNPs(fragment, freqs, min_diversity, max_diversity):
    f = list(fragment)
    diversity = minmax(min_diversity, max_diversity)
    ndiffs = math.ceil(diversity * len(fragment))
    for n in range(ndiffs):
        pos = randomPos(fragment, 1)
        f[pos] = randomBase(freqs)
    return ("".join(f))


def addIndels(fragment,
              freqs,
              min_n_indels,
              max_n_indels,
              min_indel_length,
              max_indel_length,
              typ='insertion'):
    f = list(fragment)
    n_indels = minmax(min_n_indels, max_n_indels)
    for i in range(n_indels):
        indel_length = minmax(min_indel_length, max_indel_length)
        pos = randomPos(fragment, indel_length)
        if typ == 'insertion':
            ins = [randomBase(freqs) for i in range(indel_length)]
            f = f[:pos] + ins + f[pos:]
        elif typ == 'deletion':
            f = f[:pos] + f[pos+indel_length:]
    return ("".join(f))


def chopSequence(sequence, nam,
                 min_n_fragments,
                 max_n_fragments,
                 min_fragment_length,
                 max_fragment_length):
    n_fragments = minmax(min_n_fragments, max_n_fragments)

    D = dict()
    for i in range(n_fragments):
        frag_nam = "%s_%s" % (nam.split(" ")[0], i+1)
        seq_len = minmax(min_fragment_length, max_fragment_length)
        seq_start = randomPos(sequence, seq_len)
        seq = sequence[seq_start:seq_start+seq_len]
        D[frag_nam] = seq
    return (D)


def main():
    parser = configargparse.ArgumentParser(
        description='''Generate pseudo fragmented transcripts with different \
            degress of variation and insertions / deletions''')
    parser.add_argument("--infile", dest="infile",
                        type=str,
                        help="path to input fasta file")
    parser.add_argument("--outfile", dest="outfile",
                        type=str,
                        help="path to output fasta file")
    parser.add_argument("--config", dest="config",
                        type=str,
                        help='path to config file',
                        is_config_file=True)

    parser.add_argument("--min_fragment_length", dest="min_fragment_length",
                        type=int,
                        help="minimum fragment length")
    parser.add_argument("--max_fragment_length", dest="max_fragment_length",
                        type=int,
                        help="maximum fragment length")

    parser.add_argument("--min_n_fragments", dest="min_n_fragments", type=int,
                        help="minimum number of fragments")
    parser.add_argument("--max_n_fragments", dest="max_n_fragments", type=int,
                        help="maximum number of fragments")

    parser.add_argument('--min_diversity', dest='min_diversity',
                        type=float,
                        help='minimum proportion of diversity to introduce')
    parser.add_argument('--max_diversity', dest='max_diversity',
                        type=float,
                        help='maximum proportion of diversity to introduce')

    parser.add_argument("--min_n_insertions", dest='min_n_insertions',
                        type=int,
                        help="minimum number of insertions to introduce")
    parser.add_argument("--max_n_insertions", dest='max_n_insertions',
                        type=int,
                        help="maximum number of insertions to introduce")

    parser.add_argument("--min_insertion_size", dest="min_insertion_size",
                        type=int,
                        help='minimum size insertion to introduce')
    parser.add_argument("--max_insertion_size", dest="max_insertion_size",
                        type=int,
                        help='maximum size insertion to introduce')

    parser.add_argument("--min_n_deletions", dest='min_n_deletions',
                        type=int,
                        help="minimum number of deletions to introduce")
    parser.add_argument("--max_n_deletions", dest='max_n_deletions',
                        type=int,
                        help="maximum number of deletions to introduce")

    parser.add_argument("--min_deletion_size", dest="min_deletion_size",
                        type=int,
                        help='minimum size deletion to introduce')
    parser.add_argument("--max_deletion_size", dest="max_deletion_size",
                        type=int,
                        help='maximum size deletion to introduce')

    args = parser.parse_args()

    F = UtilityFunctions.FastaToDict(args.infile)
    Fnew = dict()

    for nam, seq in F.items():
        freqs = getNucFreqs(seq)
        fragments = chopSequence(seq, nam,
                                 args.min_n_fragments,
                                 args.max_n_fragments,
                                 args.min_fragment_length,
                                 args.max_fragment_length)
        for fnam, fragment in fragments.items():
            fragment = addSNPs(fragment,
                               freqs,
                               args.min_diversity,
                               args.max_diversity)
            fragment = addIndels(fragment,
                                 freqs,
                                 args.min_n_insertions,
                                 args.max_n_insertions,
                                 args.min_insertion_size,
                                 args.max_insertion_size,
                                 typ='insertion')
            fragment = addIndels(fragment,
                                 freqs,
                                 args.min_n_deletions,
                                 args.max_n_deletions,
                                 args.min_deletion_size,
                                 args.max_deletion_size,
                                 typ='deletion')

            rc = np.random.choice([True, False])
            if rc:
                fragment = UtilityFunctions.reverseComplement(fragment)
            Fnew[fnam] = fragment
    out = open(args.outfile, "w")
    for key, val in Fnew.items():
        out.write(">%s\n%s\n" % (key, val))
    out.close()


if __name__ == "__main__":
    main()
