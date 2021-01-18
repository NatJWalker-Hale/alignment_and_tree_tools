#! /usr/bin/python3

import sys
import argparse
import itertools
from parse_fasta import parse_fasta


def get_pairwise_sequence_differences(seqDict):
    """for every pair of sequences in seqDict, populate a list of tuples
    (seq1, seq2, pos, state1, state2)"""
    diffList = []
    for i in itertools.combinations(seqDict, 2):
        for pos, state in enumerate(seqDict[i[0]]):
            if state == seqDict[i[1]][pos]:
                continue
            else:
                diffList.append((i[0], i[1], str(pos), state,
                                seqDict[i[1]][pos]))
    return diffList


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sequences", help="FASTA formatted file of \
                        aligned sequences to tabulate pairwise differences")
    args = parser.parse_args()

    seqs = dict([x for x in parse_fasta(args.sequences)])
    diffs = get_pairwise_sequence_differences(seqs)

    print("from,to,pos,from_state,to_state")
    for i in diffs:
        print(",".join(i))
