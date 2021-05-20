#! /usr/bin/python3

import sys
import os
import argparse
from parse_fasta import parse_fasta


def count_diff(seqDict, gaps=False):
    diffs = 0
    diffDict = {}
    seqNames = [x for x in seqDict.keys()]
    seqs = [x for x in seqDict.values()]
    if len(seqs[0]) != len(seqs[1]):
        print("sequences are not the same length, should be aligned")
        sys.exit()
    for i in range(0, len(seqs[0])):
        if not gaps:
            if seqs[0][i] == "-" or seqs[1][i] == "-":
                continue
            elif seqs[0][i] != seqs[1][i]:
                diffs += 1
                diffDict[i] = (seqs[0][i], seqs[1][i])
        else:
            if seqs[0][i] != seqs[1][i]:
                diffs += 1
                diffDict[i] = (seqs[0][i], seqs[1][i])
    return seqNames, diffs, diffDict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sequence", help="Pairwise alignment to count \
                        differences (FASTA)")
    parser.add_argument("-g", "--gapmode", help="Count gaps (default False)",
                        type=bool, default=False)
    args = parser.parse_args()

    seqs = dict([x for x in parse_fasta(args.sequence)])

    names, diffC, diffD = count_diff(seqs, args.gapmode)
    if args.gapmode:
        print("There are " + str(diffC) + " differences, counting gaps")
    else:
        print("There are " + str(diffC) + " differences, excluding gaps")
    print("")
    print("S1: "+names[0])
    print("S2: "+names[1])
    print("")
    print("pos\tS1\tS2")
    for k, v in sorted(diffD.items(), key=lambda x: x[0]):
        print(str(k+1) + "\t" + v[0] + "\t" + v[1])
