import sys
import os
import argparse
from parse_fasta import parse_fasta


def count_diff(seqDict, gaps=False):
    diffs = 0
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
        else:
            if seqs[0][i] != seqs[1][i]:
                diffs += 1
    return diffs


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
    print(seqs)
    diffs = count_diff(seqs, args.gapmode)
    print(str(diffs))
