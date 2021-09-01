#! /usr/bin/python3

import os
import sys
from parse_fasta import parse_fasta

# extracts alignment columns given in args.


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python " + sys.argv[0] + " aln col1 col2 col3 ...")
        sys.exit(0)

    seqDict = dict([x for x in parse_fasta(sys.argv[1])])

    collist = [int(i)-1 for i in sys.argv[2:]]
    # print collist
    for k in seqDict.keys():
        print(">"+k)
        print("".join([seqDict[k][i] for i in collist]))
