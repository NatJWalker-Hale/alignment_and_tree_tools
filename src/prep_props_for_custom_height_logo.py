#! /usr/bin/python3

from ast import parse
import sys
import argparse
from collections import Counter
from parse_fasta import parse_fasta


def get_columns(seqDict):
    colDict = {}
    pos = 0
    for k, v in seqDict.items():
        for i in v:
            try:
                colDict[pos].append(i)
            except KeyError:
                colDict[pos] = []
                colDict[pos].append(i)
            pos += 1
        pos = 0
    return colDict


def calc_col_prop(colDict):
    colPropDict = {}
    aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F",
          "P", "S", "T", "W", "Y", "V"]
    for k, v in colDict.items():
        tot = sum([x[1] for x in Counter(v).items() if x[0] != "-"])
        props = []
        for i in aa:
            try:
                c = Counter(v)[i]
            except KeyError:
                c = 0
            props.append(c / tot)
            colPropDict[k] = props
    return colPropDict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("aln", help="FASTA alignment of sites to make logo")
    args = parser.parse_args()

    seqs = dict([x for x in parse_fasta(args.aln)])

    cols = get_columns(seqs)

    props = calc_col_prop(cols)

    for v in props.values():
        print(",".join([str(x) for x in v]))
