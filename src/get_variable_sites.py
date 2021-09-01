#! /usr/bin/python3

import sys
import argparse
from collections import Counter
from parse_fasta import parse_fasta


def get_site_dict(seqDict):
    siteDict = {}
    nSite = len([v for v in seqDict.values()][0])
    s = 0
    while s < nSite:
        states = [x[s] for x in seqDict.values()]
        siteDict[s] = states
        s += 1
    return siteDict


def get_variable(siteDict):
    var = []
    for k, v in siteDict.items():
        if len(set(v)) == 1:
            continue
        else:
            var.append(k + 1)
    return var


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help="FASTA-formatted alignment \
                        to extract variable sites")
    args = parser.parse_args()

    seqs = dict([x for x in parse_fasta(args.alignment)])
    sites = get_site_dict(seqs)
    # print(sites)
    variable = get_variable(sites)
    print(variable)

    for k, v in seqs.items():
        print(">" + k)
        outSeq = ""
        for i, j in enumerate(v):
            if i in variable:
                outSeq += j
        print(outSeq)
