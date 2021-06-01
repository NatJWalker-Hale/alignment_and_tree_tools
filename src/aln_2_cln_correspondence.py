#! /usr/bin/python3

import sys
import argparse
from parse_fasta import parse_fasta
from collections import Counter


def check_aligned(seqDict):
    lens = []
    for v in seqDict.values():
        lens.append(len(v))
    if len(set(lens)) > 1:
        return False
    else:
        return True


def get_site_dict(seqDict):
    siteDict = {}
    nSite = len([v for v in seqDict.values()][0])
    s = 0
    while s < nSite:
        states = [x[s] for x in seqDict.values()]
        siteDict[s] = states
        s += 1
    return siteDict


def is_cleaned(siteList, prop):
    amb = ["X", "-", "?"]  # aa, do alph det later
    unamb = sum([v for k, v in Counter(siteList).items() if k not in amb])
    if unamb >= len(siteList) * prop:
        return False
    else:
        return True


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help="uncleaned alignment")
    parser.add_argument("proportion", type=float, default=0.5,
                        help="proportion of missing or ambiguous \
                              data to clean")
    args = parser.parse_args()

    seqs = dict([x for x in parse_fasta(args.alignment)])
    sites = get_site_dict(seqs)

    print("original\tclean")
    c = 0
    for k, v in sites.items():
        if is_cleaned(v, args.proportion):
            continue
        else:
            print(str(k+1) + "\t" + str(c+1))
            c += 1
