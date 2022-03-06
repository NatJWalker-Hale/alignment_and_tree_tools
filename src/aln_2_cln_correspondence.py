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
    parser.add_argument("-r", "--ref", help="reference sequence \
                        to find correspondence")
    args = parser.parse_args()

    seqs = dict([x for x in parse_fasta(args.alignment)])
    sites = get_site_dict(seqs)

    if args.ref is not None:
        isCleanedDict = {}

        for k, v in sites.items():
            if is_cleaned(v, args.proportion):
                isCleanedDict[k] = True
            else:
                isCleanedDict[k] = False
        # print(isCleanedDict)

        corrDict = {}
        alnPos = 0
        seqPos = 0
        for i in seqs[args.ref]:
            if i in ["X", "-", "?"]:
                corrDict[alnPos] = 0
                alnPos += 1
            else:
                corrDict[alnPos] = seqPos
                alnPos += 1
                seqPos += 1
        # print(corrDict)

        print("clean\tsequence")
        c = 0
        for k, v in corrDict.items():
            if isCleanedDict[k]:
                continue
            else:
                if v == 0:
                    print(str(c + 1) + "\t" + str(v))
                else:
                    print(str(c + 1) + "\t" + str(v + 1))
                c += 1
    else:
        print("original\tclean")
        c = 0
        for k, v in sites.items():
            if is_cleaned(v, args.proportion):
                continue
            else:
                print(str(k+1) + "\t" + str(c+1))
                c += 1
