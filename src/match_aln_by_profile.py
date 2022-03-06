#! /usr/bin/python3

import sys
import argparse
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


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("a1", help="master alignment, in FASTA")
    parser.add_argument("a2", help="alignment to compare, in FASTA")
    parser.add_argument("pa", help="profile alignment of a1 and a2, \
                        in FASTA")
    args = parser.parse_args()

    aln1 = dict([x for x in parse_fasta(args.a1)])
    aln2 = dict([x for x in parse_fasta(args.a2)])
    proAln1 = dict([x for x in parse_fasta(args.pa)][:len(aln1)])
    proAln2 = dict([x for x in parse_fasta(args.pa)][len(aln1):])
    proAln1Sites = get_site_dict(proAln1)
    proAln2Sites = get_site_dict(proAln2)

    corres = []
    posAln1 = 1
    posAln2 = 1
    for i in range(len(proAln1Sites)):
        if all([x == "-" for x in proAln1Sites[i]]):  # gap in aln1
            corres.append((0, posAln2))
            posAln2 += 1
        else:  # char in aln1
            if all([x == "-" for x in proAln2Sites[i]]):  # gap in aln2
                corres.append((posAln1, 0))
            else:  # char in aln2
                corres.append((posAln1, posAln2))
                posAln2 += 1
            posAln1 += 1

    print("pos1\tpos2")
    for i in corres:
        print(str(i[0]) + "\t" + str(i[1]))
