#! /usr/bin/python3

import sys
import argparse


def parse_blast_table(path):
    hitDict = {}
    with open(path, "r") as inf:
        for line in inf:
            fields = line.strip().split("\t")
            q = fields[0]  # query
            s = fields[2]  # subject
            e = fields[-2]  # e-value
            b = fields[-1]  # bitscore
            if q in hitDict:
                hitDict[q].append((s, e, b))
            else:
                hitDict[q] = [(s, e, b)]
    return hitDict


def get_top_hits(hitDict, eThresh, bitThresh):
    topHits = {}
    for k, v in hitDict.items():
        hitsSorted = sorted(v, key=lambda x: x[1])
        top = hitsSorted[0]
        if float(top[1]) <= eThresh and float(top[2]) >= bitThresh:
            topHits[k] = top[0]
        else:
            topHits[k] = ""
    return topHits


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("hits", help="BLAST output, outfmt 6, all columns")
    parser.add_argument("-e", "--evalue", help="E-value threshold, default \
                        0.001", type=float, default=0.001)
    parser.add_argument("-b", "--bitscore", help="bitscore threshold, default \
                        100", type=float, default=100)
    args = parser.parse_args()

    hits = parse_blast_table(args.hits)
    tops = get_top_hits(hits, args.evalue, args.bitscore)
    for k, v in tops.items():
        print(k+"\t"+v)

