#!/usr/bin/env python3


"""
Script to call Reciprocal Best Hits from 1-2 2-1 blast results in outfmt6
"""


import sys
import argparse


def get_rbh(in1: str, in2: str) -> list:
    res1 = {}
    with open(in1, "r", encoding="utf-8") as inf:
        for line in inf:
            line = line.strip().split()
            if line[0] in res1:
                continue
            res1[line[0]] = line[1:]
    res2 = {}
    with open(in2, "r", encoding="utf-8") as inf:
        for line in inf:
            line = line.strip().split()
            if line[0] in res2:
                continue
            res2[line[0]] = line[1:]
    rbh = []
    for query, tophit in res1.items():
        try:
            if query == res2[tophit[0]][0]:
                rbh.append((query, tophit[0]))
        except KeyError:
            continue

    return rbh


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("blast_hits1")
    parser.add_argument("blast_hits2")
    args = parser.parse_args()

    rbhs = get_rbh(args.blast_hits1, args.blast_hits2)
    for i in rbhs:
        print(f"{i[0]}\t{i[1]}")