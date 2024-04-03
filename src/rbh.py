#!/usr/bin/env python3


"""
Script to call Reciprocal Best Hits from 1-2 2-1 blast results in outfmt6
"""


import sys
import argparse
from typing import Iterator


def get_rbh(in1: str, in2: str) -> list[tuple[str, str]]:
    """
    returns a list of reciprocal best hits between query-reference blast in in1 and
    reference-query blast in in2
    """
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


def get_bh(in1: str) -> Iterator[tuple]:
    """
    yields a tuple of top hit per query in reference-query blast results in in1
    """
    res1 = {}
    with open(in1, "r", encoding="utf-8") as inf:
        for line in inf:
            line = line.strip().split()
            if line[0] in res1:
                continue
            res1[line[0]] = line[1:]
            yield line[0], line[1]


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-rbh", "--reciprocal_best_hits", help="query-reference and \
                        reference-query blast results in outfmt6", type=str, nargs=2)
    group.add_argument("-bh", "--best_hits", help="query-reference blast results in outfmt6",
                       type=str)
    args = parser.parse_args()

    if args.reciprocal_best_hits is not None:
        f1, f2 = args.reciprocal_best_hits
        rbhs = get_rbh(f1, f2)
        for i in rbhs:
            print(f"{i[0]}\t{i[1]}")
    
    if args.best_hits is not None:
        for i in get_bh(args.best_hits):
            print(f"{i[0]}\t{i[1]}")