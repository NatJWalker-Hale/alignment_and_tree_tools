#! /usr/bin/python3

import sys
import argparse
import tree_reader
from parse_fasta import parse_fasta
from treenode import Node


def count_diffs(seq1, seq2, ignoreGap=True):
    """iterate over two aligned sequences and return differences
    as state1posstate2, e.g. A235V"""
    if len(seq1) != len(seq2):
        sys.stderr.write(seq1 + "not equal in length to" + seq2)
    rawDiffs = [(j[0], i+1, j[1]) for i, j in
                enumerate(zip(seq1, seq2)) if j[0] != j[1]]
    if ignoreGap:
        diffs = [j[0] + str(j[1]) + j[2] for j in rawDiffs if "-" not in j]
    else:
        diffs = [j[0] + str(j[1]) + j[2] for j in rawDiffs]
    return diffs


def get_anc_desc(node):
    """iterate over all nodes in a node-labelled tree and create
    a dictionary of ancestor descendant relationships, where each
    entry is one branch"""
    brDict = {}
    for n in node.iternodes(order="preorder"):
        for c in n.children:
            brDict[(n.label, c.label, c.length)] = []
    return brDict


def add_subs(brDict, seqDict, ignoreGap=True):
    for k in brDict.keys():
        brDict[k] += count_diffs(seqDict[k[0]], seqDict[k[1]], ignoreGap)


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="newick formatted tree, with ancestral \
                        node labels")
    parser.add_argument("sequences", help="FASTA formatted \
                        alignment, including ancestral \
                        sequences")
    parser.add_argument("-g", "--gaps", help="ignore gaps (default True)",
                        type=bool, default=True)
    args = parser.parse_args()

    with open(args.tree, "r") as t:
        for s in t:
            s = s.strip()
            nwkString = s

    curroot = tree_reader.read_tree_string(nwkString)
    branches = get_anc_desc(curroot)

    seqs = dict([x for x in parse_fasta(args.sequences)])
    # print(seqs)

    add_subs(branches, seqs, args.gaps)

    print("child\tparent\tlength")
    print(curroot.label)
    for k, v in branches.items():
        print(k[1] + "\t" + k[0] + "\t" + str(k[2]) + 2*"\t" + ",".join(v) +
              "\t" + ",".join(v))
