#! /usr/bin/python3

import sys
import argparse
from parse_fasta import parse_fasta
from ete3 import Tree


def get_node_label(t, tip1, tip2, tip3):
    n = t.get_common_ancestor(tip1, tip2, tip3)
    lab = n.name
    return lab


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tree", help="newick tree file with ancestral \
                        node labels")
    parser.add_argument("-s", "--sequence", help="fasta file of ancestral \
                        sequences with names corresponding to nodes")
    parser.add_argument("-og", "--outgroup", help="tip(s) to define outgroup \
                        for rooting, space separated", nargs="+")
    parser.add_argument("-tp", "--tips", help="tips (3) to define ancestral \
                        node, space separated", nargs="+")
    args = parser.parse_args()

    og = args.outgroup
    tips = args.tips
    seqs = dict([x for x in parse_fasta(args.sequence)])
    t = Tree(args.tree, format=1)
    t.set_outgroup(t.get_common_ancestor(*og))
    name = get_node_label(t, *tips)
    print(">"+name)
    print(seqs[name])
