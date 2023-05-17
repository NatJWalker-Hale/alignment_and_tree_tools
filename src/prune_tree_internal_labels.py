#! /usr/bin/python3

import sys
import argparse
import tree_reader
from copy import deepcopy


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="newick tree to prune")
    parser.add_argument("labels", help="internal node labels \
                        to prune on, one per line")
    args = parser.parse_args()

    curroot = [i for i in tree_reader.read_tree_file_iter(args.tree)][0]

    labs = []
    with open(args.labels, "r") as inf:
        for line in inf:
            labs.append(line.strip())

    newTree = deepcopy(curroot)

    for n in curroot.iternodes(order="preorder"):
        if n.label in labs:
            for c in n.children:
                c.prune()

    print(curroot.get_newick_repr() + ";")