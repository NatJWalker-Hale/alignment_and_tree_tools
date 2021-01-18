#! /usr/bin/python3

import sys
import argparse
from ete3 import *


def parse_og(path):
    with open(path, "r") as ogfile:
        og = [line.strip() for line in ogfile.readlines()]
    return og


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("og", help="Outgroup file, one taxon per line")
    parser.add_argument("tree", help="Tree with monophyletic outgroups, for rooting")
    args = parser.parse_args()

    ogl = parse_og(args.og)

    t = PhyloTree(args.tree, format=1)

    ogseq = []
    for n in t.iter_leaves():
        if n.name.split("@")[0] in ogl:
            ogseq.append(n.name)

    try:
        t.check_monophyly(values=ogseq, target_attr="name")
    except TreeError:
        sys.stderr.write(str(args.tree) + "outgroup is non-monophyletic, exiting\n")
        sys.exit()

    anc = t.get_common_ancestor(*ogseq)
    try:
        t.set_outgroup(anc)
        print(t.write(format=1))
    except TreeError:
        sys.stderr.write(str(args.tree) + "outgroup already rooted, printing\n")
        print(t.write(format=1))


