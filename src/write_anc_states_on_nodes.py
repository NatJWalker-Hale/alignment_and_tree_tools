#! /usr/bin/python3


import sys
import argparse
import tree_reader
from copy import deepcopy
from parse_fasta import parse_fasta


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gokstad", help="write output for plotting \
                        coloured trees with gokstad", action="store_true")
    parser.add_argument("tree", help="newick tree, nodes labelled")
    parser.add_argument("aln", help="FASTA alignment, with ancestral \
                        sequences, names matching node labels")
    args = parser.parse_args()

    with open(args.tree, "r") as inf:
        nwkString = inf.readline().strip()
        curroot = tree_reader.read_tree_string(nwkString)

    labelledCurroot = deepcopy(curroot)

    seqs = dict([x for x in parse_fasta(args.aln)])

    if args.gokstad:
        sys.stderr.write("labelling for gokstad\n")

    with open("anc_state_tree.nwk", "w") as outf:
        pos = 0
        while pos < len(list(seqs.values())[0]):
            for n in curroot.iternodes():
                n.note = seqs[n.label][pos]
            labelledCurroot = deepcopy(curroot)
            if args.gokstad:
                for n in labelledCurroot.iternodes():
                    if n.istip:
                        n.label += "[&state=%s]" % n.note
                        n.note = ""
                    else:
                        n.label = "[&state=%s]" % n.note
                        n.note = ""
            else:
                for n in labelledCurroot.iternodes():
                    n.label = n.note
            outf.write(labelledCurroot.get_newick_repr() + ";" + "\n")
            pos += 1
