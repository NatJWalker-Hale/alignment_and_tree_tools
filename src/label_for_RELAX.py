#! /usr/bin/python3


"""
script to label trees for RELAX for the analysis of ant-farming mutualism in P. nagasau
"""


import sys
import argparse
import newick3


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to label trees with P. nagasau as \
                                     foreground for RELAX")
    parser.add_argument("tree", help="tree topology to label. Can have branch lengths. Internal \
                        nodes will be overwritten if present")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-r", "--relax", help="hyphy-style labels for RELAX", action="store_const",
                       dest="type", const="relax", default="relax")
    group.add_argument("-p", "--paml", help="codeml-style labels for branch models", dest="type",
                       action="store_const", const="paml")
    parser.add_argument("-b", "--brlen", help="write branch lengths", action="store_true")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    curroot = newick3.parse_from_file(args.tree)

    for n in curroot.iternodes(order=0):
        if not args.brlen:
            n.length = None
        if n.parent is None:
            continue
        if n.label == "PH01S":
            if args.type == "relax":
                n.label += "{Test}"
            else:
                n.label += "#1"
        else:
            if args.type == "relax":
                if n.istip:
                    n.label += "{Reference}"
                else:
                    n.label = "{Reference}"
            else:
                if not n.istip:
                    n.label = None

    print(newick3.to_string(curroot) + ";")
