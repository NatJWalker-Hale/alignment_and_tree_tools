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
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    curroot = newick3.parse_from_file(args.tree)

    for n in curroot.iternodes(order=0):
        if n.parent is None:
            continue
        if n.label == "PH01S":
            n.label += "{Test}"
        else:
            if n.istip:
                n.label += "{Reference}"
            else:
                n.label = "{Reference}"

    print(newick3.to_string(curroot) + ";")
    