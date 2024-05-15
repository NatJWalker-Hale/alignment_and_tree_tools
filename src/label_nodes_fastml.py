#! /usr/bin/python3


import sys
import argparse
import newick3


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="prints a newick with nodes labelled like FastML")
    parser.add_argument("tree", help="newick-formatted tree file - will overwrite internal node labels if present")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    curroot = newick3.parse_from_file(args.tree)

    count = 1
    for n in curroot.iternodes(order=0):
        if not n.istip:
            n.label = f"N{count}"
            count += 1

    print(newick3.to_string(curroot)+";")
