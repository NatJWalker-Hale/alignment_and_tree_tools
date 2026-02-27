#! /usr/bin/env python3


import sys
import argparse
import numpy as np
import newick3
from newick3 import Node


def brlen_var(tree: Node) -> float:
    """Calculate the variance of branch lengths across gene trees for a given branch of interest."""
    lens = [n.length for n in tree.iternodes()]
    if len(set(lens)) == 1 and lens[0] == 0:
        raise ValueError("No branch lengths found in the tree.")
    return np.var(lens, ddof=1)


def main():
    """main function"""
    parser = argparse.ArgumentParser(description="Calculate variance of branch lengths in a newick")
    parser.add_argument("tree", help="input tree with branch lengths in newick format")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    curroot = newick3.parse_from_file(args.tree)
    try:
        print(f"{args.tree}\t{brlen_var(curroot)}", file=sys.stdout)
    except ValueError as e:
        print(f"{args.tree}\t{e}", file=sys.stderr)


if __name__ == "__main__":
    main()
