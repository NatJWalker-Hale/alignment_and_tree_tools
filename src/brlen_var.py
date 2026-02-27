#! /usr/bin/env python3


import sys
import argparse
import statistics
import newick3
from newick3 import Node


def brlen_var(tree: Node, internal_only: bool) -> float:
    """Calculate the variance of branch lengths across gene trees for a given branch of interest."""
    lens = [n.length for n in tree.iternodes() if not internal_only or not n.istip]
    if all(l == 0 for l in lens):
        raise ValueError("No branch lengths found in the tree.")
    return statistics.variance(lens)


def main():
    """main function"""
    parser = argparse.ArgumentParser(description="Calculate variance of branch lengths in a newick")
    parser.add_argument("trees", nargs="+", help="input tree file(s)")
    parser.add_argument("-i", "--internal-only", action="store_true",
                        help="only consider internal branches")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    for tree_file in args.trees:
        curroot = newick3.parse_from_file(tree_file)
        try:
            print(f"{tree_file}\t{brlen_var(curroot, args.internal_only)}")
        except ValueError as e:
            print(f"{tree_file}\t{e}", file=sys.stderr)


if __name__ == "__main__":
    main()
