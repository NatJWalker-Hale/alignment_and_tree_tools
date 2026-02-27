#! /usr/bin/env python3


import sys
import argparse
import statistics
import newick3
from newick3 import Node


def brlen_var(tree: Node) -> float:
    """Calculate the variance of branch lengths across gene trees for a given branch of interest."""
    lens = [n.length for n in tree.iternodes()]
    if all(l == 0 for l in lens):
        raise ValueError("No branch lengths found in the tree.")
    return statistics.variance(lens)


def main():
    """main function"""
    parser = argparse.ArgumentParser(description="Calculate variance of branch lengths in a newick")
    parser.add_argument("trees", nargs="+", help="input tree file(s)")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    for tree_file in args.trees:
        curroot = newick3.parse_from_file(tree_file)
        try:
            print(f"{tree_file}\t{brlen_var(curroot)}")
        except ValueError as e:
            print(f"{tree_file}\t{e}", file=sys.stderr)


if __name__ == "__main__":
    main()
