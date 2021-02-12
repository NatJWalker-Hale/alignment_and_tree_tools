#! /usr/bin/python3

import sys
import argparse
from parse_fasta import parse_fasta





if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="newick formatted tree, with ancestral \
                        node labels")
    parser.add_argument("sequences", help="FASTA formatted \
                        alignment, including ancestral \
                        sequences")
    args = parser.parse_args()
