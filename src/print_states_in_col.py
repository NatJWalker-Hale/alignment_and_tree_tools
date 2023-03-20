#! /usr/bin/python3

import sys
import argparse
from parse_fasta import parse_fasta


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help="FASTA formatted alignment")
    parser.add_argument("reference", help="Sequence to extract states")
    args = parser.parse_args()

    seqs = dict([x for x in parse_fasta(args.alignment)])

    try:
        for s in seqs[args.reference]:
            print(s)
    except KeyError:
        print("sequence is not in alignment")
        sys.exit()