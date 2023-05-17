#! /usr/bin/python


import sys
import argparse
from parse_fasta import parse_fasta


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("seqs", help="FASTA-formatted sequences to return \
                        uppercase")
    args = parser.parse_args()

    seqDict = dict([x for x in parse_fasta(args.seqs)])

    for k, v in seqDict.items():
        print(">" + k)
        print(v.upper())
