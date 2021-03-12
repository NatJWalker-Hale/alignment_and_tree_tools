#! /usr/bin/python3

import sys
import argparse
from parse_fasta import parse_fasta


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--positions", type=int, nargs="+",
                        help="optionally, a space=separated list \
                        of positions to correspond")
    parser.add_argument("alignment", help="alignment to get \
                        corresponding positions")
    parser.add_argument("ref_seq", help="reference sequence \
                        to get correspondence to")
    args = parser.parse_args()

    seqs = dict([x for x in parse_fasta(args.alignment)])

    if args.ref_seq not in seqs.keys():
        print("reference sequence is not in alignment")
        sys.exit()

    corrDict = {}
    aln_pos = 1
    seq_pos = 1
    for i in seqs[args.ref_seq]:
        if i == "-":
            corrDict[aln_pos] = 0
            aln_pos += 1
        else:
            corrDict[aln_pos] = seq_pos
            aln_pos += 1
            seq_pos += 1

    if args.positions is not None:
        for p in args.positions:
            print(corrDict[p])
    else:
        print("aln_pos\tref_pos")
        for k, v in corrDict.items():
            print(str(k)+"\t"+str(v))
