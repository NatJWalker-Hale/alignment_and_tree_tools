#!/usr/bin/env python3


"""
simple script to use pard to calculate various physico chemical distances at substituted sites
between two aligned proteins
"""


import sys
import argparse
from pard.grantham import grantham
from pard.sneath import sneath
from pard.experimental_exchangeability import experimental_exchangeability
from pard.koshi_goldstein import koshi_goldstein, MatrixType
import sequence as sq


def count_diffs(seq1: str, seq2: str, gaps=False) -> list:
    """iterate over two aligned sequences and return differences
    as state1posstate2, e.g. A235V"""
    if len(seq1) != len(seq2):
        raise ValueError("sequences not equal in length!")
    raw_diffs = [(j[0], i+1, j[1]) for i, j in
                enumerate(zip(seq1, seq2)) if j[0] != j[1]]
    if not gaps:
        diffs = [j for j in raw_diffs if "-" not in j]
    return diffs


if __name__ == "__main__":
    # if len(sys.argv[1:]) == 0:
    #     sys.argv.append("-h")

    parser = argparse.ArgumentParser(description="Calculates physicochemical and other distances \
                                     for substitutions between two sequences using pard. Grantham \
                                     and Sneath are distances, so higher is more different. \
                                     Experimental Exchangeability (Yamolsky and Stoltzfuss) and \
                                     Koshi and Goldstein are exchangeability matrices, so higher \
                                     is more similar")
    parser.add_argument("alignment", help="FASTA-formatted alignment of two sequences")
    parser.add_argument("-d", "--distance", default="g", choices=["g", "s", "ee", "kg"],
                        help="which distance measure to choose: (g)rantham, (s)neath, \
                        (e)xperimental (e)xchangeability (Yamolsky and Stoltzfuss), \
                        (k)oshi and (g)oldstein")
    args = parser.parse_args(sys.argv[1:] or ['--help'])

    aln = dict(sq.parse_fasta(args.alignment))

    if len(aln) > 2:
        sys.stderr.write(f"Found more than two sequences in {args.alignment}\n")
        sys.exit()

    seq0, seq1 = aln.values()
    aln_diffs = count_diffs(seq0, seq1)
    print(f"pos\t{args.distance}")
    for diff in aln_diffs:
        match args.distance:
            case "g":
                dist = grantham(diff[0], diff[2])
            case "s":
                dist = sneath(diff[0], diff[2])
            case "ee":
                dist = experimental_exchangeability(diff[0], diff[2],
                                                    False, warning=False)
            case "kg":
                dist = koshi_goldstein(diff[0], diff[2], matrix_type=MatrixType.ALL_RESIDUES,
                                       symmetric=False, warning=False)
        print(f"{diff[1]}\t{dist:.4f}")
