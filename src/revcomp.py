#! /usr/bin/python3


"""
Script to rev or revcomp a sequence or specific substring of a sequence
"""


import sys
import argparse
import sequence as sq


def rev(s: str) -> str:
    """
    return a reversed string
    """
    return s[::-1]


def revcomp(s: str) -> str:
    """
    return a reverse-complemented sequence
    """
    corres = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return s.translate(s.maketrans(corres))[::-1]


def transform(s: str, r = False) -> str:
    if not r:
        return revcomp(s)
    else:
        return rev(s)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("FASTA", help="input FASTA", type=str)
    parser.add_argument("seq", help="sequence ID in FASTA to process", type=str)
    parser.add_argument("-r", "--rev", help="reverse. Otherwise revcomp", action="store_true")
    parser.add_argument("-s", "--start", help="start position. If no end specified, will \
                        go all the way to sequence end", type=int)
    parser.add_argument("-e", "--end", help="end position", type=int)
    args = parser.parse_args(sys.argv[1:] or ["--help"])


    seqs = dict(sq.parse_fasta(args.FASTA))

    if not args.start:
        start = 1
    else:
        start = args.start

    try:
        target = seqs[args.seq]
    except KeyError:
        sys.stderr.write(f"sequence {args.seq} not in {args.FASTA}\n")
        sys.exit()
    if args.end:
        out = (target[:start - 1] +
               transform(target[start - 1:args.end], r = args.rev) +
               target[args.end:])
    else:
        out = (target[:start - 1] + 
                transform(target[start - 1:], r = args.rev))
    sys.stdout.write(f">{args.seq}\n")
    sys.stdout.write(f"{out}\n")