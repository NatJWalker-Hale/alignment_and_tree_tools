#! /usr/bin/python3


import sys
import argparse
import sequence as sq


def sliding_window(s: str, step: int) -> list:
    """
    get sliding windows from strings
    """
    start = 0
    out = []
    for i in range(0, len(s), step):
        wind = s[start:i+step]
        out.append(wind)
        start += step
    return out


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("FASTA", help="input FASTA-formatted sequences")
    parser.add_argument("window", help="window size", type=int)
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    seq_dict = dict(sq.parse_fasta(args.FASTA))

    for header, seq in seq_dict.items():
        windows = sliding_window(seq, args.window)
        print(windows)
        count = 1
        for w in windows:
            with open(args.FASTA+f".{count}", "a", encoding="utf-8") as outf:
                outf.write(f">{header}\n")
                outf.write(f"{w}\n")
            count += 1

