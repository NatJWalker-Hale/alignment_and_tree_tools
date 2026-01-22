#! /usr/bin/env python3


import sys
import argparse
from collections import Counter
import sequence as sq


def filter_singletons(col_dict: dict) -> dict:
    out_dict = {}
    for header, col in col_dict.items():
        state_count = Counter(c for c in col.values() if c.upper() != "N")
        if len(state_count) == 2 and any(v == 1 for v in state_count.values()):
            continue
        out_dict |= {header: col}
    return out_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('aln', help='FASTA formatted alignment to filter')
    args = parser.parse_args(sys.argv[1:] or ['--help'])

    seq_dict = dict(sq.parse_fasta(args.aln))
    col_dict = sq.get_columns(seq_dict)

    out_dict = sq.col_dict_to_seq_dict(filter_singletons(col_dict))

    for header, seq in out_dict.items():
        print(f">{header}", file=sys.stdout)
        print(f"{seq}", file=sys.stdout)


if __name__ == "__main__":
    main()
