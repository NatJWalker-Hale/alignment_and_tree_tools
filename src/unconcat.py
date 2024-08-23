#! /usr/bin/env python3


"""
splits a supermatrix into individual constituent FASTAs based on a RAxML-style partition file
"""


import sys
import argparse
import sequence as sq


def unconcat(# col_dict: dict,
             seq_dict: dict,
             part_dict: dict):
    """
    given a sequence dictionary and a dictionary of partitions, write a FASTA per
    partition
    """
    for name, sites in part_dict.items():
        # out_cols = {}
        # start = sites[0] - 1
        # end = sites[1] - 1
        # current = start
        # while (current >= start) & (current <= end):
        #     out_cols[current] = col_dict[current]
        #     current += 1
        # out_seq = sq.col_dict_to_seq_dict(out_cols)
        # sq.write_fasta(out_seq, f"{name}.fa")
        start = sites[0] - 1
        end = sites[1]
        with open(f"{name}.fa", "w", encoding="utf-8") as outfa:
            for header, seq in seq_dict.items():
                outfa.write(f">{header}\n")
                outfa.write(f"{seq[start:end]}\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="splits a supermatrix into individual constituent \
                                     FASTAs based on a RAxML-style partition file")
    parser.add_argument("supermatrix", help="FASTA-formatted input supermatrix")
    parser.add_argument("partition_file", help="RAxML-formatted partition file")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    seqs = dict(sq.parse_fasta(args.supermatrix))
    cols = sq.get_columns(seqs)

    part = dict(sq.parse_partition_file(args.partition_file))
    # unconcat(cols, part)
    unconcat(seqs, part)
