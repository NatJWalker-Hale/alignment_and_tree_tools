#! /usr/bin/python


import sys
import argparse
from parse_fasta import parse_fasta
from random import sample


def get_ident_seqs(seq_dict: dict) -> dict[list]:
    ident_dict = {}
    for i, j in seq_dict.items():
        if True in [i in v for v in ident_dict.values()]:
            continue
        for x, y in seq_dict.items():
            if (j.upper() == y.upper()) & (i != x):
                try:
                    ident_dict[i].append(x)
                except KeyError:
                    ident_dict[i] = [x]
    ident_list = []
    for k, v in ident_dict.items():
        ident_list.append([k] + v)
    return ident_list


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")


    parser = argparse.ArgumentParser()
    parser.add_argument("FASTA", help = "FASTA to find identical sequences")
    parser.add_argument("-p", "--pick",
                        help = "Randomly pick k-1 representatives of each identical group",
                        action = "store_true")
    parser.add_argument("-r", "--remove",
                        help = "Remove all k-1 randomly picked representatives and print filtered FASTA (implies -p)",
                        action = "store_true")
    args = parser.parse_args()


    seqs = dict([x for x in parse_fasta(args.FASTA)])


    ident_seqs = get_ident_seqs(seqs)


    if args.remove:
        remove = []
        for i in ident_seqs:
            for j in [k for k in sample(i, len(i) - 1)]:
                remove.append(j)
        for k, v in seqs.items():
            if k not in remove:
                print(">" + k)
                print(v)
    else:
        if args.pick:
            for i in ident_seqs:
                for j in [k for k in sample(i, len(i) - 1)]:
                    print(j)
        else:
            for i in ident_seqs:
                print("\t".join(i))