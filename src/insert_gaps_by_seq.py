#! /usr/bin/python3

import sys
import argparse


def parse_fasta(path):  # courtesy of Jonathan Chang \
    # https://gist.github.com/jonchang/6471846
    """Given a path tries to parse a fasta file. Returns an iterator which
    yields a (name, sequence) tuple"""
    with open(path) as handle:
        name = sequence = ""
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, sequence
                name = line[1:]
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    ap = argparse.ArgumentParser()
    ap.add_argument("-a", "--alignment", help="Empirical alignment, \
                    in FASTA format.")
    ap.add_argument("-s", "--simulated", help="Simulated alignment, \
                    in FASTA format.")
    args = ap.parse_args()

    emp_dict = dict([x for x in parse_fasta(args.alignment)])
    sim_dict = dict([x for x in parse_fasta(args.simulated)])

    for k, v in sim_dict.items():
        for n, _ in enumerate(v):
            if emp_dict[k][n] == "-":
                sim_dict[k] = sim_dict[k][:n]+"-"+sim_dict[k][n+1:]

    for k, v in sim_dict.items():
        print(">"+k)
        print(v)
