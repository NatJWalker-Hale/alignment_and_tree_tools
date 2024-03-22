#! /usr/bin/python3

import sys
import argparse
from parse_fasta import parse_fasta


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
            try:
                if emp_dict[k][n] == "-":
                    sim_dict[k] = sim_dict[k][:n]+"-"+sim_dict[k][n+1:]
            except KeyError:
                sys.stderr.write(k + " not in simulated data, skipping\n")
                
    for k, v in sim_dict.items():
        print(">"+k)
        print(v)
