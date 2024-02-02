#! /usr/bin/python3


import sys
import argparse
import numpy as np
from collections import Counter


class Assembly:
    def __init__(self):
        self.file = ""
        self.length = 0
        self.n_contigs = 0
        self.lengths = []
        self.Ns = 0
        self.As = 0
        self.Cs = 0
        self.Gs = 0
        self.Ts = 0
        self.N50 = 0
        self.N90 = 0
        self.L50 = 0
        self.L90 = 0

    def print_stats(self):
        print(f"File: {self.file}")
        print(f"Sequences: {self.n_contigs}")
        print(f"Total length: {self.length}")
        print(f"Total length (without Ns): {self.length - self.Ns}")
        print(f"Average length: {self.length / self.n_contigs:.2f}")
        print(f"Shortest: {self.lengths[0]} Longest: {self.lengths[-1]}")
        print(f"As: {self.As} {self.As / self.length:.2%}")
        print(f"Cs: {self.Cs} {self.Cs / self.length:.2%}")
        print(f"Gs: {self.Gs} {self.Gs / self.length:.2%}")
        print(f"Ts: {self.Ts} {self.Ts / self.length:.2%}")
        print(f"Ns: {self.Ns} {self.Ns / self.length:.2%}")
        print(f"GC: {(self.Gs + self.Cs) / self.length:.2%}")
        print(f"N50: {self.N50} L50: {self.L50}")
        print(f"N90: {self.N90} L90: {self.L90}")


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("assembly", help="FASTA-formatted genome assembly")
    args = parser.parse_args()

    # seq_dict = dict([x for x in parse_fasta(args.assembly)])
    
    gen = Assembly()
    gen.file = args.assembly
    with open(args.assembly, "r") as f:
        length = 0
        line = f.readline().strip()  # here we read lines into mem one at a time
        while line:
            if line.startswith(">"):
                if length:
                    gen.lengths.append(length)
                    gen.length += length
                length = 0
                gen.n_contigs += 1
            else:
                length += len(line)
                counts = Counter(line.upper())
                gen.Ns += counts["N"]
                gen.As += counts["A"]
                gen.Cs += counts["C"]
                gen.Gs += counts["G"]
                gen.Ts += counts["T"]
            line = f.readline().strip()
        gen.lengths.append(length)
        gen.length += length

    # one way (slow, but works)
    # halflen = gen.length / 2
    # sums = 0
    # for i in sorted(gen.lengths):
    #     sums += i
    #     if sums >= halflen:
    #         gen.N50 = i
    #         break

    # using cumsum (still pretty slow tbh)
    gen.lengths = sorted(gen.lengths)
    csum = np.cumsum(gen.lengths)
    idx = np.where(csum == min(csum[csum >= int(gen.length / 2)]))
    # print(idx[0][0])
    gen.N50 = gen.lengths[idx[0][0]]
    gen.L50 = gen.n_contigs - idx[0][0]

    # n90 doing cumsum
    idx = np.where(csum == min(csum[csum >= int(gen.length / 10)]))
    gen.N90 = gen.lengths[idx[0][0]]
    gen.L90 = gen.n_contigs - idx[0][0]


    gen.print_stats()
