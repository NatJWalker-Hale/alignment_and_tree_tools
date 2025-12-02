#! /usr/bin/python3

import os
import sys
import argparse
from statistics import variance
from collections import Counter
import newick3

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to prepare CAFE input from a directory of \
                                     homolog trees")
    parser.add_argument("indir", help="directory to search")
    parser.add_argument("treefile_ending", help="treefile suffix to inspect")
    args = parser.parse_args()

    files = []
    for f in os.listdir(args.indir):
        if f.endswith(args.treefile_ending):
            files.append(f)

    counts = {}
    taxon_set = set()
    for f in files:  # do first pass to get full set of taxa
        print(f"working on {f}", file=sys.stderr)
        curroot = newick3.parse_from_file(f)
        taxa = [l.label.split("@")[0] for l in curroot.leaves()]
        counts[f.split(".")[0]] = Counter(taxa)
        taxon_set.update(set(taxa))

    taxon_set = list(taxon_set)
    print(f"Desc\tFamily ID\t{"\t".join(taxon_set)}\tVariance", file=sys.stdout)
    for k, v in counts.items():
        print(f"(null)\t{k}\t{"\t".join(str(v[t]) for t in taxon_set)}\t"
              f"{variance([v[t] for t in taxon_set])}")



    