#! /usr/bin/python3

import os
import sys
import argparse
import pandas as pd
from parse_fasta import parse_fasta


def parse_results_table(path):
    df = pd.read_csv(path, sep="\t", skiprows=[0, 1])
    comp = df.loc[df["Status"].isin(["Complete", "Fragmented"])]["# Busco id"].tolist()
    return comp


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("workdir", help="working directory containing all \
                        BUSCO runs, with prefix 'run_'")
    args = parser.parse_args()
    wd = args.workdir

    d = []
    for _, dirs, _ in os.walk(wd):
        d.extend([x for x in dirs if x.startswith("run_")])
        break

    buscoDict = {}
    for rd in d:
        sp = rd.split("_")[1]
        tpath = wd + "/" + rd + "/run_embryophyta_odb10/full_table.tsv"
        sc = parse_results_table(tpath)
        buscoDict[sp] = sc

    sets = [set(v) for _, v in buscoDict.items()]
    print(sets[0].intersection(*sets[1:]))
    
    

