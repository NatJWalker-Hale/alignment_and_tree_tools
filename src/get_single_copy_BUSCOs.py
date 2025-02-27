#! /usr/bin/python3


import os
import sys
import argparse
import pandas as pd
from collections import Counter
import sequence as sq


def parse_results_table(path):
    df = pd.read_csv(path, sep="\t", skiprows=2)
    comp = df.loc[df["Status"] == "Complete"]["# Busco id"].tolist()
    return comp


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("workdir", help="working directory containing all BUSCO runs, in format \
                        busco_[tag], where [tag] is a label that will replace the name of this \
                        genome's single copy BUSCO sequence")
    parser.add_argument("odb", help="string identifying the odb used for the BUSCO run, e.g. \
                        embryophyta_odb10")
    parser.add_argument("prop", help="proportion of species required to have BUSCO copy for locus \
                        to be included, default 0 (no filter)", type=float)
    args = parser.parse_args(sys.argv[1:] or ["--help"])
    wd = args.workdir

    busco_outs = []
    for f in os.listdir(wd):
        if f.startswith("busco_"):
            busco_outs.append(f)

    # print(busco_outs)
        
    all_sc = []
    for d in busco_outs:
        all_sc += parse_results_table(f"{d}/run_{args.odb}/full_table.tsv")

    keep_sc = []
    for locus, count in Counter(all_sc).items():
        if count >= args.prop * len(busco_outs):
            keep_sc.append(locus)

    # print(keep_sc)
    # print(len(keep_sc))

    for locus in keep_sc:
        with open(f"{locus}.faa", "w", encoding="utf-8") as out_aaf:
            with open(f"{locus}.fna", "w", encoding="utf-8") as out_cdsf:
                for d in busco_outs:
                    tag = d.lstrip("busco_")
                    try:
                        aas = dict(
                            sq.parse_fasta(
                            (f"{d}/run_{args.odb}/busco_sequences/single_copy_busco_sequences/"
                            f"{locus}.faa")
                                        )
                                )
                        cds = dict(
                            sq.parse_fasta(
                            (f"{d}/run_{args.odb}/busco_sequences/single_copy_busco_sequences/"
                            f"{locus}.fna")
                                        )
                                )
                        out_aas = {tag: value for value in aas.values()}
                        out_cds = {tag: value for value in cds.values()}
                        out_aaf.write(sq.get_fasta_str(out_aas))
                        out_cdsf.write(sq.get_fasta_str(out_cds))
                    except FileNotFoundError:
                        continue
    # d = []
    # for _, dirs, _ in os.walk(wd):
    #     d.extend([x for x in dirs if x.startswith("run_")])
    #     break

    # buscoDict = {}
    # for rd in d:
    #     sp = rd.split("_")[1]
    #     tpath = wd + "/" + rd + "/run_embryophyta_odb10/full_table.tsv"
    #     sc = parse_results_table(tpath)
    #     buscoDict[sp] = sc

    # sets = [set(v) for _, v in buscoDict.items()]
    # oneToOne = sets[0].intersection(*sets[1:])

    # for orth in oneToOne:
    #     seqDict = {}
    #     for dirs in d:
    #         seqs = dict([x for x in parse_fasta(wd + "/" + dirs + "/translated_proteins/" + orth + ".faa")])
    #         seqDict[next(iter(seqs.items()))[0].split(" ")[0]] = next(iter(seqs.items()))[1]
    #     with open(wd + "/" + orth + ".faa", "w") as outfa:
    #         for k, v in seqDict.items():
    #             outfa.write(">"+k+"\n")
    #             outfa.write(v+"\n")
                    
    
    

