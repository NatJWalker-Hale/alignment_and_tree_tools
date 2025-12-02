#! /usr/bin/python3


import os
import sys
import glob
import argparse
import pandas as pd
from collections import Counter
import sequence as sq


def parse_results_table(fpath: str, sc = False):
    df = pd.read_csv(fpath, sep="\t", skiprows=2)
    if sc:
        comp = df.loc[df["Status"] == "Complete"]["# Busco id"].tolist()

    else:
        comp = df.loc[df["Status"].isin(["Complete", "Duplicated"])]["# Busco id"].tolist()
    return list(set(comp))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("workdir", help="working directory containing all BUSCO runs, in format \
                        busco_[tag], where [tag] is a label that will replace the name of this \
                        genome's single copy BUSCO sequence")
    parser.add_argument("odb", help="string identifying the odb used for the BUSCO run, e.g. \
                        embryophyta_odb10")
    parser.add_argument("prop", help="proportion of species required to have BUSCO copy for locus \
                        to be included, default 0 (no filter)", type=float)
    parser.add_argument("-s", "--single_copy", help="collect only single-copy BUSCOs. Otherwise \
                        all overlapping BUSCOs will be kept", action="store_true")
    args = parser.parse_args(sys.argv[1:] or ["--help"])
    wd = args.workdir

    busco_outs = []
    for f in os.listdir(wd):
        if f.startswith("busco_"):
            busco_outs.append(f)

    # print(busco_outs)

    all = []
    for d in busco_outs:
        all += parse_results_table(f"{d}/run_{args.odb}/full_table.tsv", sc = args.single_copy)


    keep = []
    for locus, count in Counter(all).items():
        if count >= args.prop * len(busco_outs):
            keep.append(locus)

    # print(keep_sc)
    # print(len(keep_sc))

    for locus in keep:
        with open(f"{locus}.faa", "w", encoding="utf-8") as out_aaf:
            with open(f"{locus}.fna", "w", encoding="utf-8") as out_cdsf:
                for d in busco_outs:
                    tag = d.lstrip("busco_")
                    try:
                        aas = dict(
                            sq.parse_fasta(
                                (glob.glob(f"{d}/run_{args.odb}/busco_sequences/single_*/{locus}.faa") +
                                 glob.glob(f"{d}/run_{args.odb}/busco_sequences/multi_*/{locus}.faa"))[0]
                            )
                        )
                        cds = dict(
                            sq.parse_fasta(
                                (glob.glob(f"{d}/run_{args.odb}/busco_sequences/single_*/{locus}.fna") +
                                 glob.glob(f"{d}/run_{args.odb}/busco_sequences/multi_*/{locus}.fna"))[0]
                            )
                        )
                        count = 1
                        out_aas = {}
                        out_cds = {}
                        for key, value in aas.items():
                            out_aas[f"{tag}@{locus}_{count}"] = value
                            out_cds[f"{tag}@{locus}_{count}"] = cds[key]
                            count += 1
                        out_aaf.write(sq.get_fasta_str(out_aas))
                        out_cdsf.write(sq.get_fasta_str(out_cds))
                    except IndexError:
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
                    
    
    

