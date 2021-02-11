#! /usr/bin/python3

import sys
import os
import argparse
from parse_fasta import parse_fasta
from get_altall_from_fastml import main


def read_node_file(path):
    """Read node file reads a tab delimited file with node
    numbers (e.g. N8) and desired names on one line:
    N8  giDODAa2_anc
    N10 amaDODAa1_anc
    etc.
    It returns a dictionary with node as key and name as
    value"""
    nodeDict = {}
    with open(path, "r") as inf:
        for line in inf:
            n = line.split("\t")[0].strip()
            name = line.split("\t")[1].strip()
            nodeDict[n] = name
    return nodeDict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("nodefile", help="tab delimited file with node and \
                        desired name on one line")
    parser.add_argument("fastml_dir", help="directory of fastml run to \
                        summarise")
    parser.add_argument("-gp", "--gapmode", help="Mode to include gaps. Options \
                    are [prob/fixed]. prob includes the alternate gap \
                    character if the PP is >= 0.2, while fixed treats \
                    the gap reconstruction as fixed, i.e. preferring the \
                    gap state if PP >= 0.5")
    args = parser.parse_args()

    fastmlDir = args.fastml_dir
    if fastmlDir[-1] != "/":
        fastmlDir += "/"
    nodeFile = args.nodefile

    nodes = read_node_file(nodeFile)
    probFile = fastmlDir + "prob.marginal.csv"
    gapFile = fastmlDir + "Ancestral_MaxMarginalProb_Char_Indel.txt"
    indelProb = fastmlDir + "IndelsMarginalProb.txt"
    seqFile = fastmlDir + "seq.marginal_IndelAndChars.txt"
    if args.gapmode is None:
        gapMode = "prob"
    else:
        gapMode = args.gapmode

    seqs = dict([x for x in parse_fasta(seqFile)])

    for k, v in nodes.items():
        mapSeq = seqs[k]
        altAllDict = main(probFile, gapFile, indelProb,
                          k, gapMode)
        altAllSeq = [x for x in altAllDict.values()][0]
        with open(v + ".pep.fa", "w") as outf:
            outf.write(">" + v + "_MAP\n")
            outf.write(mapSeq + "\n")
            outf.write(">" + v + "_AltAll\n")
            outf.write(altAllSeq + "\n")
