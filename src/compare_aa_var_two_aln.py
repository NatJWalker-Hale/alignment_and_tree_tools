#! /usr/bin/python3


import sys
import argparse
import math
from parse_fasta import parse_fasta
from calc_site_specific_divergence_aa_n_aln import get_columns
from calc_site_specific_divergence_aa_n_aln import calc_col_prop
from calc_site_specific_divergence_aa_n_aln import calc_jsd


def calc_n_eff(freqs):
    H = -sum([x * math.log(x) for x in freqs if x > 0])
    n_eff = math.exp(H)
    return n_eff


if __name__ == "__main__":
    if len(sys.argv[0:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("aln1")
    parser.add_argument("aln2")
    args = parser.parse_args()

    seqs1 = dict([x for x in parse_fasta(args.aln1)])
    seqs2 = dict([x for x in parse_fasta(args.aln2)])

    cols1 = get_columns(seqs1)
    cols2 = get_columns(seqs2)

    props1 = calc_col_prop(cols1)
    props2 = calc_col_prop(cols2)

    nEffDict1 = {}
    for k, v in props1.items():
        nEffDict1[k] = calc_n_eff(v)

    nEffDict2 = {}
    for k, v in props2.items():
        nEffDict2[k] = calc_n_eff(v)

    diffNEff = {}
    for k, v in nEffDict1.items():
        diffNEff[k] = v - nEffDict2[k]

    divs = calc_jsd(props1, props2)
    diffNEffScaled = {}
    for k, v in diffNEff.items():
        diffNEffScaled[k] = v / 19

    print("pos\tdiff_in_n_eff")
    for k, v in diffNEff.items():
        print(str(k + 1) + "\t" + str(v))
