#!/usr/bin/env python3


"""
calculates the effective number of amino acids per site from a FASTA-formatted alignment as a
measure of site conservation
"""


import sys
import argparse
import math
from parse_fasta import parse_fasta
from calc_site_specific_divergence_aa_n_aln import get_columns, calc_col_prop


def calc_n_eff(freqs):
    """
    calculates the effective number of amino acids from a vector of amino acid frequencies
    """
    H = -sum((x * math.log(x)) for x in freqs if x > 0)
    n_eff = math.exp(H)
    return n_eff


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("aln", help = "FASTA formatted Amino Acid alignment")
    args = parser.parse_args()

    seqs = dict(parse_fasta(args.aln))

    cols = get_columns(seqs)

    col_props = calc_col_prop(cols)

    n_eff_dict = {}

    for k, v in col_props.items():
        n_eff_dict[k] = calc_n_eff(v)

    print("pos\tn_eff")
    for k, v in n_eff_dict.items():
        # 1-indexed
        print(str(k + 1) + "\t" + str(v))
