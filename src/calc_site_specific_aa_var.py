#!/usr/bin/env python3


"""
calculates the effective number of amino acids per site from a FASTA-formatted alignment as a
measure of site conservation
"""


import sys
import argparse
import math
from scipy.stats import entropy
from parse_fasta import parse_fasta
from calc_site_specific_divergence_aa_n_aln import get_columns, calc_col_prop


def calc_n_eff(freqs):
    """
    calculates the effective number of amino acids from a vector of amino acid frequencies
    """
    H = -sum((x * math.log(x)) for x in freqs if x > 0)
    n_eff = math.exp(H)
    return n_eff


def calc_n_eff_scipy(freqs):
    """
    calculates the effective number of amino acids from a vector of amino acid frequencies
    """
    H = entropy(freqs)
    n_eff = math.exp(H)
    return n_eff


def calc_scaled_entropy(freqs):
    """
    calculates 1 minus the log k-scaled Shannon entropy, such that fixed columns have value 1
    """
    H = -sum((x * math.log(x)) for x in freqs if x > 0)
    scaled_ent = abs(H / math.log(20))
    return 1-scaled_ent


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("aln", help = "FASTA formatted Amino Acid alignment")
    parser.add_argument("-ap", "--atleastp", help="proportion of characters \
                       required to be present in each group to calculate Neff \
                       (default 0, i.e. all non-empty columns will be used)",
                       type=float, default=0.0)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--neff", help="calculate the effective number of amino acids",
                       action="store_true")
    group.add_argument("--norm_ent", help="calculate the log k-scaled Shannon entropy",
                       action="store_true")
    args = parser.parse_args()

    seqs = dict(parse_fasta(args.aln))

    cols = get_columns(seqs)

    col_props = calc_col_prop(cols, atleast=args.atleastp)

    out_dict = {}

    for k, v in col_props.items():
        # try:
        #     calc_n_eff(v) == calc_n_eff_scipy(v)
        # except ValueError:
        #     sys.stderr.write("your code is wrong!")
        if args.neff:
            out_dict[k] = calc_n_eff(v)
        else:
            out_dict[k] = calc_scaled_entropy(v)

<<<<<<< HEAD
    if args.neff:
        print("pos\tn_eff")
    else:
        print("pos\tnorm_ent")
    for k, v in out_dict.items():
=======
    print("pos\tn_eff")
    for k, v in n_eff_dict.items():
        # 1-indexed
>>>>>>> 4b95d4df00b2df7c977bd7f8a661a7eca9b8711a
        print(str(k + 1) + "\t" + str(v))
