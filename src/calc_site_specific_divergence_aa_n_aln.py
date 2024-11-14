#!/usr/bin/env python3


"""
calculates site-specific Jensen-Shannon Divergences for n groups, given in input as n different
alignments corresponding to the same master alignment
"""


import sys
import argparse
import math
from collections import Counter
from scipy.spatial import distance
from parse_fasta import parse_fasta


def get_columns(seq_dict: dict) -> dict:
    """takes a dictionary of an alignment (key: name, value: sequence),
    and returns a dictionary of columns (key: position, value: column)"""
    col_dict = {}
    pos = 0
    for v in seq_dict.values():
        for i in v:
            try:
                col_dict[pos].append(i)
            except KeyError:
                col_dict[pos] = []
                col_dict[pos].append(i)
            pos += 1
        pos = 0
    return col_dict


def calc_col_prop(col_dict: dict, atleast = 0.) -> dict:
    """takes a column dictionary from get_columns and returns a dictionary
    of column amino acid state frequencies. Denominator does not count gaps."""
    col_prop_dict = {}
    aa = "ARNDCQEGHILKMFPSTWYV"
    for pos, col in col_dict.items():
        counts = Counter(char for char in col if char in aa)
        tot = sum(counts.values())
        if atleast > 0:
            if tot < len(col) * atleast:
                sys.stderr.write(f"skipping column {pos} with insufficient data\n")
                continue
        else:
            if tot == 0:
                sys.stderr.write(f"skipping empty column {pos}\n")
                continue
        props = [counts[char] / tot for char in aa]
        col_prop_dict[pos] = props
    return col_prop_dict


def calc_kl(a: list, b: list) -> float:
    """where it follows that KL is defined only if for all x, Q(x) = 0
    implies P(x) = 0. If P(x) = 0 the corresponding term is taken as 0
    due to the limit"""
    div = [x * math.log2((x / y)) if x > 0 else 0 for x, y in zip(a, b)]
    return sum(div)


def calc_jsd(a: list=None) -> float:
    """calculates JSD for n distributions, provided as a list of lists.
    By default the mixture is calculated with equal weights.

    note that because the mixture is calculated, the case where P(x) = 0
    when Q(x) > 0 is prevented from occurring, likewise the case of
    P(x) > 0 when Q(x) = 0"""
    if a is None:
        a = [[1.0], [1.0], [1.0]]
    m = [sum(x) / len(a) for x in zip(*a)]
    return sum((calc_kl(x, m) / len(a)) for x in a)


def get_jsd(col_prop_dict1: dict, col_prop_dict2: dict) -> dict:
    """calculate jsd for two dictionaries of matching columns"""
    dist_dict = {}
    for pos, col in col_prop_dict1.items():
        m = [(x + y) / 2 for x, y in zip(col, col_prop_dict2[pos])]
        jsd = (calc_kl(pos, m) / 2) + (calc_kl(col_prop_dict2[pos], m) / 2)
        dist_dict[pos] = jsd
    return dist_dict


def calc_jsd_scipy(col_prop_dict1: dict, col_prop_dict2: dict) -> dict:
    """test for scipy implementation"""
    dist_dict = {}
    for pos, col in col_prop_dict1.items():
        dist_dict[pos] = distance.jensenshannon(col, col_prop_dict2[pos]) ** 2
    return dist_dict


def calc_site_specific_divergence_aa_n_aln(alns: dict) -> dict:
    """
    main
    """
    alns_cols = {}
    for fname, aln in alns.items():
        alns_cols[fname] = get_columns(aln)

    cols_props = {}
    for pos, col in alns_cols.items():
        cols_props[pos] = calc_col_prop(col, atleast=args.atleastp)

    cols = []
    for _, v in cols_props.items():
        cols.append(set(k for k in v.keys()))
    in_all = sorted(list(set.intersection(*cols)))

    distances = {}

    for k in in_all:
        props_list = [x[k] for x in list(cols_props.values())]
        distances[k] = calc_jsd(props_list)

    return distances


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("alignments", help="alignments for each group, \
                        columns matching", nargs="+")
    parser.add_argument("-ap", "--atleastp", help="proportion of characters \
                       required to be present in each group to calculate jsd \
                       (default 0, i.e. all non-empty columns will be used)",
                       type=float, default=0.0)
    args = parser.parse_args()

    alns = {}
    for aln_str in args.alignments:
        alns[args.alignments.index(aln_str)] = dict(parse_fasta(aln_str))
        lens = []
        for aln in alns.values():
            lens.append(set(len(seq) for seq in aln.values()).pop())
            if len(set(lens)) > 1:
                print("alignments are not of the same length!")
                sys.exit()

    out = calc_site_specific_divergence_aa_n_aln(alns)

    # for k, v in colsProps[0].items():
    #     propsList = [x[k] for x in [v for v in colsProps.values()]]
    #     distances[k] = jsd(propsList)

    print("pos\tjsd")
    for col, dist in out.items():
        print(str(col + 1) + "\t" + str(dist))
