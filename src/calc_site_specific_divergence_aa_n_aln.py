#! /usr/bin/python3

import sys
import argparse
import math
from collections import Counter
from scipy.spatial import distance
from parse_fasta import parse_fasta


"""calculates site-specific Jensen-Shannon Divergences for n groups, given in
input as n different alignments corresponding to the same master alignment"""


def get_columns(seqDict: dict) -> dict:
    """takes a dictionary of an alignment (key: name, value: sequence),
    and returns a dictionary of columns (key: position, value: column)"""
    colDict = {}
    pos = 0
    for k, v in seqDict.items():
        for i in v:
            try:
                colDict[pos].append(i)
            except KeyError:
                colDict[pos] = []
                colDict[pos].append(i)
            pos += 1
        pos = 0
    return colDict


def calc_col_prop(colDict: dict) -> dict:
    """takes a column dictionary from get_columns and returns a dictionary
    of column amino acid state frequencies. Denominator does not count gaps."""
    colPropDict = {}
    aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F",
          "P", "S", "T", "W", "Y", "V"]
    for k, v in colDict.items():
        tot = sum([x[1] for x in Counter(v).items() if x[0] != "-"])
        if tot == 0:
            sys.stderr.write("skipping empty column " + str(k) + "\n")
            continue
        props = []
        for i in aa:
            try:
                c = Counter(v)[i]
            except KeyError:
                c = 0
            props.append(c / tot)
            colPropDict[k] = props
    return colPropDict


def kl(a: list, b: list) -> float:
    """where it follows that KL is defined only if for all x, Q(x) = 0
    implies P(x) = 0. If P(x) = 0 the corresponding term is taken as 0
    due to the limit"""
    div = [x * math.log2((x / y)) if x > 0 else 0 for x, y in zip(a, b)]
    return sum(div)


def jsd(a: list = [[1.0], [1.0], [1.0]]) -> float:
    """calculates JSD for n distributions, provided as a list of lists.
    By default the mixture is calculated with equal weights.

    note that because the mixture is calculated, the case where P(x) = 0
    when Q(x) > 0 is prevented from occurring, likewise the case of
    P(x) > 0 when Q(x) = 0"""
    m = [sum(x) / len(a) for x in zip(*a)]
    return sum([kl(x, m) / len(a) for x in a])


def calc_jsd(colPropDict1: dict, colPropDict2: dict) -> dict:
    """calculate jsd for two dictionaries of matching columns"""
    distDict = {}
    for k, v in colPropDict1.items():
        m = [(x + y) / 2 for x, y in zip(v, colPropDict2[k])]
        jsd = (kl(v, m) / 2) + (kl(colPropDict2[k], m) / 2)
        distDict[k] = jsd
    return distDict


def calc_jsd_scipy(colPropDict1: dict, colPropDict2: dict) -> dict:
    """test for scipy implementation"""
    distDict = {}
    for k, v in colPropDict1.items():
        distDict[k] = distance.jensenshannon(v, colPropDict2[k]) ** 2
    return distDict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("alignments", help="alignments for each group, \
                        columns matching", nargs="+")
    args = parser.parse_args()

    alns = {}
    for a in args.alignments:
        alns[args.alignments.index(a)] = dict([x for x in parse_fasta(a)])
        lens = []
        for v in alns.values():
            lens.append([len(seq) for seq in v.values()][0])
            # print(lens)
            if len(set(lens)) > 1:
                print("alignments are not of the same length!")
                sys.exit()

    alnsCols = {}
    for k, v in alns.items():
        alnsCols[k] = get_columns(v)

    colsProps = {}
    for k, v in alnsCols.items():
        colsProps[k] = calc_col_prop(v)

    cols = []
    for _, v in colsProps.items():
        cols.append(set([k for k in v.keys()]))
    inAll = sorted(list(set.intersection(*cols)))

    distances = {}

    for k in inAll:
        propsList = [x[k] for x in [v for v in colsProps.values()]]
        distances[k] = jsd(propsList)

    # for k, v in colsProps[0].items():
    #     propsList = [x[k] for x in [v for v in colsProps.values()]]
    #     distances[k] = jsd(propsList)

    print("pos\tjsd")
    for k, v in distances.items():
        print(str(k + 1) + "\t" + str(v))
