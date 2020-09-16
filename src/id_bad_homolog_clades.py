#! /usr/bin/env python3

import sys
import os
import argparse
import copy
from ete3 import Tree, PhyloTree
import numpy as np


def calc_trlen(t):
    trlen = 0
    for n in tr.iter_descendants(strategy="postorder"):
        trlen += n.dist
    return trlen


def calc_sub_trlen(t):
    sub_trlen = []
    for n in tr.iter_descendants(strategy="postorder"):
        if not n.is_leaf():
            tmp_trlen = 0
            for c in n.iter_descendants(strategy="postorder"):
                tmp_trlen += c.dist
            sub_trlen.append((n, tmp_trlen, tmp_trlen / len(n)))
    return sub_trlen


def resample(sub_trlen_l, niter):
    resample_succ = {}
    i = 0
    nnode = len(sub_trlen_l)
    tmp_l = [x[2] for x in sub_trlen_l]
    for n in sub_trlen_l:
        resample_succ[n[0]] = 0
    while i < niter:
        rs = np.random.choice(tmp_l, size=nnode, replace=False)
        q95 = np.quantile(rs, q=0.95)
        for n in sub_trlen_l:
            if n[2] > q95:
                resample_succ[n[0]] += 1
        i += 1
    for key, value in resample_succ.items():
        resample_succ[key] = value / niter
    return resample_succ


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    ap = argparse.ArgumentParser()
    ap.add_argument("-t", "--tree", help="Newick treefile to check.")
    ap.add_argument("-og", "--outgroupf",
                    help="File of outgroup taxa in tree (one per line).")
    ap.add_argument("-it", "--iterate", help="Number of resamples [1000]",
                    type=int, default=1000)
    args = ap.parse_args()

    og_list = []
    with open(args.outgroupf, "r") as ogf:
        for line in ogf:
            og_list.append(line.strip())

    tr = PhyloTree(args.tree,
                   sp_naming_function=lambda node: node.name.split("@")[0])

    og_in_tr = []
    for l in tr.iter_leaves():
        if l.name.split("@")[0] in og_list:
            og_in_tr.append(l.name)

    all_l = [b.name for b in tr.iter_leaves()]
    ing = list(set(all_l) - set(og_in_tr))

    tr.set_outgroup(tr.get_common_ancestor(*og_in_tr))
    tr.prune(ing, preserve_branch_length=True)
    all_l = list(set(all_l) - set(og_in_tr))

    trlen = calc_trlen(tr)

    sub_trlen = calc_sub_trlen(tr)

    resamp_dict = resample(sub_trlen, args.iterate)

    print([(k, v) for k, v in sorted(resamp_dict.items(), key=lambda x: x[1])][-9][0])