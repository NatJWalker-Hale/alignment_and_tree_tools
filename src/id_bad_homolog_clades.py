#! /usr/bin/env python3

import sys
import os
import argparse
import copy
from ete3 import Tree, PhyloTree
import numpy as np


def is_leaf(node):
    if node.is_leaf():
        return True
    else:
        return False


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    ap = argparse.ArgumentParser()
    ap.add_argument("-t", "--tree", help="Newick treefile to check.")
    ap.add_argument("-og", "--outgroupf",
                    help="File of outgroup taxa in tree (one per line).")
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

    all_l = [x.name for x in tr.iter_leaves()]
    ing = list(set(all_l) - set(og_in_tr))

    tr.set_outgroup(tr.get_common_ancestor(*og_in_tr))
    tr.prune(ing, preserve_branch_length=True)
    all_l = list(set(all_l) - set(og_in_tr))

    trlen = 0
    for n in tr.iter_descendants(strategy="postorder"):
        trlen += n.dist

    node_jackknife = []

    for i in tr.iter_descendants(strategy="postorder"):
        if not i.is_leaf():
            tmp_trlen = 0
            for j in tr.iter_descendants(strategy="postorder"):
                if j == i:
                    for k in j.iter_descendants(strategy="postorder"):
                        tmp_trlen += k.dist
                else:
                    continue
            node_jackknife.append((i, tmp_trlen, tmp_trlen / len(i),
                                   len(i)))

    sorted_node_jackknife = sorted(node_jackknife,
                                   key=lambda x: x[2])
    q1, q3 = np.percentile([x[2] for x in sorted_node_jackknife], [25, 75])
    iqr = q3 - q1
    upper_bound = q3 + (1.5 * iqr)

    outlier_tips = []
    outlier_node = [x[0] for x in sorted_node_jackknife if x[2] > upper_bound]
    for n in outlier_node:
        outlier_tips.extend([l.name for l in n.iter_leaves()])
    for s in set(outlier_tips):
        print(s)
