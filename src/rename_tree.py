#! /usr/bin/python

"""
script to relabel tree based on file of current names and new names
"""


import sys
import argparse
import newick3


def name_corres(curnames: str, newnames: str) -> dict:
    corres = {}
    with open(curnames, "r", encoding="utf-8") as f1, open(newnames, "r", encoding="utf-8") as f2:
        for line1, line2 in zip(f1, f2):
            corres[line1.strip()] = line2.strip()
    return corres


def relabel_tree(tree: newick3.Node, corres_dict: dict):
    for n in tree.iternodes():
        if n.istip:
            try:
                n.label = corres_dict[n.label]
            except KeyError:
                sys.stderr.write(f"{n.label} not in labels, skipping\n")
                continue


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Script to rename tree based on old names and new names")
    parser.add_argument("tree", help="newick-formatted tree file. Names will be overwritten")
    parser.add_argument("current_names", help="current labels in tree, one per line")
    parser.add_argument("new_names", help="new labels in tree, one per line. Order must match \
                        current_names")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    label_corres = name_corres(args.current_names, args.new_names)
    print(label_corres)

    curroot = newick3.parse_from_file(args.tree)

    relabel_tree(curroot, label_corres)

    print(newick3.to_string(curroot) + ";")

