#! /usr/bin/python3

import sys
import argparse
import tree_reader


def parse_og(path):
    with open(path, "r") as ogfile:
        og = [line.strip() for line in ogfile.readlines()]
    return og


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("og", help="Outgroup file, one taxon per line")
    parser.add_argument("tree", help="Tree with outgroups to remove")
    args = parser.parse_args()

    with open(args.tree, "r") as inf:
        nwkString = inf.readline().strip()
    OGs = parse_og(args.og)
    curroot = tree_reader.read_tree_string(nwkString)
    OGseqs = []
    for i in curroot.lvsnms():
        if i.split("@")[0] in OGs:
            OGseqs.append(i)
    for i in OGseqs:
        print(i)
