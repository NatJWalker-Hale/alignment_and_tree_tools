#! /usr/bin/python3

import sys
import argparse
import tree_reader


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("speciesTree", help="species tree in newick")
    parser.add_argument("geneTree", help="gene tree in newick format, \
                        with taxon names formatted species@sequence, \
                        where species matches species tree.")
    args = parser.parse_args()

    sTree = [x for x in tree_reader.read_tree_file_iter(args.speciesTree)][0]
    gTree = [x for x in tree_reader.read_tree_file_iter(args.geneTree)][0]

    species = sTree.lvsnms()
    genes = []
    for i in gTree.lvsnms():
        if "@" in i:
            i = i.split("@")
            genes.append((i[0], "@"+i[1]))
        else:
            genes.append((i, ""))

    mapping = {}
    for s in species:
        mapping[s] = []
        for g in genes:
            if s == g[0]:
                seq = g[0] + g[1]
                mapping[s].append(seq)

    for g in genes:
        if g[0] + g[1] not in mapping.values():
            sys.stderr.write(g[0] + g[1] + " is not represented by any \
                             taxon in the species tree. Please check.")
            sys.exit()

     

