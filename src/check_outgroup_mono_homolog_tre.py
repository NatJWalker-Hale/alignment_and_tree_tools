#! /usr/bin/env python3

import sys
import argparse
from ete3 import Tree, PhyloTree

"""Takes an outgroup rooted homolog tree and
checks if the outgroup sequences are monophyletic."""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-g", "--gene_tree",
                        help="Homolog tree to be assessed.",
                        required=True)
    parser.add_argument("-og", "--outgroupf",
                        help="Outgroup taxon names, one per line.",
                        required=True)

    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    args = parser.parse_args()

    og_list = []
    with open(args.outgroupf, "r") as ogf:
        for line in ogf:
            og_list.append(line.strip())

    tr = PhyloTree(args.gene_tree,
                   sp_naming_function=lambda node: node.name.split("@")[0])

    og_in_tr = []
    for l in tr.iter_leaves():
        if l.species in og_list:
            og_in_tr.append(l.species)

    print(args.gene_tree+"\t"+str(tr.check_monophyly(
                                                     values=og_in_tr,
                                                     target_attr="species")[0]
                                  )
          )
