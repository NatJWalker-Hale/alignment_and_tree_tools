#! /usr/bin/env python3


import sys
import argparse
import newick3
from treenode import Node


def is_decisive(branch: Node, tree: Node) -> bool:
    """Take two tree node objects: one an (unrooted) tree corresponding to a single branch with
    two tips/clades on either side, and the other a gene tree. Returns true if each of the sets 
    at the end of each of the branches is non-empty (see Lanfear and Hahn)"""
    b0, b1 = branch.children  # Assuming the branch tree has one internal node with two branches
    s0, s1 = (set(l.label for l in b0.children[0].leaves()),
              set(l.label for l in b0.children[1].leaves()))
    s2, s3 = (set(l.label for l in b1.children[0].leaves()),
              set(l.label for l in b1.children[1].leaves()))
    gene_leaves = set(l.label for l in tree.leaves())
    print(s0)
    print(gene_leaves)
    print(s0.intersection(gene_leaves))


def profile_gene_trees():
    """main function"""
    parser = argparse.ArgumentParser(description="Profile gene trees for decisive branches")
    parser.add_argument("branch_tree", help="Newick tree with one branch of interest")
    parser.add_argument("gene_tree", help="Newick tree file with gene trees to profile")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    # Read the branch tree and gene trees
    branch_tree = newick3.parse_from_file(args.branch_tree)
    gene_tree = newick3.parse_from_file(args.gene_tree)

    is_decisive(branch_tree, gene_tree)
    

if __name__ == "__main__":
    profile_gene_trees()