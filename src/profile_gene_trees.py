#! /usr/bin/env python3


import sys
import argparse
import newick3
from phylo3 import Node


"""
Take two tree node objects: one an (unrooted) tree corresponding to a single branch with
two tips/clades on either side, and the other a gene tree. Returns true if each of the sets 
at the end of each of the branches is non-empty (see Lanfear & Hahn)

Each branch constraint is an 'unrooted' (tritomy root) tree with a minimum of one internal
branch, if each of the defining taxon subsets is a single taxon. If each instead consists of 
a clade, then there are at least 5 internal branches. A range of values is permitted for a mix
of definitional clades and singletons.

For example, consider we are test the split A B | C D, where A, B, C, and D represent four
clades or tips. In the case of tips, the input constraint tree is (1 internal branch):

((A,B),C,D));

    |---A
|---|
|   |---B
|---C
|---D

If instead A, B, C and D represent clades that can be decisive based on the presence of any one
of the taxa in the clades, the input constraint tree is (5 internal branches):

(((A1,A2),(B1,B2)),(C1,C2),(D1,D2));

        |---A1
    |---|
    |   |---A2
|---|    
|   |   |---B1
|   |---|
|       |---B2
|   |---C1 
|---|
|   |---C2
|   |---D1
|---|
    |---D2

In parsing this, we assume that any single internal branch tree represents the focal split,
while in trees with 1 > x >= 5 internal branches, the first branch with any non-tip descendants
is the focal split
"""

# def is_decisive(branch: Node, tree: Node) -> bool:
    
#     if branch.n_int_branch == 1:
#         if set(branch.lvsnms()) <= set(tree.lvsnms()):
#             return True
#     else:
#         g_taxa = set(tree.lvsnms())  # only do this once
#         clade_taxa = []
#         for c in branch.children:
#             if c.istip:
#                 clade_taxa.append([c.label])
#             elif all([n.istip for n in c.children]):
#                 clade_taxa.append(c.lvsnms())
#             else:  # this is the focal branch
#                 for n in c.children:
#                     if n.istip:
#                         clade_taxa.append([n.label])
#                     else:
#                         clade_taxa.append(n.lvsnms())
#         if all([set(labs) & g_taxa for labs in clade_taxa]):
#             return True
#     return False


def _get_partitions(branch: Node) -> list[list[str]]:
    """
    Extract the defining taxon partitions from a constraint tree.
    Each partition is a list of taxon names; decisiveness requires
    at least one representative from each.
    """
    partitions = []
    for c in branch.children:
        if c.istip:
            partitions.append([c.label])
        elif all(n.istip for n in c.children):
            partitions.append(c.lvsnms())
        else:  # focal branch: descend one level
            for n in c.children:
                if n.istip:
                    partitions.append([n.label])
                else:
                    partitions.append(n.lvsnms())
    return partitions


def is_decisive(branch: Node, tree: Node) -> bool:
    """
    True if the gene tree contains at least one taxon from every
    partition defined by the constraint branch (Lanfear & Hahn).
    """
    g_taxa = set(tree.lvsnms())

    if branch.n_int_branch == 1:
        return set(branch.lvsnms()) <= set(tree.lvsnms())
    else:
        partitions = _get_partitions(branch)

    return all(set(labs) & g_taxa for labs in partitions)


def profile_gene_trees():
    """main function"""
    parser = argparse.ArgumentParser(description="Profile gene trees for decisive branches")
    parser.add_argument("branch_tree", help="Newick tree with one branch of interest")
    parser.add_argument("gene_tree", help="Newick tree file with gene trees to profile")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    # Read the branch tree and gene trees
    branch_tree = newick3.parse_from_file(args.branch_tree)
    gene_tree = newick3.parse_from_file(args.gene_tree)

    print(f"{args.gene_tree}\t{is_decisive(branch_tree, gene_tree)}")
    

if __name__ == "__main__":
    profile_gene_trees()