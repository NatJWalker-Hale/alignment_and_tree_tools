#! /usr/bin/env python3


import sys
import argparse
import newick3
import phylo3
from phylo3 import Node


"""
Written by NWH with assistance from Claude Opus 4.6

Take two tree node objects: one an (unrooted) tree corresponding to a single branch with
two tips/clades on either side, and the other a gene tree. Returns true if each of the sets
at the end of each of the branches is non-empty (see Lanfear & Hahn).

Each branch constraint is an 'unrooted' (tritomy root) tree with a minimum of one internal
branch (not counting root), if each of the defining taxon subsets is a single taxon. If each
instead consists of a clade, then there are at least 5 internal branches. A range of values is
permitted for a mix of definitional clades and singletons.

For example, consider we are test the split A B | C D, where A, B, C, and D represent four
clades or tips. In the case of tips, the input constraint tree is (1 internal branch):

((A,B),C,D));

    |---A
|---|
|   |---B
|
|---C
|
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
|
|   |---C1
|---|
|   |---C2
|
|   |---D1
|---|
    |---D2

In parsing this, we assume that any single internal branch tree represents the focal split,
while in trees with 1 > x >= 5 internal branches, the first branch with any non-tip descendants
is the focal split.

If a given gene tree is decisive for the constraint, then all we need to do is check for the
existence of the split defined by the focal branch - we do not require e.g. that the clades are
monophyletic in the gene trees. That is to say, the decisiveness check determines if they could
theoretically reproduce the exact relationships in the constraint, but we do not require it to.

There is also a quartet mode to get around the problem of incidental tips causing many gene trees
to be entirely non-concordant.
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
        return set(branch.lvsnms()) <= g_taxa
    partitions = _get_partitions(branch)
    return all(not set(labs).isdisjoint(g_taxa) for labs in partitions)


def get_split(branch: Node) -> list[set[str]]:
    """
    Extract the defining taxon partitions from a constraint tree.
    Use the focal branch to define two sets of taxa in a bipartition.
    Return the bipartition as a list of sets of taxa.
    """
    for c in branch.children:
        if not c.istip and not all(n.istip for n in c.children):
            s1 = set(c.lvsnms())
            s2 = set(branch.lvsnms()) - s1
            return [s1, s2]
    raise ValueError("No focal branch found in constraint tree")


def compare_split(bp1: list[set], bp2: list[set]) -> bool:
    """True if concordant (at least one empty intersection)."""
    return (bp1[0].isdisjoint(bp2[0]) or
            bp1[0].isdisjoint(bp2[1]) or
            bp1[1].isdisjoint(bp2[0]) or
            bp1[1].isdisjoint(bp2[1]))


def conflicts(bp2_side: set, bp1_side: dict) -> bool:
    seen = set()
    for taxon in bp2_side:
        side = bp1_side.get(taxon)
        if side is not None:
            seen.add(side)
            if len(seen) == 2:
                return True
    return False


# def quartet_concordance(partitions: list[list[str]], gene_taxa: set, gene_splits: list[set]):
#     parts = [set(p) & gene_taxa for p in partitions]
#     # print(parts)

#     # skip if not decisive (any partition unrepresented)
#     if any(len(p) == 0 for p in parts):
#         return None

#     A, B, C, D = parts
#     focal = 0
#     alt1 = 0
#     alt2 = 0

#     # print(gene_splits)

#     for below, above in gene_splits:
#         a0, a1 = len(A & below), len(A & above)
#         b0, b1 = len(B & below), len(B & above)
#         c0, c1 = len(C & below), len(C & above)
#         d0, d1 = len(D & below), len(D & above)

#         focal += a0 * b0 * c1 * d1 + a1 * b1 * c0 * d0
#         alt1  += a0 * c0 * b1 * d1 + a1 * c1 * b0 * d0
#         alt2  += a0 * d0 * b1 * c1 + a1 * d1 * b0 * c0

#     # total = len(A) * len(B) * len(C) * len(D)
#     # total = sum([focal, alt1, alt2])   # this doesn't work because it undercounts that every
#     # internal branch below the split contributes the same quartet

#     return focal / total, alt1 / total, alt2 / total


def quartet_concordance(partitions: list[list[str]],
                        gene_tree: Node,
                        ) -> tuple[float, float, float, float] | None:
    """
    Quartet concordance for a focal branch defined by four taxon partitions
    A, B, C, D. Enumerates every (a, b, c, d) one from each partition and
    determines its unrooted topology in the gene tree by LCA depth: the
    pair with the deepest pairwise LCA is sister, and the remaining two
    are the partner pair. Returns (focal, alt1, alt2, unresolved) as
    proportions of |A|*|B|*|C|*|D|, or None if any partition is empty
    after restricting to the gene tree's taxa.
    """
    # one DFS pass: node depth (root = 0) + lookup of leaves by label
    depth_of = {gene_tree: 0}
    leaf_by_name = {}
    stack = []
    for n in gene_tree.descendants():
        if n.istip:
            leaf_by_name[n.label] = n
        depth_of[n] = depth_of[n.parent] + 1

    # while stack:
    #     node = stack.pop()
    #     if node.istip:
    #         leaf_by_name[node.label] = node
    #     for child in node.children:
    #         depth_of[child] = depth_of[node] + 1
    #         stack.append(child)

    parts = [[leaf_by_name[t] for t in p if t in leaf_by_name]
            for p in partitions]
    if any(not p for p in parts):
        return None

    # ancestor set per quartet leaf - LCA = deepest node in the intersection
    anc = {}
    for leaf in {n for p in parts for n in p}:
        s, n = set(), leaf
        while n is not None:
            s.add(n)
            n = n.parent
        anc[leaf] = s

    A, B, C, D = parts
    focal = alt1 = alt2 = unresolved = 0
    for a in A:
        for b in B:
            d_ab = max(depth_of[n] for n in anc[a] & anc[b])
            for c in C:
                d_ac = max(depth_of[n] for n in anc[a] & anc[c])
                d_bc = max(depth_of[n] for n in anc[b] & anc[c])
                for d in D:
                    d_ad = max(depth_of[n] for n in anc[a] & anc[d])
                    d_bd = max(depth_of[n] for n in anc[b] & anc[d])
                    d_cd = max(depth_of[n] for n in anc[c] & anc[d])
                    pair_depth = {"ab": d_ab, "cd": d_cd,
                                "ac": d_ac, "bd": d_bd,
                                "ad": d_ad, "bc": d_bc}
                    top = max(pair_depth.values())
                    sisters = {p for p, dp in pair_depth.items() if dp == top}
                    if sisters <= {"ab", "cd"}:
                        focal += 1
                    elif sisters <= {"ac", "bd"}:
                        alt1 += 1
                    elif sisters <= {"ad", "bc"}:
                        alt2 += 1
                    else:
                        unresolved += 1

    total = len(A) * len(B) * len(C) * len(D)
    return focal / total, alt1 / total, alt2 / total, unresolved / total


def profile_gene_trees():
    """main function"""
    parser = argparse.ArgumentParser(description="Profile gene trees for decisive branches")
    parser.add_argument("branch_trees", help="Newick trees with one branch of interest")
    parser.add_argument("gene_trees", nargs="+", help="input gene tree file(s)")
    parser.add_argument("-s", "--strict", action="store_true",
                        help="Require input trees to be decisive for all constraints")
    parser.add_argument("-q", "--quartets", action="store_true",
                        help="Use exhaustive quartet sampling instead of bipartitions per gene \
                        tree")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    branch_trees = [t for t in newick3.read_tree_file_iter(args.branch_trees)]

    # precompute splits and lookups once
    constraints = {}
    branch = 0
    for branch_tree in branch_trees:
        nni1, nni2 = branch_tree.nni()
        sys.stderr.write(f"branch\n")
        sys.stderr.write(f"{newick3.to_string(branch_tree, length_fmt="")};\n")
        sys.stderr.write(f"branch {branch}, NNI 1\n")
        sys.stderr.write(f"{newick3.to_string(nni1, length_fmt="")};\n")
        sys.stderr.write(f"branch {branch}, NNI 2\n")
        sys.stderr.write(f"{newick3.to_string(nni2, length_fmt="")};\n")
        # print(newick3.to_string(branch_tree))
        for b in [branch_tree, nni1, nni2]:
            bp1 = get_split(b)
            bp1_side = {taxon: 0 for taxon in bp1[0]}
            bp1_side.update({taxon: 1 for taxon in bp1[1]})
            try:
                constraints[branch].append((b, bp1_side))
            except KeyError:
                constraints[branch] = [(b, bp1_side)]
        branch += 1
    # print(constraints)

    for tree_file in args.gene_trees:
        gene_tree = newick3.parse_from_file(tree_file)
        if gene_tree.is_rooted():
            gene_tree.unroot()
        gene_taxa = set(gene_tree.lvsnms())
        gene_splits = list(phylo3.get_gene_tree_splits(gene_tree))
        # print(gene_splits)

        if args.strict:
            if not all(is_decisive(entries[0][0], gene_tree)
                       for entries in constraints.values()):
                for branch, entries in constraints.items():
                    for con in range(len(entries)):
                        print(f"{tree_file}\t{branch}\t{con}\tundecisive")
                continue

        if args.quartets:
            for branch, entries in constraints.items():
                branch_tree, _ = entries[0]
                result = quartet_concordance(_get_partitions(branch_tree), gene_tree)
                if result is None:
                    print(f"{tree_file}\t{branch}\tundecisive")
                else:
                    qf, q1, q2, qu = result
                    print(f"{tree_file}\t{branch}\t"
                            f"{qf:.4f}\t{q1:.4f}\t{q2:.4f}\t{qu:.4f}")
            continue

        for branch, constraint in constraints.items():
            for con, (branch_tree, bp1_side) in enumerate(constraint):
                if not args.strict and not is_decisive(branch_tree, gene_tree):
                    print(f"{tree_file}\t{branch}\t{con}\tundecisive")
                    continue
                concordant = True
                for bp2 in gene_splits:
                    if conflicts(bp2[0], bp1_side) and conflicts(bp2[1], bp1_side):
                        concordant = False
                        break
                print(f"{tree_file}\t{branch}\t{con}\t{concordant}")

if __name__ == "__main__":
    profile_gene_trees()
    # print(compare_split([{'b', 'c', 'a'}, {'e', 'd'}],
    #                     [{'b', 'a'}, {'e', 'c', 'd'}]))
    # True
    # print(compare_split([{'b', 'c', 'a'}, {'e', 'd'}],
    #                     [{'b', 'e'}, {'a', 'c', 'd'}]))
    # False
    # print(compare_split([{'b', 'c', 'a'}, {'e', 'd'}],
    #                     [{'a', 'b'}, {'c', 'd'}]))
    # True