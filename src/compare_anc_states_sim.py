#! /usr/bin/python3


"""
Script to summarise the results of ASR on simulated alignments. Takes a tree 
and a simulated alignment including ancestral nodes, as well as a reconstructed
alignment, and compares them for a specified node or set of nodes.
"""


import sys
import argparse
import tree_reader
from parse_fasta import parse_fasta


def count_diffs(seq1: str, seq2: str, ignoreGap=True) -> list:
    """iterate over two aligned sequences and return differences
    as state1posstate2, e.g. A235V"""
    if len(seq1) != len(seq2):
        sys.stderr.write("sequences not equal in length!")
        sys.exit()
    rawDiffs = [(j[0], i+1, j[1]) for i, j in
                enumerate(zip(seq1, seq2)) if j[0] != j[1]]
    if ignoreGap:
        diffs = [j for j in rawDiffs if "-" not in j]
    else:
        diffs = [j for j in rawDiffs]
    p = 1 - (len(diffs) / len(seq1))
    return diffs, p


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="Newick tree with internal nodes \
                        labelled")
    parser.add_argument("true_aln", help="FASTA-formatted known alignment, \
                        with ancestral sequences matching node labels in tree")
    parser.add_argument("recon_aln", help="FASTA-formatted reconstructed \
                        sequences, with names matching node labels in tree")
    parser.add_argument("nodes", help="One or more (space-separated) nodes to \
                        analyse", nargs="+")
    args = parser.parse_args()

    with open(args.tree, "r") as t:
        for s in t:
            s = s.strip()
            nwkString = s
    t = tree_reader.read_tree_string(nwkString)
    true = dict([x for x in parse_fasta(args.true_aln)])
    recon = dict([x for x in parse_fasta(args.recon_aln)])

    # first let's get raw p-distance between true and est

    with open("sim_err.tsv", "a") as outf:
        for n in args.nodes:
            true_seq = true[n]
            rec_seq = recon[n]
            _, p = count_diffs(true_seq, rec_seq)
            outf.write(n + "\t" + str(p) + "\n")

    # now we have to assess success of substitutions

    with open("sub_err.tsv", "a") as outf:
        for n in t.iternodes(order="preorder"):
            if n.label in args.nodes:
                for c in n.children:
                    if c.label in args.nodes:  # this is branch of interest
                        true_par = true[n.label]
                        true_ch = true[c.label]
                        recon_par = recon[n.label]
                        recon_ch = recon[c.label]
                        true_subs, _ = count_diffs(true_par, true_ch)
                        recon_subs, _ = count_diffs(recon_par, recon_ch)
                        cor = 0
                        cor_wrong_par = 0
                        cor_wrong_desc = 0
                        cor_wrong_par_wrong_desc = 0
                        fp = 0
                        fn = 0
                        print(true_subs)
                        print(recon_subs)
                        for sr in recon_subs:
                            sitematch = False
                            for s in true_subs:
                                if sr[1] == s[1]:  # sites match
                                    sitematch = True
                                    if sr[0] == s[0] and sr[1] == s[1]:
                                        # perfect match
                                        cor += 1
                                    elif sr[0] == s[0] and sr[1] != s[1]:
                                        # wrong desc
                                        cor_wrong_desc += 1
                                    elif sr[0] != s[0] and sr[1] == s[1]:
                                        # wrong par
                                        cor_wrong_par += 1
                                    elif sr[0] != s[0] and sr[1] != s[1]:
                                        # wrong par and desc
                                        cor_wrong_par_wrong_desc += 1
                            if not sitematch:
                                fp += 1  # false +ve
                        for s in true_subs:  # have to do converse pass to get
                            # false -ve
                            sitematch = False
                            for sr in recon_subs:
                                if s[1] == sr[1]:
                                    sitematch = True
                            if not sitematch:
                                fn += 1

                        outf.write("\t".join([n.label, c.label,
                                              str(len(true_subs)),
                                              str(len(recon_subs)),
                                              str(cor),
                                              str(cor_wrong_par),
                                              str(cor_wrong_desc),
                                              str(cor_wrong_par_wrong_desc),
                                              str(fp),
                                              str(fn)+"\n"]))
