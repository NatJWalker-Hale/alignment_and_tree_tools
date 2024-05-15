#! /usr/bin/python

"""
Script to write FastML-formatted ancestral state probabilities (prob.marginal.txt)
from raxml-ng-style .ancestralProbs. Requires raxml-ng reconstruction to be on a 
tree with the same node labelling scheme as FastML - use the script label_nodes_fastml.py
to produce this
"""

import sys
import argparse
import numpy as np
import pandas as pd
import newick3
import sequence as sq


AA = "YWVTSRQPNMLKIHGFEDCA"
AAORD = "ARNDCQEGHILKMFPSTWYV"


def tip_state_to_prob_str(state: str) -> str:
    out = ""
    if state == "-":
        out = " ".join(f"p({char})=1" for char in AA)
    else:
        out = f"p({state})=1"
    return out


# def node_row_to_prob_str(probs: pd.DataFrame) -> str:
#     out = ""
#     prob_dict = {j: probs.iloc[i] for i, j in enumerate(AAORD)}
#     for char in AA:
#         if prob_dict[char] != 0.:
#             out += f" p({char})={prob_dict[char]}"
#     return out


def node_row_to_prob_str(probs: pd.Series) -> str:
    out = ""
    char_probs = {char: probs.iloc[idx] for idx, char in enumerate(AAORD)}
    char_probs_ord = {char: char_probs[char] for char in AA}
    out = " ".join(f"p({char})={prob:.6f}" for char, prob in char_probs_ord.items() if prob != 0.)
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser("write FastML-formatted prob.marginal.txt from raxml-ng")
    parser.add_argument("probs", help="raxml-ng .ancestralProbs")
    parser.add_argument("alignment", help="input alignment for ASR")
    parser.add_argument("tree", help="input tree for ASR")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    curroot = newick3.parse_from_file(args.tree)

    aln = dict(sq.parse_fasta(args.alignment))
    cols = sq.get_columns(aln)

    probs = pd.read_csv(args.probs, sep="\t", header=0)
    # out = {col: {} for col in cols}
    # for n in curroot.iternodes(order=1):
    #     if n.istip:
    #         for col in out:
    #             out[col][n.label] = tip_state_to_prob_str(cols[col][n.label])
    #     else:
    #         for col in out:
    #             probs_series = probs[(probs['Node'] == n.label) &
    #                                  (probs['Site'] == col + 1)].iloc[:, 3:].squeeze()
    #             out[col][n.label] = node_row_to_prob_str(probs_series)

    with open("prob.marginal.txt", "w", encoding="utf-8") as outf:
        for col in cols:
            if col+1 == 1:
                outf.write(f"marginal probabilities at position: {col+1}\n")
            else:
                outf.write(f"\nmarginal probabilities at position: {col+1}\n")
            for n in curroot.iternodes(order=1):
                if n.istip:
                    outf.write(f"of node: {n.label}: {tip_state_to_prob_str(cols[col][n.label])}\n")
                else:
                    probs_series = probs[(probs['Node'] == n.label) &
                                         (probs['Site'] == col + 1)].iloc[:, 3:].squeeze()
                    outf.write(f"of node: {n.label}: {node_row_to_prob_str(probs_series)}\n")
        outf.write("\n\n++++++++++++++++++++++++ marginal probs +++++++++++++++++++++++++++++++\n\n")
        outf.write("node,site,A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V\n")
        outf.write(probs.iloc[:, np.r_[0, 1, 3:22]].to_csv(index=False, header=False))
        
    # out_str = ""
    # for col in out:
    #     if col+1 == 1:
    #         out_str += f"marginal probabilities at position: {col+1}\n"
    #     else:
    #         out_str += f"\nmarginal probabilities at position: {col+1}\n"
    #     for node, prob_str in out[col].items():
    #         out_str += f"of node: {node}: {prob_str}\n"

    # out_str += "\n\n++++++++++++++++++++++++ marginal probs +++++++++++++++++++++++++++++++\n\n"

    # out_str += "node,site,A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V\n"

    # out_str += probs.iloc[:, np.r_[0, 1, 3:22]].to_csv(index=False, header=False)

    # print(out_str)



