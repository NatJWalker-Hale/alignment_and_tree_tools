#! /usr/bin/python3

import sys
import argparse
import pandas as pd
import numpy as np
from numpy.random import choice


STATEDICT = {0: "A",
             1: "C",
             2: "D",
             3: "E",
             4: "F",
             5: "G",
             6: "H",
             7: "I",
             8: "K",
             9: "L",
             10: "M",
             11: "N",
             12: "P",
             13: "Q",
             14: "R",
             15: "S",
             16: "T",
             17: "V",
             18: "W",
             19: "Y"
             }


def sample_states(v):
    """with vector of probs for one pos, order
    A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y"""
    p = v / sum(v)  # since values don't sum
    draw = choice(a=list(range(20)),
                  size=1,
                  p=p)
    # since probs don't always sum
    state = STATEDICT[int(draw)]
    return state


# come back to how to do this
# def sample_states_vec(v):
#     draw = choice(a=list(range(20)),
#                   size=1,
#                   p=[round(x, 4) for x in v])
#     return draw


def sample_node(probs, node, n=1):
    nodeProbs = probs.loc[probs['Ancestral Node'] == node]
    s = [sample_states(row) for row in nodeProbs.loc[:, "A":"Y"].to_numpy()]
    return "".join(s)


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("probs", help="marginal probabilities \
                        comma separated")
    parser.add_argument("node", help="node to sample")
    args = parser.parse_args()

    probs = pd.read_csv(args.probs)
    seq = sample_node(probs, args.node)
    print(">"+args.node)
    print(seq)

    