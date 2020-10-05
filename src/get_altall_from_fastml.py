import sys
import argparse
import pandas as pd


def parse_fastml_prob(ppath, gpath, ipath):
    p = pd.read_csv(ppath)
    g = pd.read_csv(gpath, sep="\t")
    ip = pd.read_csv(ipath, sep="\t")
    return p, g, ip


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--probfile", help="prob.marginal.csv from FastML \
                    output.")
    ap.add_argument("-g", "--gapfile", help="Ancestral_MaxMarginalProb_Char_\
                    Indel.txt from FastML output.")
    ap.add_argument("-i", "--indelprob", help="IndelsMarginalProb.txt from \
                    FastML output.")
    ap.add_argument("-n", "--node", help="Node to calculate AltAll, e.g. N14.")
    ap.add_argument("-gp", "--gapmode", help="Mode to include gaps. Options \
                    are [prob/fixed]. prob includes the alternate gap \
                    character if the PP is >= 0.2, while fixed treats \
                    the gap reconstruction as fixed, i.e. preferring the \
                    gap state if PP >= 0.5")
    args = ap.parse_args()

    outdict = {args.node: ""}

    pprobs, gprobs, iprobs = parse_fastml_prob(args.probfile, args.gapfile,
                                               args.indelprob)

    """FastML stores the necessary reconstruction data in three files:
    prob.marginal.csv, which contains the probabilities for the marginal
    sequence reconstruction; Ancestral_MaxMarginalProb_Char_Indel.txt,
    which contains the max posterior state and its probability, at each
    site; and IndelsMarginalProb.txt, which contains the probability of the
    indel state, referenced by sites in the original alignment. The coding
    of this latter file is such that it doesn't uniformly contain the
    probability of insertion or deletion, but contains the prob of the
    alternate state to that called in
    Ancestral_MaxMarginalProb_Char_Indel.txt. So if the character state in
    Ancestral_MaxMarginalProb_Char_Indel.txt is sequence, the probability
    for that site in IndelsMarginalProb.txt is the probability of gap."""

    iprob_dict = {}

    ntable = iprobs.loc[iprobs["Node"] == args.node]
    cval = max(ntable["Pos"].values.tolist())
    for i in [j for j in range(1, cval+1)]:
        if i in ntable.Pos.values:
            iprob_dict[i] = ntable.loc[ntable["Pos"] ==
                                       i]["Prob_Of_Indel"].values[0]

    for index, grow in gprobs.loc[gprobs["Node"] == args.node].iterrows():
        if args.gapmode == "prob":
            if grow["Char"] == "-":  # gap
                if grow["CharProb"] < 0.8:  # uncertain gap
                    # get sequence state and set state to max prob sequence
                    prow = pprobs.loc[(pprobs["Ancestral Node"] ==
                                      grow["Node"]) & (pprobs["Pos"] ==
                                      grow["Pos_on_MSA"])].iloc[:, 2:]
                    sprobs = prow.values.tolist()[0]
                    m = max(sprobs)
                    idx = [i for i, j in enumerate(sprobs) if j == m][0]
                    outdict[args.node] += [*prow][idx]
                else:  # certain gap
                    outdict[args.node] += "-"  # set state to gap
            else:  # not gap
                if int(grow["Pos_on_MSA"]) in iprob_dict.keys():
                    if iprob_dict[int(grow["Pos_on_MSA"])] == 1:  # certain
                        istate = 0
                    else:
                        istate = iprob_dict[int(grow["Pos_on_MSA"])]
                if istate > 0.2:  # certain (like for sequence)
                    outdict[args.node] += "-"
                else:  # uncertain
                    prow = pprobs.loc[(pprobs["Ancestral Node"] ==
                                      grow["Node"]) & (pprobs["Pos"] ==
                                      grow["Pos_on_MSA"])].iloc[:, 2:]
                    sprobs = prow.values.tolist()[0]
                    vals = sorted([(i, j) for i, j in enumerate(sprobs)],
                                  key=lambda x: x[1], reverse=True)
                    idx = vals[1][0] if vals[1][1] > 0.2 else vals[0][0]
                    outdict[args.node] += [*prow][idx]
        elif args.gapmode == "fixed":
            if grow["Char"] == "-":
                outdict[args.node] += "-"
            else:
                prow = pprobs.loc[(pprobs["Ancestral Node"] ==
                                  grow["Node"]) & (pprobs["Pos"] ==
                                  grow["Pos_on_MSA"])].iloc[:, 2:]
                sprobs = prow.values.tolist()[0]
                vals = sorted([(i, j) for i, j in enumerate(sprobs)],
                              key=lambda x: x[1], reverse=True)
                idx = vals[1][0] if vals[1][1] > 0.2 else vals[0][0]
                outdict[args.node] += [*prow][idx]

for k, v in outdict.items():
    print(">"+k)
    print(v)
