#! /usr/bin/python3

import sys
import argparse
import tree_reader


def calc_shift_dist(params1, params2):
    """Calculates euclidean distance between two compositional vectors"""
    dist = sum([(a - b)**2 for a, b in zip(params1, params2)])**0.5
    return dist


def check_shift(parent, child):
    """checks if parent and child node have different models"""
    bm = parent.note.split(",")[0].split("=")[1]
    cm = child.note.split(",")[0].split("=")[1]
    if bm != cm:
        return True
    else:
        return False


def get_params(node):
    """extract model params"""
    pstring = node.note.split("{")[1].split("}")[0]
    p = [float(x) for x in pstring.split(",")]
    return p


def count_shifts(inroot, min_clade_size, min_overlap):
    countDict = {
        "Testable nodes": 0,
        "Testable orthologs": 0,
        "Orthologous shifts": 0,
        "Testable paralogs": 0,
        "Shifts on duplication": 0,
        "Paralogous shifts": 0,
        "Paralogous contrasts": 0
    }

    shiftDict = {}
    for n in inroot.iternodes(order="preorder"):
        if n.istip:
            continue
        if len(n.leaves()) >= min_clade_size:
            # only consider clades large enough to test
            countDict["Testable nodes"] += 1
            ntaxa = len(set([x.split("@")[0] for x in n.lvsnms()]))
            ntips = len(n.leaves())
            if ntaxa == ntips:  # 1-to-1 ortholog
                countDict["Testable orthologs"] += 1
                for c in n.children:
                    if check_shift(n, c):
                        shiftDict[(n, c)] = "OS"  # key is parent-child tuple, value type
            else:  # contains duplicated taxa
                child0, child1 = n.children
                names0 = set([x.split("@")[0] for x in child0.lvsnms()])
                names1 = set([x.split("@")[0] for x in child1.lvsnms()])
                if len(names0.intersection(names1)) >= min_overlap:  # n is duplication node
                    if check_shift(n.parent, n):  # shift on dup node
                        shiftDict[(n.parent, n)] = "SD"
                    dupShift = 0
                    if len(child0.leaves()) >= min_clade_size:
                        countDict["Testable paralogs"] += 1
                        if check_shift(n, child0):  # first paralog has shift
                            shiftDict[(n, child0)] = "PS"
                            dupShift += 1
                        else:  # sub-loop to check for nested shift
                            for cn in child0.iternodes(order="preorder"):
                                for ccn in cn.children:
                                    if check_shift(cn, ccn):
                                        shiftDict[(cn, ccn)] = "PS"
                                        dupShift += 1
                    if len(child1.leaves()) >= min_clade_size:
                        countDict["Testable paralogs"] += 1
                        if check_shift(n, child1):  # second paralog has shift
                            shiftDict[(n, child1)] = "PS"
                            dupShift += 1
                        else:  # sub-loop to check for nested shift
                            for cn in child1.iternodes(order="preorder"):
                                for ccn in cn.children:
                                    if check_shift(cn, ccn):
                                        shiftDict[(cn, ccn)] = "PS"
                                        dupShift += 1
                    # if dupShift == 1:  # can be 0, 1 or 2
                    #     countDict["Paralogous contrasts"] += 1
                else:  # no overlap in children taxa - still "orthologous", 1-to-many
                    countDict["Testable orthologs"] += 1
                    try:
                        if check_shift(n.parent, n):
                            if (n.parent, n) in shiftDict.keys():
                                continue
                            else:
                                shiftDict[(n.parent, n)] = "OS"
                    except AttributeError:
                        sys.stderr.write("no parent node, check\n")
    for _, v in shiftDict.items():
        if v == "OS":
            countDict["Orthologous shifts"] += 1
        elif v == "SD":
            countDict["Shifts on duplication"] += 1
        elif v == "PS":
            countDict["Paralogous shifts"] += 1
    return countDict, shiftDict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--trim", help="remove outgroup \
                        prior to profiling", type=bool, default=False)
    parser.add_argument("tree", help="janus output NEXUS tree, \
                        *.gophy.results.tre")
    args = parser.parse_args()
    tFile = args.tree

    with open(tFile, "r") as t:
        for s in t.readlines():
            if s.startswith("tree"):
                nwkString = s.split("=", 1)[1].lstrip().rstrip()

    curroot = tree_reader.read_tree_string(nwkString)
    if args.trim:
        curroot = curroot.children[1]
    cDict, shifts = count_shifts(curroot, 10, 10)
    print(tFile, cDict)