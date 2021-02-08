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
        "Paralogous shifts": 0
    }

    paramDict = {
        "Orthologous shifts": [],
        "Paralogous shifts": []
    }

    for n in inroot.iternodes(order="preorder"):
        nparam = get_params(n)
        if n.istip():
            continue
        if len(n.leaves) >= min_clade_size:
            # only consider clades large enough to test
            countDict["Testable nodes"] += 1
            ntaxa = len(set(n.lvsnms))
            ntips = len(n.leaves)
            if ntaxa == ntips:  # ortholog
                countDict["Testable orthologs"] += 1
                for i in n.children:
                    if check_shift(n, i):
                        countDict["Orthologous shifts"] += 1
                        cparam = get_params(c)
                        diff = calc_shift_dist(nparam, cparam)
                        paramDict["Orthologous shifts"].append((
                            nparam,
                            cparam,
                            diff
                        ))
            else:  # contains duplicated taxa
                child0, child1 = n.children
                names0 = set([x.split("@")[0] for x in child0.lvsnms])
                names1 = set([x.split("@")[0] for x in child1.lvsnms])
                if len(names0.intersection(names1)) >= min_overlap:
                    





                        
            


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="janus output NEXUS tree, \
                        *.gophy.results.tre")
    args = parser.parse_args()
    tFile = args.tree

    with open(tFile, "r") as t:
        for s in t.readlines():
            if s.startswith("tree"):
                nwkString = s.split("=", 1)[1].lstrip().rstrip()

    curroot = tree_reader.read_tree_string(nwkString)
    print(curroot.note)