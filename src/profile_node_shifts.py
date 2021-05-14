#! /usr/bin/python3

import sys
import argparse
import tree_reader


def parse_og(path):
    with open(path, "r") as ogfile:
        og = [line.strip() for line in ogfile.readlines()]
    return og


def calc_shift_dist(params1, params2):
    """Calculates euclidean distance between two compositional vectors"""
    dist = sum([(a - b)**2 for a, b in zip(params1, params2)])**0.5
    return dist


def check_shift(parent, child, format="NEX"):
    """checks if parent and child node have different models"""
    if format == "NEX":
        bm = parent.note.split(",")[0].split("=")[1]
        cm = child.note.split(",")[0].split("=")[1]
    elif format == "NWK":
        bm = parent.label
        cm = child.label
    else:
        sys.stderr.write("tree format unrecognised\n")
        sys.exit()
    if bm != cm:
        return True
    else:
        return False


def get_model(node, format="NEX"):
    if format == "NEX":
        m = node.note.split(",")[0].split("=")[1]
    elif format == "NWK":
        m = node.label
    else:
        sys.stderr.write("tree format unrecognised\n")
        sys.exit()
    return m


def get_params(node):
    """extract model params"""
    pstring = node.note.split("{")[1].split("}")[0]
    p = [float(x) for x in pstring.split(",")]
    return p


def number_tree(root):
    """iterates over nodes and numbers them.
    modifies in place"""
    c = 0
    for n in root.iternodes(order="preorder"):
        n.number = c
        c += 1


# def count_shifts(inroot, min_clade_size, min_overlap, format="NEX"):
#     countDict = {
#         "Testable nodes": 0,
#         "Testable orthologs": 0,
#         "Orthologous shifts": 0,
#         "Testable paralogs": 0,
#         "Shifts on duplication": 0,
#         "Paralogous shifts": 0,
#         "Paralogous contrasts": 0
#     }

#     shiftDict = {}
#     count = 0
#     for n in inroot.iternodes(order="preorder"):
#         if count == 0:
#             count += 1
#             continue
#         if n.istip:
#             count += 1
#             continue
#         if len(n.leaves()) >= min_clade_size:
#             # only consider clades large enough to test
#             countDict["Testable nodes"] += 1
#             ntaxa = len(set([x.split("@")[0] for x in n.lvsnms()]))
#             ntips = len(n.leaves())
#             if ntaxa == ntips:  # 1-to-1 ortholog
#                 countDict["Testable orthologs"] += 1
#                 # first, compare for shift to parent
#                 # can happen if clade is child of 1-to-many
#                 # or many-to-many
#                 if check_shift(n.parent, n, format):
#                     if (count, n.parent, n) in shiftDict.keys():
#                         continue
#                     else:
#                         shiftDict[(count, n.parent, n)] = "OS"
#                 for c in n.children:
#                     if len(c.leaves()) >= min_clade_size:
#                         if check_shift(n, c, format):
#                             shiftDict[(count, n, c)] = "OS"  # key is num-parent-child tuple, value type
#             else:  # contains duplicated taxa
#                 child0, child1 = n.children
#                 names0 = set([x.split("@")[0] for x in child0.lvsnms()])
#                 names1 = set([x.split("@")[0] for x in child1.lvsnms()])
#                 if len(names0.intersection(names1)) >= min_overlap:  # n is duplication node
#                     if check_shift(n.parent, n, format):  # shift on dup node
#                         shiftDict[(count, n.parent, n)] = "SD"
#                     dupShift = 0
#                     if len(child0.leaves()) >= min_clade_size:
#                         countDict["Testable paralogs"] += 1
#                         if check_shift(n, child0, format):  # first paralog has shift
#                             shiftDict[(n, child0)] = "PS"
#                             dupShift += 1
#                         else:  # sub-loop to check for nested shift
#                             for cn in child0.iternodes(order="preorder"):
#                                 for ccn in cn.children:
#                                     if len(ccn.leaves()) >= min_clade_size:
#                                         if check_shift(cn, ccn, format):
#                                             shiftDict[(cn, ccn)] = "PS"
#                                             dupShift += 1
#                     if len(child1.leaves()) >= min_clade_size:
#                         countDict["Testable paralogs"] += 1
#                         if check_shift(n, child1, format):  # second paralog has shift
#                             shiftDict[(n, child1)] = "PS"
#                             dupShift += 1
#                         else:  # sub-loop to check for nested shift
#                             for cn in child1.iternodes(order="preorder"):
#                                 for ccn in cn.children:
#                                     if len(ccn.leaves()) >= min_clade_size:
#                                         if check_shift(cn, ccn, format):
#                                             shiftDict[(cn, ccn)] = "PS"
#                                             dupShift += 1
#                     # if dupShift == 1:  # can be 0, 1 or 2
#                     #     countDict["Paralogous contrasts"] += 1
#                 else:  # no overlap in children taxa - still "orthologous", 1-to-many
#                     countDict["Testable orthologs"] += 1
#                     try:
#                         if len(n.leaves()) >= min_clade_size:
#                             if check_shift(n.parent, n, format):
#                                 if (n.parent, n) in shiftDict.keys():
#                                     continue
#                                 else:
#                                     shiftDict[(n.parent, n)] = "OS"
#                     except AttributeError:
#                         sys.stderr.write("no parent node, check\n")
#             count += 1
#     for _, v in shiftDict.items():
#         if v == "OS":
#             countDict["Orthologous shifts"] += 1
#         elif v == "SD":
#             countDict["Shifts on duplication"] += 1
#         elif v == "PS":
#             countDict["Paralogous shifts"] += 1
#     return countDict, shiftDict

def count_shifts(inroot, min_clade_size, min_overlap, format="NEX"):
    # different version which compares back to parent instead of forward to
    # child
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

    nodeDict = {
        "node": 0,
        "shift": False,
        "model": 0
    }

    count = 0
    for n in inroot.iternodes(order="preorder"):
        if count == 0:  # root
            count += 1
            continue
        if n.istip:
            count += 1
            continue
        ntips = len(n.leaves())
        if ntips >= min_clade_size:
            # only consider clades large enough
            # to be tested
            countDict["Testable nodes"] += 1
            ntaxa = len(set([x.split("@")[0] for x in n.lvsnms()]))
            if ntaxa == ntips:  # 1-to-1 orthology
                countDict["Testable nodes"] += 1
                if check_shift(n.parent, n, format):
                    if (n.parent, n) in shiftDict.keys():
                        # already registered for whatever reason
                        # in subloop
                        continue
                    else:
                        shiftDict[(n.parent, n)] = "OS"
                        m = get_model(n, format)
                        nodeDict["node": n.number,
                                 "shift": True,
                                 "model": m]
            else:  # contains duplicated taxa
                child0, child1 = n.children
                names0 = set([x.split("@")[0] for x in child0.lvsnms()])
                names1 = set([x.split("@")[0] for x in child1.lvsnms()])
                if len(names0.intersection(names1)) >= min_overlap:
                    # n is duplication node of acceptable size
                    if check_shift(n.parent, n, format):
                        # shift is on duplication node
                        if (n.parent, n) in shiftDict.keys():
                            continue
                        else:
                            shiftDict[(n.parent, n)] = "SD"
                            m = get_model(n, format)
                            nodeDict["node": n.number,
                                     "shift": True,
                                     "model": m]
                    if len(child0.leaves()) >= min_clade_size:
                        # paralog that can be tested
                        countDict["Testable paralogs"] += 1
                        if check_shift(n, child0, format):
                            # shift in paralog immediately
                            # below dup node
                            if (n, child0) in shiftDict.keys():
                                continue
                            else:
                                shiftDict[(n, child0)] = "PS"
                                m = get_model(child0)
                                nodeDict["node": child0.number,
                                         "shift": True,
                                         "model": m]
                        else:  # start the sub-loop
                            for c in child0.iternodes(order="preorder"):
                                if len(c.leaves()) >= min_clade_size:
                                    if check_shift(c.parent, c, format):
                                        if (c.parent, c) in shiftDict.keys():
                                            continue
                                        else:
                                            shiftDict[(c.parent, c)] = "PS"
                                            m = get_model(c)
                                            nodeDict["node": c.number,
                                                    "shift": True,
                                                    "model": m]
                    if len(child1.leaves()) >= min_clade_size:
                        # paralog that can be tested
                        countDict["Testable paralogs"] += 1
                        if check_shift(n, child0, format):
                            # shift in paralog immediately
                            # below dup node
                            if (n, child0) in shiftDict.keys():
                                continue
                            else:
                                shiftDict[(n, child0)] = "PS"
                                m = get_model(child0)
                                nodeDict["node": child0.number,
                                         "shift": True,
                                         "model": m]
                        else:  # start the sub-loop
                            for c in child1.iternodes(order="preorder"):
                                if len(c.leaves()) >= min_clade_size:
                                    if check_shift(c.parent, c, format):
                                        if (c.parent, c) in shiftDict.keys():
                                            continue
                                        else:
                                            shiftDict[(c.parent, c)] = "PS"
                                            m = get_model(c)
                                            nodeDict["node": c.number,
                                                     "shift": True,
                                                     "model": m]
                elif len(names0.intersection(names1)) == 0:
                    # not overlapping - this is node which has paralogs below
                    countDict["Testable orthologs"] += 1
                    if check_shift(n.parent, n, format):
                        if (n.parent, n) in shiftDict.keys():




            



if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--trimfile", help="file of taxa \
                        to remove sequences from prior to \
                        counting, one per line")
    parser.add_argument("tree", help="janus output NEXUS tree, \
                        *.gophy.results.tre")
    args = parser.parse_args()
    tFile = args.tree

    with open(tFile, "r") as t:
        for s in t.readlines():
            if s.startswith("tree"):
                nwkString = s.split("=", 1)[1].lstrip().rstrip()
                treForm = "NEX"
            elif s.startswith("("):
                nwkString = s.strip()
                treForm = "NWK"
                
    curroot = tree_reader.read_tree_string(nwkString)
    # for n in curroot.iternodes(order="preorder"):
    #     print(n.label)
    if args.trimfile is not None:
        OGs = parse_og(args.trimfile)
        OGseqs = []
        for i in curroot.lvsnms():
            if i.split("@")[0] in OGs:
                OGseqs.append(i)
        print(OGseqs)
        for n in curroot.iternodes(order="preorder"):
            if set(n.lvsnms()) == set(OGseqs):  # assumes OG monophyletic
                n.prune()

    cDict, shifts = count_shifts(curroot, 10, 10, treForm)
    for k, v in cDict.items():
