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


def populate_node_dict(root):
    """populates an empty nodeDict from tree.
    with numbered nodes"""
    nodeDict = {}
    for n in root.iternodes("preorder"):
        if n.istip:
            continue
        nodeDict[n.number] = {
            "shift": False,
            "model": n.label,
            "type": None
        }
    return nodeDict


def count_shifts(inroot, min_clade_size, min_overlap, format="NEX"):
    # different version which compares back to parent instead of forward to
    # child
    countDict = {
        "Testable nodes": 0,
        "Testable orthologs": 0,
        "Orthologous shifts": 0,
        "Testable paralogs": 0,
        "Shifts on duplication": 0,
        "Paralogous shifts": 0
    }

    shiftDict = {}
    nodeDict = populate_node_dict(inroot)

    for n in inroot.iternodes(order="preorder"):
        if n.number == 0:  # root
            continue
        if n.istip:
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
                        nodeDict[n.number] = {
                            "shift": True,
                            "model": m,
                            "type": "OS"
                        }
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
                            nodeDict[n.number] = {
                                "shift": True,
                                "model": m,
                                "type": "SD"
                            }
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
                                nodeDict[child0.number] = {
                                    "shift": True,
                                    "model": m,
                                    "type": "PS"
                                }
                        else:  # start the sub-loop
                            for c in child0.iternodes(order="preorder"):
                                if len(c.leaves()) >= min_clade_size:
                                    if check_shift(c.parent, c, format):
                                        if (c.parent, c) in shiftDict.keys():
                                            continue
                                        else:
                                            shiftDict[(c.parent, c)] = "PS"
                                            m = get_model(c)
                                            nodeDict[c.number] = {
                                                "shift": True,
                                                "model": m,
                                                "type": "PS"
                                            }
                    if len(child1.leaves()) >= min_clade_size:
                        # paralog that can be tested
                        countDict["Testable paralogs"] += 1
                        if check_shift(n, child1, format):
                            # shift in paralog immediately
                            # below dup node
                            if (n, child1) in shiftDict.keys():
                                continue
                            else:
                                shiftDict[(n, child1)] = "PS"
                                m = get_model(child1)
                                nodeDict[child1.number] = {
                                    "shift": True,
                                    "model": m,
                                    "type": "PS"
                                }
                        else:  # start the sub-loop
                            for c in child1.iternodes(order="preorder"):
                                if len(c.leaves()) >= min_clade_size:
                                    if check_shift(c.parent, c, format):
                                        if (c.parent, c) in shiftDict.keys():
                                            continue
                                        else:
                                            shiftDict[(c.parent, c)] = "PS"
                                            m = get_model(c)
                                            nodeDict[c.number] = {
                                                "shift": True,
                                                "model": m,
                                                "type": "PS"
                                            }
                elif len(names0.intersection(names1)) == 0:
                    # not overlapping - this is node which has paralogs below
                    countDict["Testable orthologs"] += 1
                    if check_shift(n.parent, n, format):
                        if (n.parent, n) in shiftDict.keys():
                            continue
                        else:
                            shiftDict[(n.parent, n)] = "OS"
                            m = get_model(n, format)
                            nodeDict[n.number] = {
                                     "shift": True,
                                     "model": m,
                                     "type": "OS"
                            }
                elif 0 <= len(names0.intersection(names1)) <= min_overlap:
                    # overlap, but not enough to robustly call paralog
                    # call ortholog instead
                    countDict["Testable orthologs"] += 1
                    if check_shift(n.parent, n, format):
                        if (n.parent, n) in shiftDict.keys():
                            continue
                        else:
                            shiftDict[(n.parent, n)] = "OS"
                            m = get_model(n, format)
                            nodeDict[n.number] = {
                                     "shift": True,
                                     "model": m,
                                     "type": "OS"
                            }
    for _, v in shiftDict.items():
        if v == "OS":
            countDict["Orthologous shifts"] += 1
        elif v == "SD":
            countDict["Shifts on duplication"] += 1
        elif v == "PS":
            countDict["Paralogous shifts"] += 1
    return countDict, shiftDict, nodeDict


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

    number_tree(curroot)
    cDict, sDict, nDict = count_shifts(curroot, 10, 10, treForm)
    
    with open(tFile + ".shiftcount", "w") as out1:
        for k, v in cDict.items():
            out1.write(tFile + "\t" + k + "\t" + str(v) + "\n")

    with open(tFile + ".nodes", "w") as out2:
        for k, v in nDict.items():
            out2.write(tFile + "\t" + str(k) + "\t" + "\t".join(str(x) for x in v.values()) + "\n")
