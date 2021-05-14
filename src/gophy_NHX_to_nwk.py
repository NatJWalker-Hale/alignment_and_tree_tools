#! /usr/bin/python3.9


import sys
import re
import argparse
from ete3 import PhyloTree


def convert(nwk):
    paramDict = {}
    p1 = re.compile("\\)\\[.+?\\]")  # internal branch NHX labels
    modSet = set(p1.findall(nwk))
    for i in modSet:
        model = re.search("&model=\\d*", i).group().split("=")[1]
        params = [float(x) for x in re.search("\\{.*\\}", i).group().lstrip("{").rstrip("}").split(",")]
        paramDict[model] = params
    s1 = p1.sub(lambda m: ")" + m.group().split("=")[1].split(",")[0], nwk)
    p2 = re.compile("\\[.+?\\]")  # tip NHX labels
    s2 = p2.sub("", s1)
    return s2, paramDict


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

    newNwkString, params = convert(nwkString)
    # print(newNwkString)

    t = PhyloTree(newNwkString, format=1)
    to_change = []
    for n in t.traverse(strategy="preorder"):
        if not n.is_leaf():
            for c in n.children[:2]:
                if not c.is_leaf():
                    # print(n.name,c.name)
                    if c.name != n.name:  # transition
                        to_change.append(c)

    outTreeName = ".".join(tFile.split(".")[:-3])
    sys.stderr.write(outTreeName+"\n")
    with open(outTreeName + ".godon.tre", "w") as outF:
        outF.write(t.write(format=1)+"\n")

    count = 0
    for node in to_change:
        for n in t.traverse(strategy="preorder"):
            if not n.is_leaf():
                if n == node:
                    count += 1
                    n.name = "#1"
                else:
                    n.name = ""
        with open(outTreeName + ".godon.tre" + str(count), "w") as newTreeFile:
            newTreeFile.write(t.write(format=1)+"\n")

    # count = 0
    # for n in t.traverse(strategy="preorder"):
    #     if not n.is_leaf():
    #         if n in to_change:
    #             count += 1
    #             n.name = "#" + str(count)
    #             with open(outTreeName + ".godon.tre" + str(count),
    #                       "w") as newTreeFile:
    #                 newTreeFile.write(t.write(format=1)+"\n")
    #         else:
    #             n.name = ""

    # outTreeName = tFile.rstrip(".gophy.results.tre")
    # with open(outTreeName + ".godon.tre", "w") as newTreeFile:
    #     newTreeFile.write(t.write(format=1)+"\n")

    with open(outTreeName + ".params.txt", "w") as paramFile:
        paramFile.write("model\tA\tC\tG\tT\n")
        for k, v in params.items():
            paramFile.write(str(k)+"\t"+"\t".join([str(x) for x in v])+"\n")

    
