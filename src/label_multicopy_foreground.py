#! /usr/bin/python3


"""
script to label trees multicopy gene trees for codon model analysis
"""


import sys
import argparse
from copy import deepcopy
import newick3


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to label multicopy gene trees for codon \
                                     model analysis")
    parser.add_argument("tree", help="tree topology to label. Can have branch lengths. Internal \
                        node labels will be overwritten if present. Expects tip labels of the \
                        form taxon@seqid")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-hy", "--hyphy", help="hyphy-style labels for aBS-REL, RELAX",
                       action="store_const", dest="type", const="hyphy", default="hyphy")
    group.add_argument("-p", "--paml", help="codeml-style labels for branch models", dest="type",
                       action="store_const", const="paml")
    parser.add_argument("-f", "--flag", help="taxon name to use as foreground (code before @)")
    parser.add_argument("-b", "--brlen", help="write branch lengths. Branch lengths will be \
                        omitted by default", action="store_true")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    curroot = newick3.parse_from_file(args.tree)
    if not args.brlen:
        for n in curroot.iternodes(order=0):
            n.length = None
    sing_lab = deepcopy(curroot)

    for n in curroot.iternodes(order=0):  # first pass, all foreground, rest background
        if args.type == "hyphy":
            if n.istip:
                name = n.label.split("@")[0].split("_")[0]
                if name == args.flag:
                    n.label += "{Test}"
                else:
                    n.label += "{Reference}"
            elif n.is_monophyletic(sep="@", spl="_", query=args.flag):
                n.label = "{Test}"
            else:
                n.label = "{Reference}"
        else:
            if n.istip:
                name = n.label.split("@")[0].split("_")[0]
                if name == args.flag:
                    n.label += "#1"
            if n.is_monophyletic(sep="@", spl="_", query=args.flag):
                n.label = "#1"
    with open(args.tree + ".all.label", "w", encoding="utf-8") as outf:
        outf.write(newick3.to_string(curroot) + ";\n") 
    
    count = 0
    for n in sing_lab.iternodes(order=0):  # second pass, label one by one, write, reset label
        if n.parent is None:  # skip root
            continue
        if n.istip:
            name = n.label.split("@")[0].split("_")[0]
            if name == args.flag:
                if args.type == "hyphy":
                    n.label += "{Test}"
                else:
                    n.label += "#1"
                with open(f"{args.tree}.{count}.label", "w", encoding="utf-8") as outf:
                    outf.write(newick3.to_string(sing_lab) + ";\n")
                if args.type == "hyphy":
                    n.label = n.label[:-6]
                else:
                    n.label = n.label[:-2]
                count +=1
        elif n.is_monophyletic(sep="@", spl="_", query=args.flag):
            if args.type == "hyphy":
                n.label = "{Test}"
            else:
                n.label = "#1"
            with open(f"{args.tree}.{count}.label", "w", encoding="utf-8") as outf:
                    outf.write(newick3.to_string(sing_lab) + ";\n")
            if args.type == "hyphy":
                n.label = ""
            else:
                n.label = ""
            count +=1

