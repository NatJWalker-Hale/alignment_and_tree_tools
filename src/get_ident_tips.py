#! /usr/bin/python

import sys
import argparse
import phylo3
import newick3
from random import choice


def get_identical_tips(root: newick3.Node) -> list[tuple]:
    idents = []
    for n in root.iternodes():
        if n.istip:
            if True in [n.label in i for i in idents]:  # avoid double counting
                continue
            if n.length <= 0.000001:  # this might need to change according to software
                ch2 = n.get_sisters()
                if ch2[0].istip:
                    idents.append((n.label, ch2[0].label))
    return idents


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")


    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help = "Tree to remove identical tips")
    parser.add_argument("-p", "--pick", help = "Pick one randomly out of each identical duo",
                        action = "store_true")
    args = parser.parse_args()


    curroot = newick3.parse_from_file(args.tree)


    ident_tips = get_identical_tips(curroot)
    print(len(ident_tips))
    if len(ident_tips) == 0:
        sys.stderr.write("There are no identical tips, exiting\n")
        sys.exit()
    else:
        sys.stderr.write("There are %s identical tips\n" % len(ident_tips))

    if args.pick:
        for i in ident_tips:
            print(choice(i))
    else:
        for i in ident_tips:
            print(i[0] + "\t" + i[1])



    

