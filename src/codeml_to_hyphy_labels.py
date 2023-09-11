#! /usr/bin/python3


import sys
import argparse
import tree_reader


def parse_and_relable(treefile: str) -> str:
    t = [x for x in tree_reader.read_tree_file_iter(treefile)][0]
    # print(t)
    for n in t.iternodes("preorder"):
        if n.parent is None:
            continue
        if not n.istip:
            try:
                if "$" in n.label:
                    n.label = "{Foreground}"
                    for c in n.iternodes("preorder"):
                        if c.istip:
                            c.label += "{Foreground}"
                        else:
                            c.label = "{Foreground}"
                if "#" in n.label:
                    n.label = "{Foreground}"
            except TypeError:
                continue
        else:
            if n.label.endswith("#1"):
                seqname = n.label.split("#")[0]
                n.label = seqname + "{Foreground}"
    print(t.get_newick_repr(showbl=True) + ";")


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="tree with paml labels to reformat")
    args = parser.parse_args()

    parse_and_relable(args.tree)
