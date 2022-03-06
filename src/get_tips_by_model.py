#! /usr/bin/python3

import sys
import argparse
import tree_reader


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="gophy-formatted output tree, with model \
                        label per node")
    args = parser.parse_args()

    t = [x for x in tree_reader.read_tree_file_iter(args.tree)][0]

    print("tip,model")
    for n in t.iternodes():
        if n.istip:
            model = n.parent.label
            print(n.label + "," + str(model))

        

