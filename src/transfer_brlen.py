#! /usr/bin/python3

import sys
import argparse
import tree_reader


def transfer_brlen(tr1, tr2):
    """For two trees with identical topology, transfer branch lengths
    from tree one to tree two. Modifies in place."""
    for i in tr1.iternodes():
        s = set(i.lvsnms())
        match = False
        for j in tr2.iternodes():
            x = set(j.lvsnms())
            if s == x:
                match = True
                brlen = i.length
                j.length = brlen
                break
        if match is False:
            return False
    if match is True:
        return True


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree1", help="nwk tree to transfer brlen from")
    parser.add_argument("tree2", help="nwk tree to transfer brlen to")
    args = parser.parse_args()

    with open(args.tree1, "r") as inf:
        s = inf.readlines()[0].strip()
        tree1 = tree_reader.read_tree_string(s)

    with open(args.tree2, "r") as inf:
        s = inf.readlines()[0].strip()
        tree2 = tree_reader.read_tree_string(s)

    if not transfer_brlen(tree1, tree2):
        print("topologies are not identical!")
        sys.exit()
    else:
        print(tree2.get_newick_repr(True))
