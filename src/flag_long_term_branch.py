#! /usr/bin/python3

import sys
import os
import argparse 
from ete3 import Tree,PhyloTree

if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        print("usage: python3 "+sys.argv[0]+" treefile cutoff")
        sys.exit()
    count = 0
    t = Tree(sys.argv[1],format=0)
    for l in t.iter_leaves():
        if l.dist > float(sys.argv[2]):
            count += 1

    if count > 0:
        print(sys.argv[1]+" flagged with "+str(count)+" tips greater than "+str(sys.argv[2]))