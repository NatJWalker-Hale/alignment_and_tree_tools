#! /usr/bin/python3

import sys
import os
import argparse 
import newick3
import phylo3
from tree_utils import *

if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        print("usage: python3 "+sys.argv[0]+" treefile cutoff")
        sys.exit()

    with open(sys.argv[1],"r") as infile:
        root = newick3.parse(infile.readline())

    count = 0
    for i in root.iternodes(order=1): # postorder 
        if i.nchildren == 0: # at tip
            i.data['len'] = i.length
            if i.length > float(sys.argv[2]):
                count += 1

    if count > 0:
        print(sys.argv[1]+" flagged with "+str(count)+" tips greater than "+str(sys.argv[2]))

    
