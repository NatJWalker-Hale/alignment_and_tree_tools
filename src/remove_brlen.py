#! /usr/bin/python3


"""
removes branch lengths from Newick-formatted trees
"""


import re
import sys


with open(sys.argv[1]) as inf:
    t_str = inf.readline().strip()

print(re.sub(r":[0-9.E-]+", "", t_str))