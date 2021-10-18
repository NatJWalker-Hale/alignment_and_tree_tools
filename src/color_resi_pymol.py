#! /usr/bin/python3


import pandas as pd
from pymol import cmd


def color_resi():
    # change absolute path here
    df = pd.read_csv("CgDODAa1_conv_subs_blosum62.tsv", header=0, sep="\t")
    rgbList = [[x, y, z] for x, y, z in zip(df['r'], df['g'], df['b'])]
    colDict = dict([x, y] for x, y in zip(df['pos'], rgbList))
    for k, v in colDict.items():
        colorName = "mycol" + str(k)
        cmd.set_color(colorName, v)
        cmd.color(selection="(resi " + str(k) + ")", color=colorName)


cmd.extend("color_resi", color_resi())
