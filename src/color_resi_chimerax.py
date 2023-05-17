#! /usr/bin/python3


from chimerax.core.commands import run


def color_resi(session):
    # change absolute path here
    colDict = {}
    with open("BvDODAa2_col_out.tsv", "r") as inf:
        next(inf)  # skip header
        for line in inf:
            line = line.split("\t")
            colDict[line[0]] = line[1]
    for k, v in colDict.items():
        cmd = f"color /A:{k} {v} target s"
        run(session, cmd)


color_resi(session)
#CEA500
#B3B3B3