#! /usr/bin/python3


from chimerax.core.commands import run


def color_resi(session):
    colDict = {}
    # change absolute path here
    with open("/home/nat/Dropbox/cary_projects/DODA/figures/convergence_portfolio/manuscript_versions/element_files/20230530_N224_conv_coin_JSD_col.tsv", "r") as inf:
        next(inf)  # skip header
        for line in inf:
            line = line.split("\t")
            colDict[line[0]] = line[1]
    for k, v in colDict.items():
        cmd = f"color /A:{k} {v} target sc"
        run(session, cmd)


color_resi(session)
# CEA500
# B3B3B3
