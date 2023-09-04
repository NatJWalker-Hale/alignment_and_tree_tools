#! /usr/bin/python3


from chimerax.core.commands import run


def color_resi_conv_coin(session):
    typeDict = {}
    with open("/home/nat/Dropbox/cary_projects/DODA/figures/convergence_portfolio/manuscript_versions/element_files/20230518_conv_and_coin_subs_N42.tsv", "r") as inf:
        # next(inf)  # skip header
        for line in inf:
            line = line.strip().split("\t")
            typeDict[line[1]] = line[5]
    for k, v in typeDict.items():
        if v == "CONV":
            cmd = f"color /A:{k} #D40000 target s"
        else:
            cmd = f"color /A:{k} #0055D4 target s"
        run(session, cmd)


color_resi_conv_coin(session)
