#! /usr/bin/env python3


import sys
import argparse
import obonet
from goatools.anno.gaf_reader import GafReader


def main():
    parser = argparse.ArgumentParser(description="takes a tab separated file with a field of comma \
                                     separated GO terms and writes a formatted output")
    parser.add_argument("input", help="input tsv (e.g. eggNOG)")
    parser.add_argument("obo", help="GO .obo, e.g. go-basic.obo")
    parser.add_argument("blast", help="tab separated annotations from blast, when eggNOG is \
                        missing")
    parser.add_argument("gaf", help="GO annotation file in GAF format, to get GO terms from blast")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    graph = obonet.read_obo(args.obo)

    with open(args.blast, "r", encoding="utf-8") as blast_inf:
        blast_annotations = {}
        for line in blast_inf:
            line = line.strip().split("\t")
            locus = line[0]
            blast_annotations[locus] = line[1:]

    ogaf = GafReader(args.gaf)
    ns2assc = ogaf.get_ns2assc()

    ns_map = {
    "molecular_function": "molecular_function",
    "biological_process": "biological_process",
    "cellular_component": "cellular_component",
    }

    print("gene\tGO_id\tGO_name\tGO_class", file=sys.stdout)
    with open(args.input, "r", encoding="utf-8") as inf:
        for line in inf:
            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            locus = line[0]
            if "-" in line[9]:
                GOs = []
                try:
                    if blast_annotations[locus][2] is not None:
                        acc = blast_annotations[locus][2].split("|")[0]
                        for ns, assc in ns2assc.items():
                            if acc in assc:
                                GOs += assc[acc]
                                break
                    else:
                        if blast_annotations[locus][3] is not None:
                            acc = blast_annotations[locus][3].split("|")[1]
                            for ns, assc in ns2assc.items():
                                if acc in assc:
                                    GOs += assc[acc]
                                    break
                except (KeyError, IndexError):
                    continue
            else:
                GOs = line[9].split(",")
            for go in GOs:
                if go in graph:
                    node = graph.nodes[go]
                    name = node.get("name", "unknown")
                    namespace = ns_map.get(node.get("namespace", ""), "unknown")
                    print(f"{locus.rstrip(".1")}\t{go}\t{name}\t{namespace}", file=sys.stdout)


if __name__ == "__main__":
    main()
