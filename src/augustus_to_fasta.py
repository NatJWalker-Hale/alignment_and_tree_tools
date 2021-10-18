#! /usr/bin/python3

import sys
import re
import argparse


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("augustusOut", help="augustus prediction output")
    args = parser.parse_args()

    p = re.compile("name = .*\\)")
    with open(args.augustusOut, "r") as inf:
        going = False
        curSeq = ""
        for line in inf:
            if line.startswith("# ----- "):
                res = re.search(p, line).group()
                ctig = res.split(" = ")[1].rstrip(")")
            if line.startswith("# start gene"):
                gene = line.split(" ")[-1].strip()
            if line.startswith("# protein sequence"):
                if line.endswith("]\n"):  # short, on single line
                    line = line.split("[")[1].rstrip("]\n")
                    curSeq += line
                    print(">" + ctig + "." + gene)
                    print(curSeq)
                    curSeq = ""
                else:
                    going = True
                    line = line.split("[")[1].strip()
                    curSeq += line
                    continue
            if going:
                if line.endswith("]\n"):
                    going = False
                    line = line.lstrip("# ").rstrip("]\n")
                    curSeq += line
                    print(">" + ctig + "." + gene)
                    print(curSeq)
                    curSeq = ""
                else:
                    line = line.lstrip("# ").strip()
                    curSeq += line
