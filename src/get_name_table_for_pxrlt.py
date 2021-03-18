import sys
import argparse


def parse_name_file(path):
    nameDict = {}
    with open(path, "r") as inf:
        for line in inf.readlines():
            line = line.strip()
            nameDict[line.split("\t")[0]] = line.split("\t")[1]
    return nameDict


if __name__ == "__main__":
    if len(sys.argv[:1]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("nameKey", help="tab delimited file of code\tname, \
                        1 per line")
    parser.add_argument("names", help="file of sequence names to change, \
                        formatted code@locus, 1 per line")
    args = parser.parse_args()

    names = parse_name_file(args.nameKey)
    print(names)

    seqs = {}
    with open(args.names, "r") as inseq:
        for line in inseq.readlines():
            if "@" not in line:
                continue
            line = line.strip()
            taxa = line.split("@")[0]
            locus = line.split("@")[1]
            try:
                seqs[taxa+"@"+locus] = names[taxa]+"@"+locus
            except KeyError:
                seqs[taxa+"@"+locus] = taxa+"@"+locus

    with open("old_names.txt", "w") as outOld:
        for k, _ in seqs.items():
            outOld.write(k+"\n")
    with open("new_names.txt", "w") as outNew:
        for _, v in seqs.items():
            outNew.write(v+"\n")
   
