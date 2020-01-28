#! /usr/bin/python3

import sys
import argparse

def parse_fasta(path): # courtesy of Jonathan Chang https://gist.github.com/jonchang/6471846
    """Given a path tries to parse a fasta file. Returns an iterator which
    yields a (name, sequence) tuple"""
    with open(path) as handle:
        name = sequence = ""
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, sequence
                name = line[1:]
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence

def get_pos(fa_dict):
    pos_dict = {}
    for key, value in fa_dict.items():
        if len(value) % 3 != 0:
            print(key+" CDS length not evenly divisible by three. Exiting.")
            break
        pos_dict[key] = (value[0::3],value[1::3],value[2::3])
    return pos_dict

if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--position",
    help="Select which codon position(s) to extract from the alignment. Any combo of 123 accepted. Defaults to 123, i.e. the original alignment reordered.",default=123)
    parser.add_argument("codon_aln",help="The input sequences in FASTA aligned by codon. Must be evenly divisible by 3.")
    args = parser.parse_args()

    fa_dict = dict([x for x in parse_fasta(args.codon_aln)])
    pos_dict = get_pos(fa_dict)
    #print(pos_dict) # debug
    pos = sorted([int(x)-1 for x in str(args.position)])
    for key,value in pos_dict.items():
        print(">"+key)
        for p in pos:
            print(value[p],end="")
        print("\n",end="")
    