#!/usr/bin/env python3


"""
simple script to remove stop codons from phylip-formatted alignments. Will not work for interleaved
"""


import re
import sys
import argparse

STOP = ["TAA","taa","TGA","tga","TAG","tag"] # for case insensitivity


def check_aligned(seq_dict: dict) -> bool:
    """
    checks if sequences in a sequence dictionary are the same length
    """
    if len(set(len(s) for s in seq_dict.values())) > 1:
        return False
    return True


def parse_phylip(path: str):
    """
    parse relaxed PHYLIP from file. Will fail with interleaved
    """
    with open(path, "r", encoding="utf-8") as inf:
        next(inf)
        header = seq = ""
        for line in inf:
            line = line.strip()
            if bool(re.search(r"\s", line)):
                line = line.split()
                if header:
                    yield header, seq
                    header = line[0]
                    seq = line[1]
                else:
                    header = line[0]
                    seq = line[1]
            else:
                seq += line.strip()
        yield header, seq


def write_phylip_str(seq_dict: dict) -> str:
    """
    writes a PHYLIP-formatted string from an aligned sequence dictionary {header: sequence}
    """
    if not check_aligned(seq_dict):
        raise ValueError("sequences are not aligned, write to FASTA instead")
    nseq = len(seq_dict)
    seql = set(len(v) for v in seq_dict.values()).pop()
    out = ""
    out += f" {nseq} {seql}\n"
    for header, seq in seq_dict.items():
        out += f"{header}\t{seq}\n"
    return out


def remove_stop(seq):
    i = 0
    while i < len(seq):
        codon = seq[i:i+3]
        if codon in STOP:
            if i+3 < len(seq):
                seq = seq[:i]+"---"+seq[i+3:]
            elif i+3 == len(seq):
                seq = seq[:i]+"---"
        i += 3
    return seq

if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("sequence",help="Codon sequence file.")
    args = parser.parse_args()

    seqs = dict(parse_phylip(args.sequence))
    for key, value in seqs.items():
        if len(value) % 3 != 0:
            print("Sequence length not evenly divisible by 3. Check file.")
            sys.exit()
        seqs[key] = remove_stop(value)
    
    print(write_phylip_str(seqs))
    