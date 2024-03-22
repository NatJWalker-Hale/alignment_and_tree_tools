#!/usr/bin/env python3


"""
Parses phylip-formatted sequence data
"""


import re


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
        next(inf, (None, None))
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


def parse_phylip_str(phy_str: str):
    """
    parse relaxed PHYLIP from string. Will fail with interleaved
    """
    lines = phy_str.strip().split("\n")
    header = seq = ""
    for line in lines[1:]:
        line = line.strip()
        if bool(re.search(r"\s", line)):  # sequence header
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


def write_phylip_str(seq_dict) -> str:
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


if __name__ == "__main__":
    s = """
         4 8
        1   ATGC
        ATGC
        2   ATGC
        ATGC
        3   ATGC
        ATGC
        4   ATGC
        ATGC
        """
    print(s)
    seqs = dict(parse_phylip_str(s))
    print(seqs)
    # print(write_phylip_str(seqs))

    seqs = dict(parse_phylip("tmp_multiline_phylip.phy"))
    print(seqs)
