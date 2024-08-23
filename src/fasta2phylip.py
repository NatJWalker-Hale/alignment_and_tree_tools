#!/usr/bin/env python3


"""
converts FASTA-formatted alignments torelaxed non-interleaved phylip
"""


import sys
import argparse
# from Bio import SeqIO
import sequence as sq


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="converts FASTA-formatted MSAs to PHYLIP format")
    parser.add_argument("aln", help="FASTA-formatted input alignment")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    seqs = dict(sq.parse_fasta(args.aln))
    print(sq.get_phylip_str(seqs))
