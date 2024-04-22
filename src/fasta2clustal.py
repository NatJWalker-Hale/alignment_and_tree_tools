#!/usr/bin/env python3


"""
converts FASTA-formatted alignments to clustal
"""


import sys
import argparse
# from Bio import SeqIO
import sequence as sq


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="converts FASTA-formatted MSAs to CLUSTAL format")
    parser.add_argument("aln", help="FASTA-formatted input alignment")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    # seq = SeqIO.parse(args.aln, format="fasta")
    # outstr = f"{'.'.join(args.aln.split('.')[:-1])}.clustal"
    # count = SeqIO.write(seq, outstr, format="clustal")
    # sys.stderr.write(f"wrote {count} sequences to {outstr}\n")

    seqs = dict(sq.parse_fasta(args.aln))
    print(sq.write_clustalw(seqs))
    