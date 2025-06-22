#! /usr/bin/python3


import sys
import argparse
import sequence as sq


if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("sequence", help="FASTA-formatted sequences for upload to NCBI. Should be \
                        formatted according to BankIt submission: >Seq1 [organism=genus species] \
                        sequence description, complete cds. Currently hardcoded to write complete \
                        CDS features, where the product field will be the sequence description")
    parser.add_argument("-t", "--table", help="translation table", type=int, default=1)
    parser.add_argument("-n", "--note", help="string to add as note", type=str)

    args = parser.parse_args(sys.argv[1:] or ["--help"])

    seqs = dict(sq.parse_fasta(args.sequence))

    for header, seq in seqs.items():
        seqid = header.split("[")[0].strip()
        desc = header.split("]")[1].split(",")[0].strip()
        out_str = (f">Feature\t{seqid}\n"
                   f"1\t{len(seq)}\tCDS\n"
                   f"\t\t\tProduct\t{desc}\n"
                   "\t\t\tcodon_start\t1\n"
                   f"\t\t\tnote\t{args.note} {seqid}\n"
                   f"\t\t\ttransl_table\t{args.table}\n\n"
                  )
        
        print(out_str)