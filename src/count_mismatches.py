#! /usr/bin/python3

import sys
import argparse
from parse_fasta import parse_fasta


def count_mismatches(seq1, seq2):
    matches = sum([1 if x[0] == x[1] else 0 for x in zip(seq1, seq2)])
    return matches


