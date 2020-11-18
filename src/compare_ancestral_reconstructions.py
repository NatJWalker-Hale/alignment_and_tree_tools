#! /usr/bin/python3

import sys
import argparse
from parse_fasta import parse_fasta

"""This script is designed to take alignments of sequences for the same (or comparable) node from different 
ancestral sequence reconstructions, and compute summary statistics of the differences. It presently expects 
two sequences per node - a MAP sequence and an AltAll sequence, and expects the FASTA to be ordered
accordingly."""

def compute_diffs(seqDict):
    """This function takes a sequence dictionary of four sequences from two different analyses and determines
    the number and character of differences between them. It outputs two dictionaries: one with summary
    statistics, and one with each position that differs, showing states for each of the sequences."""
    summDict = {"Total diffs": 0, "Single seq diffs": 0, "Fixed (analysis) diffs": 0,
                "Fixed (type) diffs": 0, "All seq diffs": 0}
    diffDict = {}
    seqs = list(seqDict.values())
    # relies on python3 dicts being insertion-ordered, be wary
    for pos, _ in enumerate(seqs[0]):
        states = [x[pos] for x in seqs]
        if len(set(states)) <= 1:  # all same
            continue
        else:
            summDict["Total diffs"] += 1
            diffDict[pos] = (states, "")
            if len(set(states[0:2])) <= 1 or len(set(states[2:4])) <= 1:  # fixed analysis difference
                summDict["Fixed (analysis) diffs"] += 1
                diffDict[pos][1] = "Fixed (analysis) diffs"
            elif len(set(states[0]+states[2])) <= 1 or len(set(states[1]+states[3])) <= 1:  # fixed type difference
                summDict["Fixed (type) diffs"] += 1
                diffDict[pos][1] = "Fixed (type) diffs"
            elif len(set(states)) == 4:  # all diff
                summDict["All seq diffs"] += 1
                diffDict[pos][1] = "All seq diffs"
            else:
                summDict["Single seq diffs"] += 1
                diffDict[pos][1] = "Single seq diffs"
    return summDict, diffDict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sequence", help="FASTA formatted multiple alignment of two MAP \
                        and two AltAll, ordered MAP AltAll MAP AltAll")

