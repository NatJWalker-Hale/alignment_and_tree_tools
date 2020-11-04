#! /usr/bin/python3

import sys
import argparse
from parse_fasta import parse_fasta


def get_difference(seqDict):
    """for a seqDict of two sequences, finds positional differences in the two
    and creates a dictionary of position: state in the SECOND sequence"""
    diffDict = {}
    s1 = list(seqDict.values())[0]
    s2 = list(seqDict.values())[1]
    for pos, state in enumerate(s1):
        if state == s2[pos]:
            continue
        else:
            diffDict[pos] = s2[pos]
    return diffDict


def create_seqs_from_diff(seqDict, diffDict):
    """for a seqDict of two sequences and diffDict of the two, use positional
    differences to create a synthetic sequence per difference where the state
    from the SECOND sequence is introduced to the background of the first."""
    outDict = {}  # key is name of first seq + _pos_orig_state_to_new_state
    baseName = list(seqDict.keys())[0]
    background = list(seqDict.values())[0]
    for k, v in diffDict.items():
        newName = baseName+"_"+str(k+1)+str(background[k])+"_to_"+str(v)
        outDict[newName] = background[:k]+str(v)+background[k+1:]
    return outDict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sequences", help="FASTA formatted file of \
                        two sequences. The first sequence is assumed to be \
                        the background for necessity testing.")
    args = parser.parse_args()

    seqs = dict([x for x in parse_fasta(args.sequences)])
    diffs = get_difference(seqs)
    outSeqs = create_seqs_from_diff(seqs, diffs)

    for k, v in outSeqs.items():
        print(">"+str(k))
        print(v)
