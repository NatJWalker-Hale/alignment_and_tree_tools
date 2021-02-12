import sys
import argparse
from parse_fasta import parse_fasta
from count_diff_pair_aln import count_diff

"""This script takes an alignment of any number of sequences
representing either a pairof MAP sequences from the same nodes
of two different ASRs,or a tetrad of sequences representing
the MAP and AltAll from the same nodes of two different ASRs.
The script expects that the sequence names be identical apart
from a prefix identifying which ASR the sequences are derived
from and a suffix indicating if the sequence is MAP or AltAll,
in format ASR_nodeName_MAP/AltAll"""


def parse_seq_dict_to_differences(seqDict):
    nodeSeqs = {}
    for k, v in seqDict.items():
        nodeName = "_".join(k.split("_")[1:-1])
        # okay for name to have underscores
        recon = k.split("_")[0]
        # recon prefix should not have underscores
        seqType = k.split("_")[-1]
        # MAP or AltAll
        if nodeName in nodeSeqs:
            if recon in nodeSeqs[nodeName]:
                nodeSeqs[nodeName][recon].append(v)
            else:
                nodeSeqs[nodeName][recon] = [v]
        else:
            nodeSeqs[nodeName] = {recon: [v]}
    
    for k, v in nodeSeqs.items():



if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    # parser.add_argument("tree", help="newick formatted tree, with ancestral \
    #                     node labels")
    parser.add_argument("sequences", help="FASTA formatted \
                        alignment with ASR prefix and MAP/AltAll \
                        suffix")
    parser.add_argument("-alt", "--altall", help="Additionally compare \
                        AltAll sequences (default False)", type=bool,
                        default=False)
    args = parser.parse_args()

    seqs = dict([x for x in parse_fasta(args.sequences)])
    nodes = parse_seq_dict_to_differences(seqs)
    print(nodes)