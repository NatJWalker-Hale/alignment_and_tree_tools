#! /usr/bin/python3

import sys
import argparse
import newick3
import logging
from phylo3 import Node
from tree_utils import get_name, remove_kink
from sequence import parse_fasta


def get_names_to_exclude(ignoref):
    ignore = [line.strip() for line in open(ignoref, "r").readlines()]
    return ignore


def sample_id_from_name(name: str) -> str:
    id = name.split("_")[1]
    return id


def is_dup_sp_ovlp(node: Node) -> bool:
    ch1, ch2 = node.children
    ch1_samps = {sample_id_from_name(n.label) for n in ch1.leaves()}
    ch2_samps = {sample_id_from_name(n.label) for n in ch2.leaves()}
    return bool(len(ch1_samps.intersection(ch2_samps)) > 0)


def is_monophyletic_sp(node: Node) -> bool:
    sp = [get_name(n.label) for n in node.leaves()]
    return len(set(sp)) == 1


def mask_monophyletic_tips(curroot, unamb_chrDICT, ignore=[]):
    going = True
    while going and curroot is not None and len(curroot.leaves()) >= 4:
        going = False
        for node in curroot.iternodes():   # walk through nodes
            if not node.istip:
                continue   # only look at tips
            name = get_name(node.label)
            if name in ignore:
                continue   # do not mask the genomes
            if node.parent.dup:
                # print("parent dup")
                continue
            for sister in node.get_sisters():
                if sister.istip and name == get_name(sister.label):  # mask
                    if sample_id_from_name(node.label) == sample_id_from_name(sister.label):
                        continue
                    if ("_ptg" in node.label) & ("_ptg" in sister.label):
                        continue
                    if "_ptg" in sister.label:
                        # print(f"from {node.label} found hifi pruning {node.label}")
                        node = node.prune()
                    elif "_ptg" in node.label:
                        # print(f"from {node.label} found hifi pruning {sister.label}")
                        node = sister.prune()
                    else:
                        if unamb_chrDICT[node.label] > unamb_chrDICT[sister.label]:
                            # print(f"from {node.label} I'm longer than {sister.label} pruning {sister.label}")
                            node = sister.prune()
                        else:
                            # print(f"from {node.label} I'm shorter than {sister.label} pruning {node.label}")
                            node = node.prune()
                    if len(curroot.leaves()) >= 4:
                        if (
                            node == curroot
                            and node.nchildren == 2
                            ) or (
                                node != curroot
                                and node.nchildren == 1):
                            node, curroot = remove_kink(node, curroot)
                    # print(newick3.to_string(curroot)+";")
                    going = True
                    break
    return curroot


def mask_paraphyletic_tips(curroot, unamb_chrDICT, ignore=[]):
    going = True
    while going and curroot is not None and len(curroot.leaves()) >= 4:
        going = False
        for node in curroot.iternodes():  # walk through nodes
            if not node.istip:
                continue  # only look at tips
            name = get_name(node.label)
            if name in ignore:
                continue  # do not mask the genomes
            parent = node.parent
            if node == curroot or parent == curroot or parent is None:
                continue  # no paraphyletic tips for the root
            if parent.dup:
                # print("parent dup")
                continue
            for para in parent.get_sisters():
                if para.istip and name == get_name(para.label):  # mask
                    if sample_id_from_name(node.label) == sample_id_from_name(para.label):
                        continue
                    if para.parent.dup:  # don't prune a para tip where the MRCA of focal tip and 
                        # para tip is a duplication node
                        continue
                    if ("_ptg" in node.label) & ("_ptg" in para.label):
                        continue
                    if "_ptg" in para.label:
                        # print(f"from {node.label} found hifi pruning {node.label}")
                        node = node.prune()
                    elif "_ptg" in node.label:
                        # print(f"from {node.label} found hifi pruning {para.label}")
                        node = para.prune()
                    else:
                        if unamb_chrDICT[node.label] > unamb_chrDICT[para.label]:
                            # print(f"from {node.label} I'm longer than {para.label} pruning {para.label}")
                            node = para.prune()
                        else:
                            # print(f"from {node.label} I'm shorter than {para.label} pruning {node.label}")
                            node = node.prune()
                    if len(curroot.leaves()) >= 4:
                        if (
                            node == curroot
                            and node.nchildren == 2
                            ) or (
                                node != curroot
                                and node.nchildren == 1):
                            node, curroot = remove_kink(node, curroot)
                    # print(newick3.to_string(curroot)+";")
                    going = True
                    break
    return curroot


def mask(curroot, clnfile, para=True, ignore=[]):
    chrDICT = {}  # key is seqid, value is number of unambiguous chrs
    for key, value in dict([x for x in parse_fasta(clnfile)]).items():
        for ch in ['-', 'X', "x", "?", "*"]:
            value = value.replace(ch, "")  # ignore gaps, xs and Xs
        chrDICT[key] = len(value)
    curroot = mask_monophyletic_tips(curroot, chrDICT, ignore)
    if para:
        curroot = mask_paraphyletic_tips(curroot, chrDICT, ignore)
    return curroot


def mask_paralogs(tree: Node, clnaln: str):
    """
    find individual paralogues that are monophyletic from a single species and reduce these to one
    tip. Keep HiFi genome if present, otherwise longest unambiguous sequence in alignment
    """
    for node in tree.iternodes(order=0):
        if node == tree or node.parent == tree or node.istip:
            continue
        if is_dup_sp_ovlp(node):
            node.dup = True
            # node.label = "D"

    # print(newick3.to_string(tree) + ";")

    chrDICT = {}  # key is seqid, value is number of unambiguous chrs
    for key, value in dict(parse_fasta(clnaln)).items():
        for ch in ['-', 'X', "x", "?", "*"]:
            value = value.replace(ch, "")  # ignore gaps, xs and Xs
        chrDICT[key] = len(value)

    tree = mask_monophyletic_tips(tree, chrDICT)
    # print(newick3.to_string(tree) + ";")
    # for node in tree.iternodes(order=0):
    #     if node == tree or node.parent == tree or node.istip:
    #         continue
    #     if is_dup_sp_ovlp(node):
    #         node.dup = True
    #         node.label = "D"
    # print("done mono")
    tree = mask_paraphyletic_tips(tree, chrDICT)
    
    # for n in tree.iternodes(order=1):  # do postorder so that nested duplicates are masked first
    #     if n.istip:  # don't consider tips
    #         continue
    #     if n.parent is None:  # skip root
    #         continue
    #     if is_dup_sp_ovlp(n):  # skip nodes which are duplicates
    #         print("found dup node")
    #         continue
    #     # if is_monophyletic_sp(n):  # want to mask this
    #     #     print("found monophyletic node")
    #     mask_monophyletic_tips(n, chrDICT)
    #     #print(newick3.to_string(n) + ";")
    #     mask_paraphyletic_tips(n, chrDICT)
    #         ls = n.leaves()
    #         keep = None
    #         for l in ls:
    #             if "_ptg_" in l.label:
    #                 print(f"keeping hifi {l.label}")
    #                 keep = l
    #         if keep is None:
    #             node_L = sorted([(l, chrDICT[l.label]) for l in ls], key=lambda x: x[1],
    #                            reverse=True)
    #             keep = node_L[0][0]
    #         for node in n.iternodes(order=0):
    #             if node.istip:
    #                 if node is not keep:
    #                     print(f"pruning {node.label}")
    #                     node = node.prune()
    #             if (
    #                 node == tree
    #                 and node.nchildren == 2
    #                 ) or (
    #                 node != tree
    #                 and node.nchildren == 1):
    #                     node, tree = remove_kink(node, tree)
    print(newick3.to_string(tree) + ";")

# def mask_monophyly(tre, clnaln, para=True, ignore=[]):
#     with open(tre, "r", encoding='utf-8') as inf:
#         intree = newick3.parse(inf.readline())
#     in_tips = len(intree.leaves())
#     curroot = mask(intree, clnaln, para, ignore)
#     out_tips = len(curroot.leaves())
#     masked = tre + ".mm"
#     if para:
#         logging.info("masking monophyletic and paraphyletic tips in %s", tre)
#     else:
#         logging.info("masking monophyletic tips in %s", tre)
#     logging.info("masked %s tips, writing to %s", (in_tips - out_tips), masked)
#     with open(masked, "w", encoding='utf-8') as outf:
#         outf.write(newick3.tostring(curroot)+";\n")
#     return masked


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script to collapse multiple sequences from the \
                                     same species to a single tip. Currently set up for Laura \
                                     Campbell's project on Phylidris nagasau, so some hardcoding \
                                     needs to be changed for general use")
    parser.add_argument("tree",
                        help="The tree file in newick format to mask tips")
    parser.add_argument("aln_cln",
                        help="The cleaned alignment in FASTA format \
                        corresponding to the tips of the tree")
    args = parser.parse_args(sys.argv[1:] or ["--help"])
    
    curroot = newick3.parse_from_file(args.tree)
    mask_paralogs(curroot, args.aln_cln)
