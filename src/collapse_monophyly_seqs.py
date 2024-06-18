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
            for sister in node.get_sisters():
                if sister.istip and name == get_name(sister.label):  # mask
                    if sample_id_from_name(node.label) == sample_id_from_name(sister.label):
                        continue
                    if ("_ptg" in node.label) & ("_ptg" in sister.label):
                        continue
                    if "_ptg" in sister.label:
                        node = node.prune()
                    elif "_ptg" in node.label:
                        node = sister.prune()
                    else:
                        if unamb_chrDICT[node.label] > unamb_chrDICT[sister.label]:
                            node = sister.prune()
                        else:
                            node = node.prune()
                    if len(curroot.leaves()) >= 4:
                        if (
                            node == curroot
                            and node.nchildren == 2
                            ) or (
                                node != curroot
                                and node.nchildren == 1):
                            node, curroot = remove_kink(node, curroot)
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
            for para in parent.get_sisters():
                if para.istip and name == get_name(para.label):  # mask
                    if sample_id_from_name(node.label) == sample_id_from_name(para.label):
                        continue
                    if ("_ptg" in node.label) & ("_ptg" in para.label):
                        continue
                    if "_ptg" in para.label:
                        node = node.prune()
                    elif "_ptg" in node.label:
                        node = para.prune()
                    else:
                        if unamb_chrDICT[node.label] > unamb_chrDICT[para.label]:
                            node = para.prune()
                        else:
                            node = node.prune()
                    if len(curroot.leaves()) >= 4:
                        if (
                            node == curroot
                            and node.nchildren == 2
                            ) or (
                                node != curroot
                                and node.nchildren == 1):
                            node, curroot = remove_kink(node, curroot)
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


def mask_paralogs(tree: Node, clnaln: str, para=True, ignore=[]):
    chrDICT = {}  # key is seqid, value is number of unambiguous chrs
    for key, value in dict(parse_fasta(clnaln)).items():
        for ch in ['-', 'X', "x", "?", "*"]:
            value = value.replace(ch, "")  # ignore gaps, xs and Xs
        chrDICT[key] = len(value)
    for n in tree.iternodes(order=0):
        if n.istip:  # don't consider tips
            continue
        if n.parent is None:  # skip root
            continue
        if is_dup_sp_ovlp(n):  # skip nodes which are duplicates
            continue
        if is_monophyletic_sp(n):  # want to mask this
            ls = n.leaves()
            keep = None
            for l in ls:
                if "_ptg" in l.label:
                    keep = l
            if 




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
    parser.add_argument("-p", "--mask_paraphyly",
                        help="Mask tips that are paraphyletic (default True)",
                        default=True)
    parser.add_argument("tree",
                        help="The tree file in newick format to mask tips")
    parser.add_argument("aln_cln",
                        help="The cleaned alignment in FASTA format \
                        corresponding to the tips of the tree")
    args = parser.parse_args(sys.argv[1:] or ["--help"])
    if args.mask_paraphyly:
        para = True
    else:
        para = False
    if args.exclude_file is not None:
        ignore = get_names_to_exclude(args.exclude_file)
    else:
        ignore = []
    _ = mask_monophyly(args.tree, args.aln_cln, para, ignore)
