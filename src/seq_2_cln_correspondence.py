#! /usr/bin/python3

import sys, os, subprocess
from collections import Counter

def parse_fasta(path): # courtesy of Jonathan Chang https://gist.github.com/jonchang/6471846
    """Given a path tries to parse a fasta file. Returns an iterator which
    yields a (name, sequence) tuple"""
    with open(path) as handle:
        name = sequence = ""
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, sequence
                name = line[1:]
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence

def get_site_dict(seq_dict):
    """Given a seq dict from an alignment, create a dictionary with columns as values"""
    site_dict = {}
    s = 0
    while s < len([v for v in seq_dict.values()][0]):
        site_dict[s] = []
        for v in seq_dict.values():
            site_dict[s].append(v[s])
        s += 1
    return site_dict
        
def get_if_cleaned_dict(site_dict, prop):
    """Given a site dict, replace the value of each site with True if the site is cleaned and False if otherwise"""
    if_cleaned_dict = {}
    amb_list = ["-","?","X"]
    for key, value in site_dict.items():
        unambs = sum([v for k,v in Counter(value).items() if k not in amb_list])
        #print(unambs)
        if unambs > len(value)*prop:
            if_cleaned_dict[key] = False
        else:
            if_cleaned_dict[key] = True
    return if_cleaned_dict

def get_site_correspondence(if_cleaned_dict):
    """Using if_cleaned_dict, get correspondence between sites in an alignment and sites in a cleaned alignment"""
    cln_len = sum(1 for v in if_cleaned_dict.values() if not v)
    # print(cln_len)
    site_corres_dict = {}
    s = 0
    c = 0
    while s < cln_len:
        if if_cleaned_dict[s]:
            c += 1
            site_corres_dict[s] = c
        else:
            site_corres_dict[s] = c
        s += 1
        c += 1
    return site_corres_dict

def get_site_correspondence_single_seq(seq_dict,seqn):
    corres_dict = {}
    seq = seq_dict[seqn]
    s = 0
    c = 0
    while s < len(seq):
        if seq[s] == "-":
            corres_dict[s] = ""
        else:
            corres_dict[s] = c
            c += 1
        s +=1
    return corres_dict

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: python "+sys.argv[0]+" aln seq prop[0-1]")
        sys.exit()
    
    seq_dict = dict([x for x in parse_fasta(sys.argv[1])])
    #print(seq_dict)
    site_dict = get_site_dict(seq_dict)
    #print(site_dict)
    if_cleaned_dict = get_if_cleaned_dict(site_dict,float(sys.argv[3]))
    #print(if_cleaned_dict)
    seq_corres = get_site_correspondence_single_seq(seq_dict,sys.argv[2])
    seq_cln_corres = {}
    pos = 0
    for k,v in if_cleaned_dict.items():
        if v:
            continue
        else:
            if type(seq_corres[k]) == int:
                seq_cln_corres[pos+1] = seq_corres[k]+1
            else: 
                seq_cln_corres[pos+1] = seq_corres[k] 
            pos += 1
    for k,v in seq_cln_corres.items():
        print(str(k)+"\t"+str(v))



