import sys,os
from ete3 import Tree,PhyloTree

# script to return a list of subsampled sequences based on genera in paralogs, such that each genus is only represented once. All names must be in format family_genus_species@seqid and tip names must
# match alignment names exactly. Alignment in fasta format by default, this can easily be changed.

t = PhyloTree(sys.argv[1],alignment=sys.argv[2],alg_format="fasta") # change format for different alignment formats
t.set_species_naming_function(lambda node: node.name.split("@")[0].rstrip("_SRA") )

# first, get rid of monophyletic tips from all the same species

dtip = []
species_list = list(set([ leaf.species for leaf in t.iter_leaves() ]))
for s in species_list:
    masks = [ n for n in t.get_monophyletic(values=[s],target_attr="species") if len(n) > 1 ]
    nlist = []
    for n in masks:
        for l in n.iter_leaves():
            nlist.append((l.name,len(l.sequence.replace("-",""))))
        dtip += [ x[0] for x in sorted(nlist,reverse=True,key = lambda x:x[1])[1:] ]
        nlist = []

#print len(dtip)
tips = set( [ l.name for l in t.iter_leaves() ] )
#print len(tips)
ktips = list(tips.symmetric_difference(set(dtip)))
#print len(ktips)
t.prune(ktips)
#print t 

dups_to_process = t.split_by_dups()

keep_seqs = []

for node in dups_to_process:
    genus_dict = {} #key is genus, value is list of tuples
    for leaf in node.iter_leaves():
        if not leaf.name.split("_")[1] in genus_dict.keys():
            genus_dict[leaf.name.split("_")[1]] = [(leaf.name,len(leaf.sequence.replace("-","")))]
        else: 
            genus_dict[leaf.name.split("_")[1]].append((leaf.name,len(leaf.sequence.replace("-",""))))
    for genus in genus_dict.keys():
        keep_seqs.append(max(genus_dict[genus], key=lambda x:x[1])[0])

print "\n".join(keep_seqs)     
    
