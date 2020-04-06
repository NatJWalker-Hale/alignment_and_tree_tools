from __future__ import division
import sys,os
from numpy import median 

def get_num_trees(splitfile):
    trees = 0
    with open(splitfile,"r") as inf:
        for line in inf:
            if line.startswith("("):
                trees += 1 
    return trees

if __name__ == "__main__":
    splitfilelist = [f for f in os.listdir(os.getcwd()) if "splits" in f]
    total = 0
    num_tree_list = []
    outf = open("GARD_summary.txt","w")
    for f in splitfilelist:
        num_tree_list.append(get_num_trees(f))
        outf.write(f+"\t"+str(get_num_trees(f))+"\n")
        total += get_num_trees(f)
    average = total / len(splitfilelist)
    med = median(num_tree_list)
    outf.write("Average\t"+str(average)+"\n")
    outf.write("Median\t"+str(med)+"\n")
    outf.close()
        
