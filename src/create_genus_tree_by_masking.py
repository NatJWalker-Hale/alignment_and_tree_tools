import newick3,phylo3,os,sys
from tree_utils import get_name,remove_kink
from seq import read_fasta_file

def mask_monophyletic_tips(curroot,ignore=[]):
	going = True
	while going and curroot != None and len(curroot.leaves()) >= 4:
		going = False
		for node in curroot.iternodes(): # walk through nodes
			if not node.istip: continue # only look at tips
			name = get_name(node.label).split("_")[1]	
			for sister in node.get_sisters():
				if sister.istip and name==get_name(sister.label).split("_")[1]: # mask
					node = sister.prune()
					if len(curroot.leaves()) >= 4:
						if (node==curroot and node.nchildren==2) or (node!=curroot and node.nchildren==1):
							node,curroot = remove_kink(node,curroot)
					going = True
					break
	return curroot

def mask_paraphyletic_tips(curroot,ignore=[]):
	going = True
	while going and curroot != None and len(curroot.leaves()) >= 4:
		going = False
		for node in curroot.iternodes(): #walk through nodes
			if not node.istip: continue #only look at tips
			name = get_name(node.label).split("_")[1]
			parent = node.parent
			if node == curroot or parent == curroot or parent == None:
				continue #no paraphyletic tips for the root
			for para in parent.get_sisters():
				if para.istip and name==get_name(para.label).split("_")[1]: # mask
					node = para.prune()	
					if len(curroot.leaves()) >= 4:
						if (node==curroot and node.nchildren==2) or (node!=curroot and node.nchildren==1):
							node,curroot = remove_kink(node,curroot)
					going = True
					break
	return curroot

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "usage: python "+sys.argv[0]+" treefile para(y/n)"
        sys.exit()

    intree = newick3.parse(open(sys.argv[1],"r").readline())
    masked = mask_monophyletic_tips(intree,ignore=[])
    if sys.argv[2] == "y":
        masked = mask_paraphyletic_tips(masked,ignore=[])
    print newick3.tostring(masked)+";\n" 
