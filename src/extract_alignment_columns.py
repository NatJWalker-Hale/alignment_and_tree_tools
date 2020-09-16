import sys,os

## extracts alignment columns given in args.

def get_seqdict(onelnaln):
    seqdict = {} # key is sequence, value is string of characters at site given by collist
    with open(onelnaln,"r") as f:
        for line in f:
            if line.startswith(">"):
                seqdict[line.lstrip(">").rstrip("\n")] = f.next()
    return seqdict


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: python "+sys.argv[0]+ "aln col1 col2 col3 ..."
        sys.exit(0)

## first, put each sequence on one line

    with open(sys.argv[1],"r") as alignment:
        with open("aln_oneline.tmp","w") as temp_aln:
            currentline=""
            for line in alignment:
                if line.startswith(">"):
                    line = line.rstrip("\n")
                    if currentline != "": temp_aln.write(currentline + "\n")
                    temp_aln.write(line + "\n")
                    currentline=""
                else:
                    line = line.rstrip("\n")
                    currentline = currentline + line
            temp_aln.write(currentline)
                
    collist = [int(i)-1 for i in sys.argv[2:]]
    #print collist
    seqdict = get_seqdict("aln_oneline.tmp")
    for k in seqdict.keys():
        print ">"+k
        print "".join([seqdict[k][i] for i in collist])

    os.remove("aln_oneline.tmp")               
