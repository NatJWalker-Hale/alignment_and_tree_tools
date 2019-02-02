import sys,os

def drop_third_pos(seq_str):
    seq_list = [x for x in seq_str]
    del seq_list[2::3]
    seq_dropped_str = "".join(seq_list)
    return seq_dropped_str

if __name__ == "__main__": 
    if len(sys.argv) != 2:
        print "usage: python "+sys.argv[0]+" codon_alignment"
        sys.exit(0)

    with open(sys.argv[1],"r") as aln:
        currentline = ""
        for line in aln:
            if line.startswith(">"):
                if currentline != "":
                    print drop_third_pos(currentline)
                print line.rstrip("\n")
                currentline = ""
            else:
                currentline += line.rstrip("\n")
    print drop_third_pos(currentline)
