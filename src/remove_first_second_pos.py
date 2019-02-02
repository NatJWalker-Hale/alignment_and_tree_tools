import sys,os

def drop_first_second_pos(seq_str):
    seq_list = [x for x in seq_str]
    seq_dropped_str = "".join(seq_list[2::3])
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
                    print drop_first_second_pos(currentline)
                print line.rstrip("\n")
                currentline = ""
            else:
                currentline += line.rstrip("\n")
    print drop_first_second_pos(currentline)
