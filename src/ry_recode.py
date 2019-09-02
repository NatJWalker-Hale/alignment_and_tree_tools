import sys, re

# def get_seq_dict(inf): # goes from fasta to python dict for arbitrary number of lines in record
#     seq_dict = {} # key is seq label, value is sequence
#     aln = open(inf,"r")
#     currname = None
#     currseq = None
#     for line in aln:
#         line = line.strip()
#         if line.startswith(">"):
#             currname = line.lstrip(">")
#             seq_dict[currname] = None
#             currseq = None
#         elif currseq is not None:
#             seq_dict[currname] += line
#         else:
#             currseq = line
#             seq_dict[currname] = currseq
#     return seq_dict

def get_seq_dict(inf): # goes from fasta to python dict for arbitrary number of lines in record
    seq_dict = {} # key is seq label, value is sequence
    aln = open(inf,"r")
    currname = None
    for line in aln:
        line = line.strip()
        if line.startswith(">"):
            currname = line.lstrip(">")
            seq_dict[currname] = ""
        else:
            seq_dict[currname] += line
    return seq_dict

def ry_recode(seq_dict):
    ry_seq_dict = {} # key is seq label, value is ry-recoded sequence
    for key in seq_dict:
        seqr = re.sub("[AGag]","R",seq_dict[key])
        seqry = re.sub("[CTct]","Y",seqr)
        ry_seq_dict[key] = seqry
    return ry_seq_dict

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python "+sys.argv[0]+" aln")
        sys.exit()

    seq_dict = get_seq_dict(sys.argv[1])
    ry_dict = ry_recode(seq_dict)
    for key in ry_dict:
        print(">"+key)
        print(ry_dict[key])

        
             
                 
