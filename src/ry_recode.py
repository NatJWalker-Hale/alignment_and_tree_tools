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

# def get_seq_dict(inf): # goes from fasta to python dict for arbitrary number of lines in record
#     seq_dict = {} # key is seq label, value is sequence
#     aln = open(inf,"r")
#     currname = None
#     for line in aln:
#         line = line.strip()
#         if line.startswith(">"):
#             currname = line.lstrip(">")
#             seq_dict[currname] = ""
#         else:
#             seq_dict[currname] += line
#     return seq_dict

# def fasta_iter(inf):
#     f = open(inf,"r")
#     fiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
#     for header in fiter:
#         header = next(fiter)[1:].strip()
#         seq = "".join(s.strip() for s in next(fiter))
#         yield header, seq

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
        
def ry_recode(string): # takes sequence string and recodes A and G to R and C and T to Y
    seqr = re.sub("[AGag]","R",string)
    seqry = re.sub("[CTct]","Y",seqr)
    return seqry

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python "+sys.argv[0]+" aln")
        sys.exit()

    fiter = parse_fasta(sys.argv[1])
    for key, value in dict([x for x in parse_fasta(sys.argv[1])]).items():
        print(">"+key)
        print(ry_recode(value))

    # seq_dict = get_seq_dict(sys.argv[1])
    # ry_dict = ry_recode(seq_dict)
    # for key in ry_dict:
    #     print(">"+key)
    #     print(ry_dict[key])