import re, argparse

def nhx_to_new(treestring):
    p = re.compile("\[.*?\]")
    return re.sub(p,"",treestring)

parser = argparse.ArgumentParser()
parser.add_argument("intree",help="the input New Hampshire Extended tree")
args = parser.parse_args()

intree = open(args.intree,"r").read()
outtree = nhx_to_new(intree)
print(outtree)
