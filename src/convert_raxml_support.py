import sys,re

# script to convert default raxml "branch" support values (in square brackets following branch lengths) into normal ones. Should, but may not, work even if you have identical brlen support combos

def convert(nwk_str):
    p = re.compile(":\d*.\d*\[\d*\]")
    while p.search(nwk_str) is not None:
        m = p.search(nwk_str)
        r = m.group().split("[")[1].rstrip("]")+m.group().split("[")[0]
        nwk_str = nwk_str[:m.start()]+r+nwk_str[m.end():]
    return nwk_str

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "usage: python "+sys.argv[0]+" treefile"
        sys.exit()

    nwk_str = open(sys.argv[1],"r").read()
    #print nwk_str
    print convert(nwk_str)
