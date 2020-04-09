import sys

with open(sys.argv[-1]) as f:
    currentline = ""
    for line in f:
        if line.startswith('>'):
            line = line.rstrip('\n')
            if currentline != "": print currentline
            print line
            currentline = ""
        else:
            line = line.rstrip('\n')
            currentline = currentline + line
print currentline
f.close()
