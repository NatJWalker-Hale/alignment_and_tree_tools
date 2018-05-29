import sys

## Usage: python get_columns_from_stockholm.py alignment.fasta alignment.stockholm cutoff 

## This extracts columns from a fasta format alignment from FSA matching a specified accuracy cutoff in the corresponding stockholm alignment. Currently, cutoff must be an integer 1 <= cutoff <= 9, as
## the stockholm format produced by fsa only gives accuracy to 1 decimal place (e.g. 1 corresponds to column accuracy of 0.1 or less).  

## first, put all sequences on one line

cutoff = int(sys.argv[3])

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

## extract positions corresponding to threshold sites
        
with open(sys.argv[2],"r") as stockholm:
    for line in stockholm:
        if line.startswith("#=GC"):
            line = line.lstrip("#=GC Accuracy").rstrip()
            sites = [site for site, value in enumerate(line) if int(value) >= cutoff]

with open("aln_oneline.tmp","r") as alignment:
    for line in alignment:
        if line.startswith(">"):
            line = line.rstrip()
            print(line)
        else:
            for index,site in enumerate(sites):
                if index == len(sites) - 1:
                    sys.stdout.write(line[site] + "\n")
                else:
                    sys.stdout.write(line[site])
                    
            
