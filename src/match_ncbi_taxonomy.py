import sys,os,subprocess
import xml.etree.ElementTree as ET

def match_ncbi_taxonomy(species_name):
    genus = species_name.split("_")[0].strip()
    species = species_name.split("_")[1].strip()
    taxon = " ".join( [genus,species] )
    search = " ".join( [ 'esearch -db taxonomy -query "',genus,species,'"'] )
    searchoutput = subprocess.Popen(search, shell=True, stdout=subprocess.PIPE).communicate()[0]
    fetch = "efetch -format native -mode xml"
    fetchp = subprocess.Popen(fetch, shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    fetchoutput = fetchp.communicate(input=searchoutput)[0]
    ncbi_taxon = ET.fromstring(fetchoutput).find('Taxon').find('ScientificName').text
    ncbi_id = ET.fromstring(fetchoutput).find('Taxon').find('TaxId').text
    if ncbi_taxon == taxon:
        istaxa = True
    else:
        istaxa = False
    return taxon,ncbi_taxon,ncbi_id,istaxa
    
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "You must specify a file of genus and species epithets, underscore delimited, each on a new line."
        sys.exit(0)

    with open(sys.argv[1],"r") as infile:
        for line in infile:
            try:
                taxon,ncbi_taxon,ncbi_id,istaxa = match_ncbi_taxonomy(line.strip())
                if istaxa:
                    print taxon+" matches the NCBI taxonomy. TaxID: "+ncbi_id
                else:
                    print taxon+" does not match the NCBI taxonomy. The NCBI taxonomy is: "+ncbi_taxon+" TaxID: "+ncbi_id
            except:
                print line.strip()+" is not in the NCBI taxonomy. Check orthography."
                continue
