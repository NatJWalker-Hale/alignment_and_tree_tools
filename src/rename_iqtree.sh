#! /usr/bin/bash

# positional input is bash rename_iqtree.sh alignment iqtree.tre
# depends on phyx

Help()
{
    echo "Script to reconvert IQ-TREE names to @ codes"
    echo "Usage: rename_iqtree.sh alignment tree"
    echo "Required arguments:"
    echo "alignment: original FASTA alignment used for inference"
    echo "     tree: IQ-TREE output .treefile"
}

if [ $# -eq 0 ] ; then
Help
else
pxlssq -s $1 -i | cut -f1 -d" " > orig_names.txt
sed 's/@/_/g' orig_names.txt > iq_names.txt
pxrlt -t $2 -c iq_names.txt -n orig_names.txt -o "$2".cn
rm orig_names.txt iq_names.txt
fi