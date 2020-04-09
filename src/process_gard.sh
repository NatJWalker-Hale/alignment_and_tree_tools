#! /usr/bin/bash

PREF=$1
OG=$2
BP_SCRIPT=$HOME/Scripts/alignment_and_tree_tools/src/get_gard_bp_loc.sh
T_SCRIPT=$HOME/Scripts/alignment_and_tree_tools/src/get_gard_trees.sh

for i in $PREF.{0..2}
do
    bash $BP_SCRIPT $i > "$i".bp
    bash $T_SCRIPT $i | pxrr -g $OG > "$i".t
    paste -d' ' "$i".bp "$i".t > "$i".bp.t 
done