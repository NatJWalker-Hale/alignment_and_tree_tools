# Note - gets trees in same order as bp from get_gard_bp_loc.sh

grep '"tree":"' $1 | sed 's/"//g;s/tree://g;s/$/;/g;s/^[ \t]*//' 