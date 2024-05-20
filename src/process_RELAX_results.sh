#! /usr/bin/bash

# run this script in a working dir with multiple RELAX jsons. Requires jq.
# writes RELAX output from alternative or partitioned descriptive model to tsv.
# usage: bash process_RELAX_output.sh [MODEL]
#        MODEL: Ralt (RELAX alternative) or PD (partitioned descriptive)

Help()
{
    echo "Run this script in a working dir with multiple RELAX jsons. Requires jq."
    echo "writes RELAX output from alternative or partitioned descriptive model to tsv"
    echo "Usage: bash process_RELAX_output.sh"
}

while :; do
    case $1 in
        -h|-\?|--help)
            Help
            exit
    esac
    break
done

if compgen -G "$(pwd)/*.RELAX.json" > /dev/null; then
    echo "Found RELAX results in $(pwd)"
else
    echo "Did not find RELAX results in $(pwd)"
    exit
fi


echo $(ls *.RELAX.json) 1>&2
 
echo "Processing above files" 1>&2

echo -e "gene\tmodel\tlnL\taicc\tlrt\tpval\tk\tpartition\tw0\tw1\tw2\tp0\tp1\tp2"
for i in *.RELAX.json
do
gene=$(echo $i | cut -f1 -d".")
lnLalt=$(cat $i | jq '.["fits"]["RELAX alternative"]["Log Likelihood"]')
lnLnull=$(cat $i | jq '.["fits"]["RELAX null"]["Log Likelihood"]')
lnLpd=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Log Likelihood"]')
lrt=$(cat $i | jq '.["test results"]["LRT"]')
pval=$(cat $i | jq '.["test results"]["p-value"]')
k=$(cat $i | jq '.["test results"]["relaxation or intensification parameter"]')
aiccnull=$(cat $i | jq '.["fits"]["RELAX null"]["AIC-c"]')
aiccalt=$(cat $i | jq '.["fits"]["RELAX alternative"]["AIC-c"]')
aiccpd=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["AIC-c"]')
nw0=$(cat $i | jq '.["fits"]["RELAX null"]["Rate Distributions"]["Reference"]["0"]["omega"]')
nw1=$(cat $i | jq '.["fits"]["RELAX null"]["Rate Distributions"]["Reference"]["1"]["omega"]')
nw2=$(cat $i | jq '.["fits"]["RELAX null"]["Rate Distributions"]["Reference"]["2"]["omega"]')
np0=$(cat $i | jq '.["fits"]["RELAX null"]["Rate Distributions"]["Reference"]["0"]["proportion"]')
np1=$(cat $i | jq '.["fits"]["RELAX null"]["Rate Distributions"]["Reference"]["1"]["proportion"]')
np2=$(cat $i | jq '.["fits"]["RELAX null"]["Rate Distributions"]["Reference"]["2"]["proportion"]')
rw0=$(cat $i | jq '.["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["0"]["omega"]')
rw1=$(cat $i | jq '.["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["1"]["omega"]')
rw2=$(cat $i | jq '.["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["2"]["omega"]')
tw0=$(cat $i | jq '.["fits"]["RELAX alternative"]["Rate Distributions"]["Test"]["0"]["omega"]')
tw1=$(cat $i | jq '.["fits"]["RELAX alternative"]["Rate Distributions"]["Test"]["1"]["omega"]')
tw2=$(cat $i | jq '.["fits"]["RELAX alternative"]["Rate Distributions"]["Test"]["2"]["omega"]')
rp0=$(cat $i | jq '.["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["0"]["proportion"]')
rp1=$(cat $i | jq '.["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["1"]["proportion"]')
rp2=$(cat $i | jq '.["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["2"]["proportion"]')
tp0=$(cat $i | jq '.["fits"]["RELAX alternative"]["Rate Distributions"]["Test"]["0"]["proportion"]')
tp1=$(cat $i | jq '.["fits"]["RELAX alternative"]["Rate Distributions"]["Test"]["1"]["proportion"]')
tp2=$(cat $i | jq '.["fits"]["RELAX alternative"]["Rate Distributions"]["Test"]["2"]["proportion"]')
pdrw0=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Reference"]["0"]["omega"]')
pdrw1=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Reference"]["1"]["omega"]')
pdrw2=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Reference"]["2"]["omega"]')
pdtw0=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Test"]["0"]["omega"]')
pdtw1=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Test"]["1"]["omega"]')
pdtw2=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Test"]["2"]["omega"]')
pdrp0=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Reference"]["0"]["proportion"]')
pdrp1=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Reference"]["1"]["proportion"]')
pdrp2=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Reference"]["2"]["proportion"]')
pdtp0=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Test"]["0"]["proportion"]')
pdtp1=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Test"]["1"]["proportion"]')
pdtp2=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Test"]["2"]["proportion"]')
echo -e "$gene\tnull\t$lnLnull\t$aiccnull\t\t\t1\tReference\t$nw0\t$nw1\t$nw2\t$np0\t$np1\t$np2"
echo -e "$gene\talt\t$lnLalt\t$aiccalt\t$lrt\t$pval\t$k\tReference\t$rw0\t$rw1\t$rw2\t$rp0\t$rp1\t$rp2"
echo -e "$gene\talt\t\t\t\t\t\tTest\t$tw0\t$tw1\t$tw2\t$tp0\t$tp1\t$tp2"
echo -e "$gene\tPD\t$lnLpd\t$aiccpd\t\t\t\tReference\t$pdrw0\t$pdrw1\t$pdrw2\t$pdrp0\t$pdrp1\t$pdrp2"
echo -e "$gene\tPD\t\t\t\t\t\tTest\t$pdtw0\t$pdtw1\t$pdtw2\t$pdtp0\t$pdtp1\t$pdtp2"
done
