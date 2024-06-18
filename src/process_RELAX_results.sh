#! /usr/bin/bash

# run this script in a working dir with multiple RELAX jsons. Requires jq.
# writes RELAX output from alternative or partitioned descriptive model to tsv.
# usage: bash process_RELAX_output.sh [MODEL]
#        MODEL: Ralt (RELAX alternative) or PD (partitioned descriptive)

Help()
{
    echo "Run this script in a working dir with multiple RELAX jsons. Requires jq."
    echo "writes RELAX output from alternative or partitioned descriptive model to tsv"
    echo "Usage: bash process_RELAX_output.sh [MODEL]"
    echo "  MODEL: Ralt (RELAX alternative) or PD (partitioned descriptive)"
}

while :; do
    case $1 in
        -h|-\?|--help)
            Help
            exit
    esac
    break
done

echo $(ls *.RELAX.json) 1>&2

echo "Processing above files" 1>&2

if [[ "$1" == "Ralt" ]]
then
    echo -e "gene\tpartition\tlnLnull\tlnLalt\tlrt\tpval\tk\taicc\tw0\tw1\tw2\tp0\tp1\tp2"
    for i in *.RELAX.json
    do
    gene=$(echo $i | cut -f1 -d".")
    lnLalt=$(cat $i | jq '.["fits"]["RELAX alternative"]["Log Likelihood"]')
    lnLnull=$(cat $i | jq '.["fits"]["RELAX null"]["Log Likelihood"]')
    lrt=$(cat $i | jq '.["test results"]["LRT"]')
    pval=$(cat $i | jq '.["test results"]["p-value"]')
    k=$(cat $i | jq '.["test results"]["relaxation or intensification parameter"]')
    aicc=$(cat $i | jq '.["fits"]["RELAX alternative"]["AIC-c"]')
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
    echo -e "$gene\treference\t$lnLnull\t$lnLalt\t$lrt\t$pval\t$k\t$aicc\t$rw0\t$rw1\t$rw2\t$rp0\t$rp1\t$rp2"
    echo -e "$gene\ttest\t$lnLnull\t$lnLalt\t$lrt\t$pval\t$k\t$aicc\t$tw0\t$tw1\t$tw2\t$tp0\t$tp1\t$tp2"
    done
elif [[ "$1" == "PD" ]]
then
    echo -e "gene\tpartition\tlnL\taicc\tw0\tw1\tw2\tp0\tp1\tp2"
    for i in *.RELAX.json
    do
    gene=$(echo $i | cut -f1 -d".")
    lnL=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Log Likelihood"]')
    aicc=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["AIC-c"]')
    rw0=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Reference"]["0"]["omega"]')
    rw1=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Reference"]["1"]["omega"]')
    rw2=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Reference"]["2"]["omega"]')
    tw0=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Test"]["0"]["omega"]')
    tw1=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Test"]["1"]["omega"]')
    tw2=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Test"]["2"]["omega"]')
    rp0=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Reference"]["0"]["proportion"]')
    rp1=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Reference"]["1"]["proportion"]')
    rp2=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Reference"]["2"]["proportion"]')
    tp0=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Test"]["0"]["proportion"]')
    tp1=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Test"]["1"]["proportion"]')
    tp2=$(cat $i | jq '.["fits"]["RELAX partitioned descriptive"]["Rate Distributions"]["Test"]["2"]["proportion"]')
    echo -e "$gene\treference\t$lnL\t$aicc\t$rw0\t$rw1\t$rw2\t$rp0\t$rp1\t$rp2"
    echo -e "$gene\ttest\t$lnL\t$aicc\t$tw0\t$tw1\t$tw2\t$tp0\t$tp1\t$tp2"
    done
fi