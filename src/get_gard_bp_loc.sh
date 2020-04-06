grep "\[[0-9]*, [0-9]*\]" $1 | sed 's/\[//g;s/\]//g;s/,//g;s/ $//g' > tmp1.txt
grep '"[0-9]*":{' $1 | sort | uniq | sed 's/"//g;s/:{//g;s/^[ \t]*//' > tmp2.txt
paste -d" " tmp1.txt tmp2.txt 
rm tmp1.txt tmp2.txt
