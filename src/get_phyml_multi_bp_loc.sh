grep "^<" $1 | cut -f1 -d":" | sed 's/>/,/g;s/-/,/g;s/<//g'
