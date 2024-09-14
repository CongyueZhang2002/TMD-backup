#!/bin/bash

> params.dat

for i in {100..300};do
for file in replica_$i/result.txt; do
    replica_number=$(echo "$file" | cut -d'/' -f1 | sed 's/replica_//')

    if [ "$replica_number" == "0" ] || [ "$replica_number" -gt 300 ]; then
        continue
    fi

    awk -v replica="$replica_number" 'NR==8 { print $0, replica }' "$file" >> params.dat
done
done



 

