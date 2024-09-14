#!/bin/bash

filename="params.dat"
cp "$filename" "$filename.bak"

awk -F'\t' -v OFS='\t' 'NR==1 {print; next} { $13 = NR-2; print }' $filename.bak > $filename.bak

#rm "$filename.bak"
#rm "$filename.tmp"
#rm "$filename.tmp1"

#cp "$filename" "$filename.bak"

#awk '{$10="0.0"; $11="0.0"; $12="0.0"; print}' "$filename.bak" > "$filename.tmp"
#awk '{temp = $10; $10 = $13; $13 = temp; print}' "$filename.tmp" > "$filename.tmp1"

#mv "$filename.tmp" "$filename"

#rm "$filename.bak"
#rm "$filename.tmp"
#rm "$filename.tmp1"
