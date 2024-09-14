#!/bin/bash

# Create or overwrite the new file
> chi2.dat

# Loop over all files with a specific pattern (e.g., *.txt)
for file in replica_*/chi2.dat; do
    awk 'NR==2 {print $1}' "$file" >> chi2.dat
    #awk 'NR==2 {print $2}' "$file" >> chi2.dat
    #awk 'NR==2 {print $3}' "$file" >> chi2.dat
    #awk 'NR==2 {print $4}' "$file" >> chi2.dat
    #awk 'NR==2 {print $5}' "$file" >> chi2.dat
    #awk 'NR==2 {print $6}' "$file" >> chi2.dat
    #awk 'NR==2 {print $7}' "$file" >> chi2.dat
    #awk 'NR==2 {print $8}' "$file" >> chi2.dat
    #awk 'NR==2 {print $9}' "$file" >> chi2.dat
    #awk 'NR==2 {print $10}' "$file" >> chi2.dat
    #awk 'NR==2 {print $11}' "$file" >> chi2.dat
    #awk 'NR==2 {print $12}' "$file" >> chi2.dat
    #awk 'NR==2 {print $13}' "$file" >> chi2.dat
    #awk 'NR==2 {print $14}' "$file" >> chi2.dat
done

