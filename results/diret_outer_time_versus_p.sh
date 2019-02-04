#!/bin/bash
# The above shows that it is a bash file

# Create a local variable with a list of K values
PS="1 2 3 4"
# Iterate over them and print (echo) them
for P in $PS; do
    # Run the program with the chosen K, and save to dump_K.csv
    bin/time_fourier_transform hpce.abp14.direct_fourier_transform_parfor_inner ${P} 30 "${P}, " > dump_${P}.csv
done

