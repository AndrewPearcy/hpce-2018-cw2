#!/bin/bash
# The above shows that it is a bash file

# Create a local variable with a list of K values
KS="2 4 8 16 32 64"
# Iterate over them and print (echo) them
for K in $KS; do
    # Run the program with the chosen K, and save to dump_K.csv
    bin/time_fourier_transform hpce.abp14.fast_fourier_transform_taskgroup 0 10 "${K}, " > dump_${K}.csv
done

