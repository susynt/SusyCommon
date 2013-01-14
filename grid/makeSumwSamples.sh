#!/bin/bash

if [[ $# > 0 ]]; then
        input=$1
else
        input=mcSamples.txt
fi

# Loop over dataset names in the input file
# The awk command just helps in case it's a multi-column file, only grabbing the first column
for ds in `cat $input | awk '{print $1}'`; do

        sumw=`ami dataset info $ds | grep totalEvents | awk '{print $2}'`
        echo -e "$ds\t$sumw"

done

