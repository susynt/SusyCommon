#!/bin/bash

if [[ $# > 0 ]]; then
        input=$1
else
        input=mcSamples.txt
fi

for line in `cat $input | tr '\t' ','`; do

        info=(`echo $line | tr ',' ' '`)

        #sample=${info[0]}
        ds=${info[0]}
        sumw=`ami dataset info $ds | grep totalEvents | awk '{print $2}'`
        echo -e "$ds\t$sumw"

done

