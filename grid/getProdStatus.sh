#!/bin/bash

# This script uses the Tier3 AMI software interface to print the prodSysStatus of a list of datasets in an input file
# Setup AMI and voms proxy prior to running this script

# Pass input file with list of datasets
if [[ $# < 1 ]]; then
        echo "Usage:"
        echo " > getProdStatus.sh <DATASET-FILE> <ONELINE-OPT>"
        exit 1
fi

# To print result all on one line for easy grepping, pass any 2nd argument
doOneLine=false
if [ $# -gt 1 ]; then doOneLine=true; fi


for ds in `cat $1 | awk '{print $1}'`; do

        prodStatus=`ami dataset info $ds | grep prodsysStatus | awk '{print $2,$3,$4,$5,$6,$7}'`
        if $doOneLine; then
                echo -e "$ds     \t$prodStatus"
        else
                echo $ds; echo -e "        $prodStatus"
        fi

done
