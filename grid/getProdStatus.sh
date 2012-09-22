#!/bin/bash

# This script uses the Tier3 AMI software interface to print the prodSysStatus of a list of datasets in an input file
# Setup AMI and voms proxy prior to running this script

# Pass input file with list of datasets
if [[ $# < 1 ]]; then
        echo "Usage:"
        echo " > getProdStatus.sh <DATASET-FILE>"
        exit 1
fi

for ds in `cat $1`; do

        echo $ds
        prodStatus=`ami dataset info $ds | grep prodsysStatus | awk '{print $2,$3,$4,$5,$6,$7}'`
        echo -e "        $prodStatus"

done
