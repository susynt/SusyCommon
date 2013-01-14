#!/bin/bash

# This script takes a file input with a list of datasets, one per line.
# It will spit out DATASET SUMW

# Setup pyAMI before using this script.  You might also need a voms proxy.
# On lxplus, do
# source /afs/cern.ch/atlas/software/tools/pyAMI/setup.sh

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

