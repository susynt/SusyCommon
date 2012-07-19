#!/bin/bash

iteration=i3

#
# Submit counter job for all mc datasets in input file passed as argument
#

if [[ $# < 1 ]]; then
        echo "Pass DS file as input argument"
        exit 1
fi

# set it up manually I guess
#source ${ATLAS_LOCAL_ROOT_BASE}/packageSetups/atlasLocalPandaClientSetup.sh

# Loop over samples
for inDS in `cat $1`; do

        sample=${inDS#mc12_8TeV.*.}
        sample=${sample%.merge.*/}

	outDS="user.Steve.$iteration.$sample.eventCounter/"

        command="./eventCounter.py %IN"

	echo 
	echo "__________________________________________________________________________________________________________"
        echo "INPUT   $inDS"
        echo "OUTPUT  $outDS"
        echo "sample  $sample"
        echo "command:"
        echo "    $command"
	
	# prun command
	prun --exec "$command" --tmpDir /tmp --noBuild \
                --excludedSite=RHUL,OX,SARA,SHEF,LRZ \
		--outputs "sumWeights.root" \
                --athenaTag=17.3.1.1 \
		--inDS  $inDS \
		--outDS $outDS
        

done
