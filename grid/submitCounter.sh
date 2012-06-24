#!/bin/bash

iteration=i1

# get the samples of interest
if [[ $# = 0 ]]; then
        echo "submit all samples"
        pattern="mc12"
else
        pattern=$1
fi
matches=(`cat newMcSamples.txt | grep $pattern | tr '\t' ','`)
echo "${#matches[@]} matches"

# set it up manually I guess
#source ${ATLAS_LOCAL_ROOT_BASE}/packageSetups/atlasLocalPandaClientSetup.sh

# Loop over samples
for line in ${matches[@]}; do

        info=(`echo $line | tr ',' ' '`)
        #sample=${info[0]}
        inDS=${info[0]}
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
                --excludedSite=RHUL,OX,SARA,SHEF \
		--outputs "sumWeights.root" \
                --athenaTag=17.3.1.1 \
		--inDS  $inDS \
		--outDS $outDS
        

done
