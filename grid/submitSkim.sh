#!/bin/bash

outDS="user.Steve.group.phys-susy.data12_8TeV.eventSkim.NTUP_SUSY.t0pro13_v01_p1181_i0001/"
exclude=RHUL,OX,SARA,SHEF,LRZ,GOEGRID,GRIF,IFIC,CERN,IFAE,INFN

# get the samples of interest
pattern="data"
if [[ $# > 0 ]]; then
        pattern=$1
fi

# Loop over samples
for inDS in `cat dataSamples.txt | grep "$pattern"`; do

        command="./skimEvents.py %IN -l metEvents"

	echo 
	echo "__________________________________________________________________________________________________________"
        echo "INPUT   $inDS"
        echo "OUTPUT  $outDS"
        echo "command $command"
	
	# prun command
	prun --exec "$command" --tmpDir /tmp --noBuild \
                --nGBPerJob=MAX \
                --excludedSite=$exclude \
		--outputs "d3pd.root" \
                --athenaTag=17.0.5.5 \
		--inDS  $inDS \
		--outDS $outDS
        
done
