#!/bin/bash

#outDS="user.Steve.group.phys-susy.data12_8TeV.eventSkim.NTUP_SUSY.t0pro13_v01_p1181_i0001/"
outDS="user.sfarrell.group.phys-susy.data12_8TeV.periodD.physics_Muons.eventSkim.NTUP_SUSY.t0pro13_v01_p1181_i0002/"
exclude=MWT2,SARA,GRIF,MANC,TRIUMF,AGLT2

# get the samples of interest
pattern="data"
if [[ $# > 0 ]]; then
        pattern=$1
fi

# Loop over samples
#for inDS in `cat dataSamples.txt | grep "$pattern"`; do
for inDS in `cat testSamples.txt | grep "$pattern"`; do

        command="./skimEvents.py %IN -l testEvents"

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
                --athenaTag=17.3.1.1 \
		--inDS  $inDS \
		--outDS $outDS
        
done
