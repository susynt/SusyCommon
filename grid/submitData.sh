#!/bin/bash

iteration="i15"
#iteration="test4"

# get the samples of interest
if [[ $# = 0 ]]; then
        echo "submit all samples"
        pattern="data"
else
        pattern=$1
fi
datasets=(`cat dataSamples.txt | grep $pattern`)
echo "${#datasets[@]} datasets"


# Setup Panda before running this script
#source ${ATLAS_LOCAL_ROOT_BASE}/packageSetups/atlasLocalPandaClientSetup.sh

# Loop over datasets
for inDS in ${datasets[@]}; do

	# get sample name from input dataset name
        sample=${inDS#*phys-susy.}
        sample=${sample%.PhysCont.*}
        sample=${sample/physics_/}

	# final output ds name
	outDS="user.Steve.$iteration.$sample.SusyNt/"

        command="./gridScript.sh %IN -s $sample"

	echo 
        echo "__________________________________________________________________________________________________________"
	echo "INPUT   $inDS"
	echo "OUTPUT  $outDS"
        echo "sample  $sample"
	
	# prun command
	prun --exec "$command" --useRootCore --tmpDir /tmp --inTarBall=area.tar \
             --excludedSite=OX,SARA,SHEF,PIC,FZK,LPSC,ARC,GLASGOW,GRIF-LAL \
             --extFile '*.so,*.root' --match "*root*" --outputs "susyNt.root" \
             --nGBPerJob=MAX \
             --athenaTag=17.0.5.5 \
	     --inDS  $inDS \
	     --outDS $outDS

done
