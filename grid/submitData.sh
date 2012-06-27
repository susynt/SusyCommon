#!/bin/bash

tag=n0022
nickname=Steve

# get the samples of interest
if [[ $# = 0 ]]; then
        echo "submit all samples"
        pattern="data"
else
        pattern=$1
fi

datasets=(`cat dataSamples.txt | grep $pattern`)
echo "${#datasets[@]} datasets"


# Loop over datasets
for inDS in ${datasets[@]}; do

	# get sample name from input dataset name
        sample=${inDS#*phys-susy.}
        sample=${sample%.PhysCont.*}
        sample=${sample/physics_/}

	# final output ds name
	#outDS="user.Steve.$sample.SusyNt.$iteration/"
        outDS="user.$nickname.${inDS%/}_$tag/"
        outDS=${outDS/NTUP_SUSY/SusyNt}

        command="./gridScript.sh %IN -s $sample"

	echo 
        echo "__________________________________________________________________________________________________________"
	echo "INPUT    $inDS"
	echo "OUTPUT   $outDS"
        echo "sample   $sample"
        echo "command  $command"
	
	# prun command
	prun --exec "$command" --useRootCore --tmpDir /tmp --inTarBall=area.tar \
             --excludedSite=OX,SARA,SHEF,PIC,FZK,LPSC,ARC,GLASGOW,GRIF-LAL \
             --extFile '*.so,*.root' --match "*root*" --outputs "susyNt.root" \
             --athenaTag=17.3.1.1 \
             --nGBPerJob=MAX \
	     --inDS  $inDS \
	     --outDS $outDS

done
