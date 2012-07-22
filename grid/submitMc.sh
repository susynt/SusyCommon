#!/bin/bash

# Settings
tag=n0036
nickname=Steve
nGBPerJob=MAX
#nGBPerJob=14
athenaTag=17.3.1.1


# get the samples of interest
if [[ $# = 0 ]]; then
        echo "submit all samples"
        pattern="mc12"
else
        pattern=$1
fi

# User blacklist of sites
blackList=`cat blacklist.txt`

# A trick to parse the text, first separate columns by commas
matches=(`cat mcSamples.txt | grep $pattern | tr '\t' ','`)
echo "${#matches[@]} matches"


# Loop over samples
for line in ${matches[@]}; do

        # Now replace the commas with normal whitespace
        info=(`echo $line | tr ',' ' '`)
        inDS=${info[0]}
        sumw=${info[1]}
        sample=${inDS#mc12_8TeV.*.}
        sample=${sample%.merge.*/}

        # Build the output dataset name
        outDS="user.$nickname.${inDS%/}_$tag/"
        outDS=${outDS/merge.NTUP_SUSY/SusyNt}

        command="./gridScript.sh %IN -s $sample -w $sumw"

	echo 
	echo "__________________________________________________________________________________________________________"
        echo "INPUT    $inDS"
        echo "OUTPUT   $outDS"
        echo "sample   $sample"
        echo "sumw     $sumw"
        echo "command  $command"

	
	# prun command
	prun --exec "$command" --useRootCore --tmpDir /tmp --inTarBall=area.tar \
             --extFile '*.so,*.root' --match "*root*" --outputs "susyNt.root" \
             --excludedSite=$blackList \
             --nGBPerJob=$nGBPerJob \
             --athenaTag=$athenaTag \
             --inDS  $inDS \
             --outDS $outDS

done
