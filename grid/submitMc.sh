#!/bin/bash

# Settings
tag=n0044
nickname=Steve
nGBPerJob=MAX
#nGBPerJob=14
athenaTag=17.3.1.1


# dataset file and grep pattern can be provided as arguments
pattern="mc12"
dsFile="mcSamples.txt"
while [[ $# > 0 ]]; do
        if [[ $1 == "-f" ]]; then
                shift
                dsFile=$1
        else
                pattern=$1
        fi
        shift
done

# User blacklist of sites
blackList=`cat blacklist.txt`

# A trick to parse the text, first separate columns by commas
matches=(`cat $dsFile | grep $pattern | tr '\t' ','`)
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

        # If output DS name is too long, need to trim it down.
        # For now, treat this on a case by case basis
        echo $outDS
        if [[ ${#outDS} > 131 ]]; then
                # For Powheg WZ samples
                outDS=${outDS/2LeptonFilter/2L}
        fi

        command="./gridScript.sh %IN --saveTau -s $sample -w $sumw"

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
