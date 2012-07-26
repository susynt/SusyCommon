#!/bin/bash

# Settings
tag=n0036
nickname=Steve
nGBPerJob=MAX
#nGBPerJob=14
athenaTag=17.3.1.1


# dataset file and grep pattern can be provided as arguments
pattern="mc12"
dsFile="susySamples.txt"
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
        xsec=${info[2]}
        sample=${inDS#mc12_8TeV.*.}
        sample=${sample%.merge.*/}

        # Build the output dataset name
        outDS="user.$nickname.${inDS%/}_$tag/"
        outDS=${outDS/merge.NTUP_SUSY/SusyNt}

        command="./gridScript.sh %IN -s $sample -w $sumw -x $xsec"

	echo 
	echo "__________________________________________________________________________________________________________"
        echo "INPUT    $inDS"
        echo "OUTPUT   $outDS"
        echo "sample   $sample"
        echo "sumw     $sumw"
        echo "xsec     $xsec"
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
