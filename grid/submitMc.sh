#!/bin/bash

# Settings
tag=n0101
nickname=Steve
nGBPerJob=MAX
athenaTag=17.3.1.1
metFlav="Egamma10NoTau_STVF"


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

        # Extract dataset info
        # By default, no cross section is needed and the xsec from SUSYTools will be used
        # However, if the xsec column exists in the input file, that xsec will be used
        inDS=${info[0]}
        sumw=${info[1]}
        xsec=-1
        if [ ${#info[@]} -gt 2 ]; then
                echo "Setting xsec"
                xsec=${info[2]}
        fi

        sample=${inDS#mc12_8TeV.*.}
        sample=${sample%.merge.*/}

        # Build the output dataset name
        outDS="user.$nickname.${inDS%/}_$tag/"
        outDS=${outDS/merge.NTUP_SUSY/SusyNt}

        # AF2 flag
        if [[ $inDS = *_a???_* ]]; then
                af2Opt="--af2"
        fi

        # If output DS name is too long, need to trim it down.
        # For now, treat this on a case by case basis
        if [[ ${#outDS} > 131 ]]; then
                # For Powheg WZ samples
                outDS=${outDS/2LeptonFilter/2L}
        fi

        command="./gridScript.sh %IN --saveTau -s $sample -w $sumw -x $xsec --metFlav $metFlav $af2Opt"

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
