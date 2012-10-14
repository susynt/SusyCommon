#!/bin/bash

# Settings
tag=n0107
nickname=Steve
nGBPerJob=MAX
athenaTag=17.3.1.1
metFlav="Egamma10NoTau_STVF"
sysOpt="--sys"


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
# Do this by "translating" tabs into commas and spaces into commas
matches=(`cat $dsFile | grep $pattern | tr '\t' ',' | tr ' ' ','`)
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

        # AF2 flag
        if [[ $inDS = *_a???_* ]]; then
                af2Opt="--af2"
        else
                af2Opt=""
        fi

        # Build the output dataset name
        outDS="user.$nickname.${inDS%/}_$tag/"
        outDS=${outDS/merge.NTUP_SUSY/SusyNt}

        # If output DS name is too long, need to trim it down.
        # For now, treat this on a case by case basis
        if [[ ${#outDS} > 131 ]]; then
                # For simplified model samples
                outDS=${outDS/UEEE3_CTEQ6L1_/}
        fi

        command="./gridScript.sh %IN --saveTau -s $sample -w $sumw -x $xsec --metFlav $metFlav $af2Opt $sysOpt"

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
