#!/bin/bash

tag=n0023
nickname=Steve

# get the samples of interest
if [[ $# = 0 ]]; then
        echo "submit all samples"
        pattern="mc12"
else
        pattern=$1
fi

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
	prun --exec "$command" --tmpDir /tmp --inTarBall=area.tar --useRootCore \
                --excludedSite=ECDF,WEIZMANN,OX,SARA,SHEF,PIC,LPSC,ARC,GLASGOW,GRIF-LAL \
                --extFile '*.so,*.root' --match "*root*" --outputs "susyNt.root" \
                --athenaTag=17.3.1.1 \
		--inDS  $inDS \
		--outDS $outDS
                #--nGBPerJob=MAX \

done
