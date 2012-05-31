#!/bin/bash

iteration="i16"
#iteration="test15"

# get the samples of interest
if [[ $# = 0 ]]; then
        echo "submit all samples"
        pattern="mc11"
else
        pattern=$1
fi
matches=(`cat susySamples.txt | grep $pattern | tr '\t' ','`)
echo "${#matches[@]} matches"

# set it up manually I guess
#source ${ATLAS_LOCAL_ROOT_BASE}/packageSetups/atlasLocalPandaClientSetup.sh

# Loop over samples
for line in ${matches[@]}; do
        info=(`echo $line | tr ',' ' '`)
        sample=${info[0]}
        inDS=${info[1]}
        sumw=${info[2]}
        xsec=1
        lumi=1

	outDS="user.Steve.$iteration.$sample.SusyNt/"

        command="./gridScript.sh %IN -s $sample -w $sumw -x $xsec -l $lumi"

	echo 
	echo "__________________________________________________________________________________________________________"
        echo "INPUT   $inDS"
        echo "OUTPUT  $outDS"
        echo "sample  $sample"
        echo "lumi    $lumi"
        echo "sumw    $sumw"
        echo "xsec    $xsec"

	
	# prun command
	prun --exec "$command" --tmpDir /tmp --inTarBall=area.tar --useRootCore \
		--match "*root*" --outputs "susyNt.root" \
                --nGBPerJob=MAX \
                --extFile '*.so,*.root' \
                --athenaTag=17.0.5.5 \
		--inDS  $inDS \
		--outDS $outDS
                #--excludedSite=ECDF,WEIZMANN,OX,SARA,SHEF,PIC,LPSC,ARC,GLASGOW,GRIF-LAL \

done
