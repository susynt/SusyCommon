#!/bin/bash

# This script will calculate the lumi per sample for all D3PD datasets in a provided sample list

if [ $# -lt 1 ]; then
        echo "Provide input sample list for lumi calculation"
        exit 1
fi

susyXsecFile=$ROOTCOREDIR/data/SUSYTools/susy_crosssections_8TeV.txt
dsFile=$1

# Optional grep pattern
pattern="mc12"
if [ $# -gt 1 ]; then
        pattern=$2
fi

# Function to extract cross section from SUSYTools file
function getSusyXsec {
        # Extract the run number
        mcID=${1#mc12_8TeV.}
        mcID=${mcID%%.*}
        xsec=`cat $susyXsecFile | grep "^$mcID" | awk '{print $3*$4*$5}'`
        echo $xsec
}


# A trick to parse the text, first separate columns by commas
matches=(`cat $dsFile | grep "$pattern" | tr '\t' ',' | tr ' ' ','`)
echo "${#matches[@]} matches"

echo
echo -e "\tsumw     \txsec     \tlumi [1/fb]"

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
                xsec=${info[2]}
        fi

        if [[ $xsec == "-1" ]]; then
                xsec=`getSusyXsec $inDS`
        fi

        sumw=${sumw/e+/*10^}
        lumi=`echo "scale=1; $sumw / ($xsec * 1000.)" | bc`

        echo $inDS
        echo -e "\t$sumw    \t$xsec    \t$lumi"

done
