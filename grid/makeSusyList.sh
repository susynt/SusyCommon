#!/bin/bash

#
# Generate susy simplified model sample list for grid submission from a list of datasets
# You need to have tier3 ami setup in order to use this
#

if [[ $# < 1 ]]; then
        echo "Supply list of samples"
        exit 1
fi

function getXsec {

        ds=$1
        # extract dataset ID
        id=${ds#mc12_8TeV.}
        id=${id%.Herwigpp*}
        #grep $id $TestArea/SusyCommon/data/modeA_lightslep_MC1eqMN2.txt | awk '{print $5}'
        #cat $TestArea/SusyCommon/data/modeA_lightslep_MC1eqMN2.txt $TestArea/SusyCommon/data/modeC_lightslep_MC1eqMN2.txt | grep $id | awk '{print $5}'
        cat $TestArea/SusyCommon/data/mode*_lightslep_MC1eqMN2.txt | grep $id | awk '{print $5}'

}


input=$1

for ds in `cat $input`; do

        # Get the sum of weights from AMI (ok for herwigpp samples)
        sumw=`ami dataset info $ds | grep totalEvents | awk '{print $2}'`
        # Get the cross section from the text file (for now)
        xsec=`getXsec $ds`

        echo -e "$ds\t$sumw\t$xsec"

done

