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
        # slep and noslep files have different format
        if [[ $ds = *_noslep_* ]]; then
                cat $TestArea/SusyCommon/data/mode*W_MC1eqMN2.txt | grep $id | awk '{print $4}' | tr '*' ' ' | awk '{print $1*$2*$3}'
        else
                cat $TestArea/SusyCommon/data/mode*lightslep_MC1eqMN2.txt | grep $id | awk '{print $5}'
        fi

}
function getXsecUnc {

        ds=$1
        # extract dataset ID
        id=${ds#mc12_8TeV.}
        id=${id%.Herwigpp*}
        # slep and noslep files have different format
        if [[ $ds = *_noslep_* ]]; then
                cat $TestArea/SusyCommon/data/mode*W_MC1eqMN2.txt | grep $id | awk '{print $5}' | awk '{print 0.01*$1}'
        else
                cat $TestArea/SusyCommon/data/mode*lightslep_MC1eqMN2.txt | grep $id | awk '{print $6}' | awk '{print 0.01*$1}' 
        fi
}


input=$1

for ds in `cat $input`; do

        # Get the sum of weights from AMI (ok for herwigpp samples)
        sumw=`ami dataset info $ds | grep totalEvents | awk '{print $2}'`
        # Get the cross section from the text file (for now)
        xsec=`getXsec $ds`
        # Get the cross section uncertainty from the text file (for now)
        xsecunc=`getXsecUnc $ds`

        echo -e "$ds\t$sumw\t$xsec\t$xsecunc"

done

