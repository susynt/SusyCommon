#!/bin/bash

#
# Extract cross sections from AMI
# Using tier3 software interface, setup AMI before using
#

# Pass input file with list of datasets
if [[ $# < 1 ]]; then
        echo "Usage:"
        echo " > makeXsec.sh <DATASET-FILE>"
        exit 1
fi

for ds in `cat $1 | awk '{print $1}'`; do

        eff=`ami dataset evtinfo $ds | grep GenFiltEff_mean | awk '{print $2}'`
        xsec=`ami dataset evtinfo $ds | grep crossSection_mean | awk '{print $2*1000.}'`
        sumw=`ami dataset info $ds | grep totalEvents | awk '{print $2}'`

        #echo $ds
        #echo "xsec $xsec"
        #echo "eff  $eff"
        #effXsec=`echo "$xsec*$eff" | bc`
        effXsec=`echo "$xsec $eff" | awk '{print $1*$2}'`
        #echo "xsec*eff $effXsec"

        echo -e "$ds\t$sumw\t$effXsec"

done
