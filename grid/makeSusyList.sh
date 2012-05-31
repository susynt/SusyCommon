#!/bin/bash

for ds in `dq2-ls "mc11_7TeV.*.simplifiedModel_wA_slep_noWcascade_*.merge.NTUP_SUSY.e1043_a131_s1353_a145_r2993_p935/" | sort`
do

        # Get the grid point number
        id=${ds#*wA_slep_}
        id=${id%.merge*}

        sample="wA_sl_$id"
        #sumw=`ami dataset info $ds | grep totalEvents | awk '{print $2}'`
        sumw=1

        echo -e "$sample\t$ds\t$sumw"

done
