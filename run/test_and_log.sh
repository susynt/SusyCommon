


#!/usr/bin/env bash


# Commands to run NtMaker on a test xaod file, and log the output
#
# davide.gerbaudo@gmail.com
# 2014-08-16

function prepare_filelist {
    local dest_file="xaod_filelist.txt"
    local input_files=""
    #See /scratch/gerbaudo/xaod_example_input/README

    #input_files+=" root://eosatlas//eos/atlas/user/j/jpoveda/r5625_test/AOD.01507244._011801.pool.root.1"
    #input_files+="/scratch/gerbaudo/xaod_example_input/mc14_8TeV.117050.PowhegPythia_P2011C_ttbar.merge.AOD.e1727_s1933_s1911_r5591_r5625/AOD.01522855._015079.pool.root.1"

    #input_files+="/scratch/gerbaudo/xaod_example_input/data12_8TeV.00204442.physics_Egamma.merge.AOD.r5723_p1751_tid01534464_00/AOD.01534464._000002.pool.root.1"  #data xAOD
    #input_files+="/scratch/gerbaudo/xaod_example_input/mc14_8TeV.117050.PowhegPythia_P2011C_ttbar.merge.AOD.e1727_s1933_s1911_r5591_r5625_tid01522855_00/AOD.01522855._015010.pool.root.1" #ttbar xAOD
    
    #input_files+="/scratch/gerbaudo/xaod_example_input/AOD.01507244._011801.pool.root.1"  #Ximo original test file 
    #input_files+="/scratch/gerbaudo/xaod_example_input/mc14_13TeV.204422.Herwigpp_UEEE4_CTEQ6L1_Tt_T800_L1.merge.AOD.e3064_s1982_s2008_r5787_r5853/AOD.01604167._000005.pool.root.1" #ttbar Ximo
    
    #DC14
    #input_files+="/scratch/gerbaudo/xaod_example_input/mc14_8TeV.117050.PowhegPythia_P2011C_ttbar.merge.AOD.e1727_s1933_s1911_r5591_r5625/AOD.01522855._015010.pool.root.1"  #ttbar DxAOD
    #input_files+="/scratch/gerbaudo/xaod_example_input/mc14_8TeV.117050.PowhegPythia_P2011C_ttbar.merge.AOD.e1727_s1933_s1911_r5591_r5625/AOD.01522855._015010.pool.root.1"  #ttbar xAOD
    #input_files+="/scratch/gerbaudo/xaod_example_input/data12_8TeV.00204442.physics_Egamma.merge.AOD.r5723_p1751/AOD.01534464._000002.pool.root.1" #data xAOD
    #input_files+="/scratch/gerbaudo/xaod_example_input/mc14_valid.164288.Herwigpp_UEEE3_CTEQ6L1_simplifiedModel_wA_noslep_lepW_15.recon.AOD.e2658_s1967_s2008_r5787/AOD.04307072._000025.pool.root.1" #signal xAOD


    #DC14 19.1.4.8
    # DATA
    #input_files+="/gdata/atlas/dantrim/SusyAna/xaod/samples/data12_8TeV.00202660.physics_Muons.merge.AOD.r5724_p1751_p1814_tid04559347_00/AOD.04559347._000001.pool.root.1"
    #input_files+="/gdata/atlas/dantrim/SusyAna/xaod/samples/data12_8TeV.00203602.physics_Egamma.merge.DAOD_SUSY1.r5723_p1751_p1784/data12_8TeV/DAOD_SUSY1.04269564._000001.pool.root.1"

    #AOD
    #input_files+="/gdata/atlas/dantrim/SusyAna/xaod/samples/mc14_13TeV.204422.Herwigpp_UEEE4_CTEQ6L1_Tt_T800_L1.merge.AOD.e3064_s1982_s2008_r5787_r5853/AOD.01604167._000001.pool.root.1"
    #input_files+="/data7/atlas/dantrim/SusyAna/xaod/samples/aod_p1872/mc14_13TeV.174830.MadGraphPythia_AUET2BCTEQ6L1_ttbarWjExcl.merge.AOD.e3214_s1982_s2008_r5787_r5853_tid01604171_00/AOD.01604171._000045.pool.root.1"
    #DOAD
    #input_files += "/data7/atlas/dantrim/SusyAna/xaod/samples/topApr30/mc14_8TeV.174830.MadGraphPythia_AUET2BCTEQ6L1_ttbarWjExcl.merge.DAOD_SUSY1.e1672_s1933_s1911_r5591_r5625_p1785_tid04554118_00/DAOD_SUSY1.04554118._000009.pool.root.1"
    input_files += "/gdata/atlas/dantrim/SusyAna/xaod/samples/mc14_13TeV.167757.Sherpa_CT10_ZtautauMassiveCBPt0_CVetoBVeto.merge.DAOD_SUSY1.e2798_s1982_s2008_r5787_r5853_p1872/DAOD_SUSY1.05248867._000003.pool.root.1"

    #input_files+="/gdata/atlas/dantrim/SusyAna/xaod/samples/stop2L/mc14_13TeV.110401.PowhegPythia_P2012_ttbar_nonallhad.merge.DAOD_SUSY1.e2928_s1982_s2008_r5787_r5853_p1846/DAOD_SUSY1.04964196._000044.pool.root.1"

    input_files+="/data7/atlas/dantrim/SusyAna/xaod/samples/mc14_13TeV.174830.MadGraphPythia_AUET2BCTEQ6L1_ttbarWjExcl.merge.DAOD_SUSY1.e3214_s1982_s2008_r5787_r5853_p1872_tid05249011_00/DAOD_SUSY1.05249011._000005.pool.root.1"

    
    if [ ! -f ${dest_file} ]
    then
        touch ${dest_file}
        for X in ${input_files}
        do
            echo ${X} >> ${dest_file}
        done
    fi
    echo ${dest_file}
}

function first_last () {
    head -n 10 -- "$1"
    tail -n 10 -- "$1"
}

 #awk '{a[NR]=$0} NR<=n {print} END{for (i=NR-n+1;i<=NR;i++) print a[i]}' n=02
function main {
    local tag=$(git describe --dirty --tags)
    local stamp=$(date +%F)
    local input_file=$(prepare_filelist)
    local tmp_log="${stamp}_${tag}.log"
    local tmp_detailed_log="${stamp}_${tag}_detailed.log"
    local dest_dir="${WORKAREA}/susynt_xaod_timing/"
    # show only first and last few lines, see http://unix.stackexchange.com/questions/48777/command-to-display-first-few-and-last-few-lines-of-a-file
    local nlines_head=30
    local nlines_tail=60
    #if [[ "${tag}" == *dirty* ]]
    #then
    #    echo "dirty status, commit"
        #exit 1
    #else
        echo "running"
#        NtMaker -f ${input_file} -p mc12a --saveContTau --nLepFilter 1 --nLepTauFilter 2 --filterTrig --sys -d 10 \
#            2>&1 | \
#            tee >(tail -n ${nlines_tail}) >(head -n ${nlines_head}; cat >/dev/null) >/dev/null | \
#            tee ${tmp_log}

	#to run on data
    #NtMaker -f ${input_file} -s data_muons  -d 10 -n 100 2>&1 |tee ${tmp_detailed_log}
	#NtMaker -f ${input_file} -s data_egamma -d 10 -n 100 2>&1  |tee ${tmp_detailed_log}
	#to run on MC
        #NtMaker -f ${input_file} -p mc12a --saveContTau --savePh --nLepFilter 1 --nLepTauFilter 2 --filterTrig -d 0 --sys > ${tmp_detailed_log} 2>&1
	#-- sys -- broken
    NtMaker -f ${input_file} -s mc14 --nLepFilter 1 -n 100 -d 10 -n 500 --sys \
            2>&1 |tee ${tmp_detailed_log} 

	more ${tmp_detailed_log} |head -n  ${nlines_head} > ${tmp_log} 
	more ${tmp_detailed_log} |tail -n  ${nlines_tail} >> ${tmp_log} 
	
        mv ${tmp_log} ${dest_dir}
        mv ${tmp_detailed_log} ${dest_dir}
        echo "Logged to ${dest_dir}/${tmp_log}"
    #fi
}

main
