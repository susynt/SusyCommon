#!/usr/bin/env bash


# Commands to run NtMaker on a test xaod file, and log the output
#
# davide.gerbaudo@gmail.com
# 2014-08-16

function prepare_filelist {
    local dest_file="xaod_filelist.txt"
    local input_files=""
    input_files+=" root://eosatlas//eos/atlas/user/j/jpoveda/r5625_test/AOD.01507244._011801.pool.root.1"
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
    local dest_dir="/afs/cern.ch/user/g/gerbaudo/tmp/susynt_xaod_timing/"
    # show only first and last few lines, see http://unix.stackexchange.com/questions/48777/command-to-display-first-few-and-last-few-lines-of-a-file
    local nlines_head=30
    local nlines_tail=60
    if [[ "${tag}" == *dirty* ]]
    then
        echo "dirty status, commit"
        exit 1
    else
        echo "running"
        NtMaker -f ${input_file} -p mc12a --saveContTau --nLepFilter 1 --nLepTauFilter 2 --filterTrig --sys -d 10 \
            2>&1 | \
            tee >(tail -n ${nlines_tail}) >(head -n ${nlines_head}; cat >/dev/null) >/dev/null | \
            tee ${tmp_log}
        mv ${tmp_log} ${dest_dir}
        echo "Logged to ${dest_dir}/${tmp_log}"
    fi
}

main
