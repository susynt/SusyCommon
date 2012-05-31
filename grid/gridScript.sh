#!/bin/bash

echo "**** GRID SCRIPT BEGIN ****"

# Convert the comma separated list of files into an array
list=(`echo $1 | tr ',' ' '`)
nFiles=${#list[@]}
echo $nFiles input files
shift

# Write the list of files to the fileList
if [ -f gridFileList.txt ]; then
    rm gridFileList.txt
fi

for file in ${list[@]}; do
    echo $file >> gridFileList.txt
done

# Remaining options are passed directly to NtMaker
echo "Passing options to job: $@"

# Now, run the executable
NtMaker -f gridFileList.txt $@

exitcode=$?
if [[ $exitcode != 0 ]]; then
    echo "**** GRID SCRIPT FAIL ****"
    exit 1
fi

echo "**** GRID SCRIPT END ****"
