#!/bin/bash

tmpDir=/tmp/sfarrell
mkdir -p $tmpDir

tarFile=area.tar

if [ -e $tarFile ]; then
    echo "Removing old tar file"
    rm $tarFile
fi

# We don't need to actually submit a job, but still need to provide some dummy arguments to prun
prun --exec "echo" --tmpDir=$tmpDir --outTarBall=$tarFile --useRootCore --noSubmit --athenaTag=17.0.5.5 --outDS=user.Steve.DummyDoesNothing

