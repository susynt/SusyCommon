#!/bin/bash

nickname=$USER
tmpDir=/tmp/$USER
mkdir -p $tmpDir

tarFile=area.tar

if [ -e $tarFile ]; then
    echo "Removing old tar file"
    rm $tarFile
fi

# Create a snapshot list of the SVN tags used to make this tarball, useful for when you forget
svn info $ROOTCOREDIR/../SusyNtuple $ROOTCOREDIR/../SusyCommon $ROOTCOREDIR/../MultiLep $ROOTCOREDIR/../SUSYTools > tarBallTags 

# We don't need to actually submit a job, but still need to provide some dummy arguments to prun
prun --exec "echo" --tmpDir=$tmpDir --outTarBall=$tarFile --useRootCore --noSubmit --athenaTag=17.0.5.5 --outDS=user.$nickname.DummyDoesNothing

