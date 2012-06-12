#!/bin/bash

# Run this from your work dir, not within any package

# Checkout SusyNtuple
svn co svn+ssh://svn.cern.ch/reps/atlasinst/Institutes/UCIrvine/SUSYAnalysis/SusyNtuple/trunk SusyNtuple

# All this does is checkout the common code and run the install script there
# This will take a while!
svn co svn+ssh://svn.cern.ch/reps/atlasphys/Physics/SUSY/Analyses/Multilepton/trunk/fiveinvfb/MultiLep MultiLep
source MultiLep/installscripts/install_script.sh
