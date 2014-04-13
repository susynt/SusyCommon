#!/bin/bash

# Run this from your work dir, not within any package.
# This script will check out the corresponding tags of MultiLep, SUSYTools, etc.
# Checkout SusyNtuple tag separately

# All this does is checkout the common code and run the install script there
# This will take a while!
svn co svn+ssh://svn.cern.ch/reps/atlasphys/Physics/SUSY/Analyses/WeakProduction/MultiLep/tags/MultiLep-01-06-01-02 MultiLep
source MultiLep/installscripts/install_script.sh

