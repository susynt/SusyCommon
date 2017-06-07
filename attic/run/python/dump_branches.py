#!/bin/env python

# Script used to explore the available branched in the different NTUP_* formats
#
# davide.gerbaudo@gmail.com
# Jun 2014

import sys
import ROOT as r
r.gROOT.SetBatch(1)

filename = sys.argv[1]
treename = 'physics' if len(sys.argv)<2 else sys.argv[2] # alternatively, 'susy' for NTUP_SUSY

inputFile = r.TFile.Open(filename)
inputTree = inputFile.Get(treename)

branchNames = sorted([b.GetName() for b in inputTree.GetListOfBranches()])
print '\n'.join(branchNames)

