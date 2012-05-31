#!/bin/env python

import sys
from ROOT import *

##
## eventCounter - a tool for summing MC generator weights
##

def main():

    if len(sys.argv) < 2:
        print "usage:  eventCounter file1,file2,..."
        return

    print "eventCounter begin"

    # parse files
    files = sys.argv[1].split(',')
    print "files:", files

    # build chain
    chain = TChain("susyMeta/CutFlowTree")
    for file in files:
        chain.Add(file)
    nEntries = chain.GetEntries()
    chain.ls()
    print "total entries:", nEntries

    # output root file and histo
    f = TFile("sumWeights.root", "recreate")
    h = TH1F("sumw", "sumw", 1, 0, 2)

    # loop over susy metadata tree
    for entry in chain:

        names        = [name for name in entry.name]
        descriptions = [desc for desc in entry.description]
        nEvents      = [nEvent for nEvent in entry.nWeightedAcceptedEvents]

        # loop over filters
        for name, desc, nEvent in zip(names, descriptions, nEvents):
            
            if name == "AllExecutedEvents" and desc == "Nb of executed events before any cut":
            
                print name, nEvent
                h.Fill(1, nEvent)


    # Write out the histo
    h.Write()
    f.Close()

    print "eventCounter end"


if __name__ == "__main__":
    main()
