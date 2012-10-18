#!/bin/env python

import sys
from ROOT import *

##
## eventCounter - a tool for summing MC generator weights
##

def main():

    if len(sys.argv) < 2:
        print 'usage:  eventCounter file1,file2,...'
        return

    print 'eventCounter begin'

    # parse files
    files = sys.argv[1].split(',')
    print 'files:', files

    # d3pd chain
    d3pdChain = TChain('susy')
    # meta data chain
    metaChain = TChain('susyMeta/CutFlowTree')

    for file in files:
        d3pdChain.Add(file)
        metaChain.Add(file)
    d3pdEntries = d3pdChain.GetEntries()
    metaEntries = metaChain.GetEntries()
    metaChain.ls()
    print 'total d3pd entries:', d3pdEntries
    print 'total meta entries:', metaEntries

    # output root file and histos
    f = TFile('sumWeights.root', 'recreate')
    hRaw = TH1F('rawEntries', 'rawEntries', 1, 0, 2)
    h = TH1F('sumw', 'sumw', 1, 0, 2)

    # Take raw entries directly from chain
    hRaw.Fill(1, d3pdEntries)

    # loop over susy metadata tree
    for entry in metaChain:

        names        = [name for name in entry.name]
        descriptions = [desc for desc in entry.description]
        nEvents      = [nEvent for nEvent in entry.nWeightedAcceptedEvents]

        # loop over filters
        for name, desc, nEvent in zip(names, descriptions, nEvents):
            
            print name, nEvent
            if name == 'AllExecutedEvents' and desc == 'Nb of executed events before any cut':
            
                h.Fill(1, nEvent)

                # At least in data, there can be more than one 'AllExecutedEvents', so only save one
                break


    # Print summary
    print 'count results:'
    print '  d3pd', hRaw.GetBinContent(1)
    print '  sumw', h.GetBinContent(1)

    # Write out the histo
    h.Write()
    f.Close()

    print 'eventCounter end'


if __name__ == '__main__':
    main()
