#!/bin/env python
#
# Count the number of events processed producing SusyNtuple for a given set of samples
# Run without arguments to get the help message.
#
# davide.gerbaudo@gmail.com
# 2013-03

import glob, os, re, sys

import ROOT as r
r.gROOT.SetBatch(1)
r.gErrorIgnoreLevel = 9999 # suppress messages about missing dict
                           # (can't get rid of the 'duplicate' ones?)
treeName = 'susyNt'
binLabel='Initial'
printBinLabels=False       # usually we count the 'Initial' numbers, but there are others available

if len(sys.argv) < 2 :
    print "Usage : %s targetDir [-r regexp-filter]" % sys.argv[0]
    print "Example :"
    print "%s /gdata/atlas/ucintprod/SusyNt/mc12_n0127 -r 'Sherpa_CT10.*(WW|VV)\'" % sys.argv[0]
    sys.exit()
target = sys.argv[1]
regex = sys.argv[3] if len(sys.argv)>3 and sys.argv[2]=='-r' else None

sampleDirs = sorted([t for t in glob.glob(target+'/*')
                     if os.path.isdir(t) and (re.search(regex, t) if regex else True)])

def getProcessedEvents(filename, histoName='rawCutFlow', binLabel=binLabel, printBinLabels=printBinLabels) :
    f = r.TFile.Open(filename)
    histoName = 'genCutFlow'
    histo = f.Get(histoName)
    h = histo
    if printBinLabels : print [h.GetXaxis().GetBinLabel(i) for i in range(1,1+h.GetNbinsX())]
    processedEvents = histo.GetBinContent(histo.GetXaxis().FindBin(binLabel))
    f.Close()
    return processedEvents

print '-'*60
for sd in sampleDirs :
    container = os.path.basename(sd)
    files = glob.glob(sd+'/'+'*.root*')
    nEntries = sum([getProcessedEvents(f) for f in files])
    print container + ' : '+str(nEntries)
    continue
