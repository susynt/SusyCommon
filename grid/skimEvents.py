#!/bin/env python

#
# skimEvents - a script to skim d3pd events
#

import sys

# Always run in batch mode
sys.argv.append('-b')

from ROOT import TChain, TFile
from argparse import ArgumentParser

# Command line options
parser = ArgumentParser(description='D3PD event skimmer')
add_arg = parser.add_argument
# Batch mode
add_arg('-b', '--batch', action='store_true')
# Event list
add_arg('-l', '--event-list', help='txt file of run event')
# Input files
add_arg('inputs', help='comma separated list of input files') 
args = parser.parse_args()

print 'skimEvents begin'

# parse files
files = sys.argv[1].split(',')
#print 'files:', files

# build chain
chain = TChain('susy')
for file in files:
    chain.Add(file)
nEntries = chain.GetEntries()
chain.ls()
print 'total entries:', nEntries

# output root file and histo
outFile = TFile('d3pd.root', 'recreate')

# Build the list of selections
runEventList = []

if args.event_list:
    with open(args.event_list) as eventFile:
        for line in eventFile:
            run, event = line.split()
            runEventList.append('(RunNumber=='+run+' && EventNumber=='+event+')')

sel = '||'.join(runEventList)
print 'selection:'
print sel

# copy tree
tree = chain.CopyTree(sel)
tree.Write()

print 'output entries:', tree.GetEntries()
tree.Scan('RunNumber:EventNumber')

# Write out the histo
outFile.Close()

print 'skimEvents end'

