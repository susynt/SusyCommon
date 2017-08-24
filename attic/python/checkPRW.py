#!/usr/bin/env python

###############################################
#
# checkPRW.py
#
# open up a PRW config file to check if a
# (set of) mcChannels are represented or
# not
#
# daniel.joseph.antrim@cern.ch
# January 2017
#
##############################################
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal cmd-line options
r.gROOT.SetBatch(True)
import sys

from optparse import OptionParser

def look_for_dsid(dsid, prwfile) :

    rfile = r.TFile.Open(prwfile)
    if not rfile.cd("PileupReweighting") :
        print "look_for_dsid    ERROR Could not find directory 'PileupReweighting' in \
                input file %s"%prwfile
        print "look_for_dsid    ERROR  > exiting"
        sys.exit()

    tdir = rfile.Get("PileupReweighting")
    keys = tdir.GetListOfKeys()

    prw_channels = []
    for key in keys :
        if "pileup_chan" in str(key.GetName()) :
            prw_channels.append(str(key.GetName()).split("_")[1].replace("chan",""))

    n_occurrences = 0
    for channel in prw_channels :
        if str(channel) == str(dsid) :
            n_occurrences += 1

    print 50*"-"
    print " Found %d occurences of channel %s"%(n_occurrences, str(dsid))
    print 50*"-"
    

#_________________________________________
if __name__ == "__main__" :

    input_prw_file = None

    parser = OptionParser()
    parser.add_option("-d", "--dsid", default="",help="The DSID (mc channel) you want to look for in the file")
    parser.add_option("-i", "--input", default="",help="Provide a PRW configuration file")
    (options, args) = parser.parse_args()
    dsid_to_look_for = options.dsid
    input_prw_file   = options.input

    print 50*"-"
    print " checkPRW"
    print ""
    print "\t > dsid        : %s"%str(dsid_to_look_for)
    print "\t > prw file    : %s"%str(input_prw_file)
    print ""
    print 50*"-"

    if not input_prw_file :
        print "checkPRW    ERROR You did not provide an input file!"
        print "checkPRW    ERROR  > exiting" 
        sys.exit()

    if dsid_to_look_for != "" :
        look_for_dsid(dsid_to_look_for, input_prw_file)

