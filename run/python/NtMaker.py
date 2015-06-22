#!/bin/env python

# A configurable python version of util/NtMaker.cxx
#
# davide.gerbaudo@gmail.com

import sys
import os
import optparse

sys.path.append(os.environ['ROOTCOREBIN']+'/../SusyNtuple/')
from susyntuple import utils
r = utils.import_root()

def main():
    print 'load_packages'
    utils.load_packages()
    r.xAOD.Init()
    print 'generate_dicts'
    utils.generate_dicts()
    print 'import_SUSYDefs_enums'
    utils.import_SUSYDefs_enums()
    print 'foo'
    usage="""
%prog foo baz
todo fill info here
"""


    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-n', '--num-events', default=-1, help='default: process all events')
    parser.add_option('-k', '--skip', default=0, help='default: do not skip events')
    parser.add_option('-d', '--debug', default=0, help='debug level, 0=none, ??=??')
    parser.add_option('-f', '--filelist', default='fileList.txt', help='input filelist, default fileList.txt')
    parser.add_option('-s', '--sample', default='mc12_dummy', help='sample name, used to se isMC flag') # todo: deduce from file
    parser.add_option('-p', '--production', help='mc12a or mc12b required for simulated samples') # todo obsolete
    parser.add_option('-w', '--sumw', default=1.0, help='sum of mc weights for norm') # obsolete
    parser.add_option('-x', '--xsec', default=-1.0, help='cross section, default use susy db') # obsolete
    parser.add_option('--errXsec', default=-1, help='cross section uncert')
    parser.add_option('-l', '--lumi', default=5312, help='lumi in 1/pb')
    parser.add_option('-m', '--write-tuple', )
    parser.add_option('--grl', help='good run list')
    parser.add_option('--sys', action='store_true', help='turn on systematics, default off')
    parser.add_option('--save-photons', action='store_true', help='save photons')
    parser.add_option('--save-taus', action='store_true', help='save taus')
    parser.add_option('--save-container-taus', action='store_true', help='save container taus')
    parser.add_option('--save-truth', action='store_true', help='save truth')
    parser.add_option('--af2', action='store_true', help='AtlasFast2 simulated sample')
    parser.add_option('--met-flavor', default='Default', help='met flavor, see...')
    

    parser.add_option('--met-muons', action='store_true', help='use met muons to compute met')
    parser.add_option('--met-fix', action='store_true', help='turn on ele-jet overlap fix')
    parser.add_option('--filter-off', action='store_true', help='do not filter events')
    parser.add_option('--num-lep-filter', default=0, help='minimum number of light leptons required')
    parser.add_option('--num-lep-tau-filter', default=2, help='minimum number of light+tau leptons required')
    parser.add_option('--filter-trig', action='store_true', help='filter on trigger, usually used only for data') # todo guess
    parser.add_option('--triggers', default=2, help='triggers to be used 0=run1, 1=run2 (default)') # todo update, use string


    (opts, args) = parser.parse_args()
    debug = int(opts.debug)
    if debug:
        print "working from {0}".format(os.getcwd())
        print "being called as : {0}".format(' '.join(os.sys.argv))
        allOptions = [x.dest for x in parser._get_all_options()[1:]]
        print "options parsed:\n"+'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions)



    fileName=opts.filelist # Set up the input files
    treeName = "CollectionTree" # default when making transient tree anyway
 
    f = r.TFile.Open(fileName)
    t = r.xAOD.MakeTransientTree( f, treeName) # Make the "transient tree"; is this needed to load the dictionaries for collections?
 
    for i, event in enumerate(t):
        print "event ",i
        if i>10:break

    chain =  r.TChain("CollectionTree");
    fileErr = r.ChainHelper.addInput(chain, opts.filelist, debug);
    if fileErr: return 1
    event = r.xAOD.TEvent (r.xAOD.TEvent.kClassAccess)
    chain.ls()

    susyAna = r.Susy.SusyNtMaker()
    susyAna.setDebug(debug)
    susyAna.setTriggerSet(opts.triggers)
    susyAna.setSample(opts.sample if opts.sample else 'unknown')
    susyAna.setLumi(opts.lumi)
    susyAna.setSumw(opts.sumw)
    if opts.sys : susyAna.setSys(opts.sys)    
    if opts.save_photons : susyAna.setSelectPhotons(opts.save_photons)
    if opts.save_taus : susyAna.setSelectTaus(opts.save_taus)
    if opts.save_container_taus : susyAna.setSaveContTaus(opts.save_container_taus)
    if opts.af2 : susyAna.setAF2(opts.af2)
    susyAna.setXsec(opts.xsec)
    susyAna.setErrXsec(opts.errXsec)
    # susyAna.setMetFlavor(opts.met_flavor)
    if opts.save_truth : susyAna.setSelectTruthObjects(opts.save_truth)
    # //susyAna.setDoMetFix(doMetFix)
    if opts.filter_off : susyAna.setFilter(opts.filter_off) # todo: negative yes or no?
    susyAna.setNLepFilter(opts.num_lep_filter)
    susyAna.setNLepTauFilter(opts.num_lep_tau_filter)
    if opts.filter_trig : susyAna.setFilterTrigger(opts.filter_trig)

    # susyAna.setGRLFile(opts.grl)

    # MCProduction mcProd = MCProd_Unknown
    # if(mcProdStr.CompareTo("mc12a", TString::kIgnoreCase)) mcProd = MCProd_MC12a
    # else if(mcProdStr.CompareTo("mc12b", TString::kIgnoreCase)) mcProd = MCProd_MC12b
    # susyAna.setMCProduction(mcProd)

    nEntries = chain.GetEntries()
    nEvt = nEntries if opts.num_events<0 else opts.num_events
    print "Total entries:   ",nEntries
    print "Process entries: ",nEvt
    chain.Process(susyAna, opts.sample, nEvt, opts.skip)

    print "SusyNtMaker job done"

    sys.exit(0)
    return


if __name__=='__main__':
    main()
