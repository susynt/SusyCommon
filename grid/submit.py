#!/bin/env python

#
# SusyNt grid job submit script
# Run with -h to see the options
# You must set up panda tools before running this script
# 
# Examples
#
#  To process all dgemt:
#  > ./submit.py susy -f dgemt.txt
#
#  To process period B of data and assign tag n9999:
#  > ./submit.py data -p periodB -t n9999
#


from argparse import ArgumentParser
import re
import subprocess

# Some grid option defaults
defaultTag='n0146'
defaultNickname='sfarrell'
defaultMet='Default'

def get_mc_prod(dataset):
    """
    Determine if a D3PD sample is MC12b by parsing the DS name to extract
    the reconstruction tag (e.g. r4485).
    """
    # Need to parse out the reconstruction tag
    recoTag = re.sub('.*[0-9].._r', '', dataset)
    recoTag = re.sub('_.*', '', recoTag)
    try:
        recoTagNum = int(recoTag)
    except ValueError as e:
        print 'ERROR - unable to parse dataset name for reconstruction tag'
        print 'Dataset:', dataset
        print 'Extracted tag:', recoTag
        raise e
    # This could be extended with other known prods
    mc12bTag = 4485
    if recoTagNum < mc12bTag:
        return 'mc12a'
    else:
        return 'mc12b'

def main():

    # Job arguments
    parser = ArgumentParser(description='SusyCommon grid submission')
    add_arg = parser.add_argument
    add_arg('job', choices=['data', 'mc', 'susy'], 
            help='specifies some default settings, like input file')
    add_arg('-f', '--input-files', nargs='*', 
            help='input file with datasets, can specify more than one')
    add_arg('-p', '--pattern', help='grep pattern to select datasets')
    add_arg('-t', '--tag', default=defaultTag, help='SusyNt tag to assign')
    add_arg('-v', '--verbose', action='store_true', help='verbose output')
    add_arg('--athenaTag', default='17.2.9.1', help='athena tag')
    add_arg('--nickname', default=defaultNickname, help='grid nickname, for naming output DS')
    add_arg('--destSE', default='SLACXRD_SCRATCHDISK', 
            help='replicate output dataset to specified site')
    add_arg('--met', default=defaultMet, help='MET flavor to use', 
            choices=['STVF', 'STVF_JVF', 'Default'])
    add_arg('--doMetFix', action='store_true', help='Turns on MET ele-jet overlap fix')
    add_arg('--contTau', action='store_true', help='Store container taus')
    add_arg('--nLepFilter', default='1', help='Number of preselected light leptons to filter on')
    add_arg('--nLepTauFilter', default='2', help='Number of preselected light+tau to filter on')
    add_arg('--filterTrig', action='store_true', help='Turn on trigger filter')
    add_arg('--sys', type=bool, default=True, help='toggle systematics')
    add_arg('--nFilesPerJob', default=None, help='prun option')
    add_arg('--nGBPerJob', default='MAX', help='prun option')
    add_arg('--noSubmit', action='store_true', help='test prun without submitting')
    add_arg('--useShortLivedReplicas', action='store_true', help='prun option')
    add_arg('--cmtConfig', default=None, help='prun option to set cmt config')
    add_arg('--saveTruth', action='store_true', help='Store truth info')
    add_arg('--filterOff', action='store_true', help='Disable event filters (GRL...TileTrip)')
    args = parser.parse_args()

    # Standard options for data
    if args.job == 'data':
        #input_files = ['dataSamples.txt']
        input_files = ['data12_Egamma.txt', 'data12_Muons.txt']
        pattern = 'data'

    # Standard options for standard model mc
    elif args.job == 'mc':
        input_files = ['mcSamples.txt']
        pattern = 'mc'

    # Standard options for susy signals
    else:
        input_files = ['susySamples.txt']
        pattern = 'mc'

    # Override standards
    if args.input_files: input_files = args.input_files
    if args.pattern: pattern = args.pattern

    # Blacklisted sites
    with open('./blacklist.txt') as f:
        blacklist = f.read()
        blacklist = blacklist.replace('\n', '')

    # Print job
    print 'Submitting', args.job, args.tag
    print 'input file:', input_files
    print 'pattern:   ', pattern

    # Loop over inputs
    for input_file in input_files:
        with open(input_file) as inputs:
            for line in inputs:

                info = line.split()
                if len(info) == 0: continue

                # Match pattern
                if re.search(pattern, line) == None: continue

                # Extract sumw and xsec, if provided
                inDS = info[0]
                sumw, xsec, errXsec = '1', '-1', '-1'
                if len(info) > 1: sumw = info[1]
                if len(info) > 2: xsec = str(eval(info[2])) if '*' in info[2] else info[2]
                if len(info) > 3: errXsec = info[3]

                # Get sample name
                sample = re.sub('.merge.*', '', inDS)
                sample = re.sub('mc12_8TeV\.[0-9]*\.', '', sample)
                sample = re.sub('.*phys-susy\.', '', sample)
                sample = re.sub('\.PhysCont.*', '', sample)
                sample = re.sub('physics_', '', sample)

                # Output dataset
                outDS = 'user.'+args.nickname+'.'+re.sub('/', '', inDS)+'_'+args.tag+'/'
                outDS = re.sub('NTUP_SUSY', 'SusyNt', outDS)
                outDS = re.sub('SKIM', '', outDS)
                outDS = re.sub('merge\.', '', outDS)

                # Filter output dataset names that are too long
                # Do this on a case by case basis
                if len(outDS) > 131:
                    outDS = re.sub('2LeptonFilter', '2L', outDS)
                    outDS = re.sub('UEEE3_CTEQ6L1_', '', outDS)
                    outDS = re.sub('AUET2CTEQ6L1_', '', outDS)
                    outDS = re.sub('AUET3CTEQ6L1_', '', outDS)
                    outDS = re.sub('AUET2BCTEQ6L1_', '', outDS)
                    outDS = re.sub('AU2CT10_', '', outDS)

                # Grid command
                gridCommand = './gridScript.sh %IN --metFlav ' + args.met
                gridCommand += ' --nLepFilter ' + args.nLepFilter
                gridCommand += ' --nLepTauFilter ' + args.nLepTauFilter
                gridCommand += ' -w ' + sumw + ' -x ' + xsec + ' -s ' + sample
                gridCommand += ' --errXsec ' + errXsec

                # MC production tag
                if args.job != 'data':
                    mcProd = get_mc_prod(inDS)
                    gridCommand += ' -p ' + mcProd

                # Container taus - forced on, for now
                #gridCommand += ' --saveContTau' if args.contTau else ' --saveTau'
                gridCommand += ' --saveContTau'

                # Met fix
                if args.doMetFix: gridCommand += ' --doMetFix'

                # Trigger filtering
                if args.filterTrig: gridCommand += ' --filterTrig'

                # Truth leptons
                gridCommand += ' --saveTruth' if args.saveTruth else ''

                # Turn off all filtering
                gridCommand += ' --filterOff' if args.filterOff else ''

                # Systematics
                if args.sys: gridCommand += ' --sys'

                # AF2 sample option
                if re.search('_a[0-9]*_', inDS): gridCommand += ' --af2'

                print '\n' + ('_'*90)
                print 'Input   ', inDS
                print 'Output  ', outDS
                print 'Sample  ', sample
                print 'Command ', gridCommand
                print ''

                # The prun command
                prunCommand = 'prun --exec "' + gridCommand + '" --useRootCore --tmpDir /tmp '
                prunCommand += ' --inDS ' + inDS + ' --outDS ' + outDS
                prunCommand += ' --inTarBall=area.tar --extFile "*.so,*.root" --match "*root*"'
                prunCommand += ' --safetySize=600'
                prunCommand += ' --outputs "susyNt.root,gridFileList.txt"'
                prunCommand += ' --destSE=' + args.destSE
                prunCommand += ' --athenaTag=' + args.athenaTag
                prunCommand += ' --excludedSite=' + blacklist

                # You can only have one of the following options
                if(args.nFilesPerJob is not None):
                    prunCommand += ' --nFilesPerJob=' + args.nFilesPerJob
                else:
                    prunCommand += ' --nGBPerJob=' + args.nGBPerJob

                # For testing
                if(args.noSubmit): prunCommand += ' --noSubmit'
                if(args.useShortLivedReplicas):
                    prunCommand += ' --useShortLivedReplicas'
                if(args.cmtConfig is not None):
                    prunCommand += ' --cmtConfig ' + args.cmtConfig

                # Execute prun command
                if args.verbose: print prunCommand
                subprocess.call(prunCommand, shell=True)

if __name__ == '__main__':
    main()
