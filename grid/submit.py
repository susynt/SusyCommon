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
defaultTag='n0111'
defaultNickname='Steve'
defaultMet='Egamma10NoTau_STVF'

def main():

    # Job arguments
    parser = ArgumentParser(description='SusyCommon grid submission')
    add_arg = parser.add_argument
    add_arg('job', choices=['data', 'mc', 'susy'], help='specifies some default settings, like input file')
    add_arg('-f', '--input-file', help='input file with datasets')
    add_arg('-p', '--pattern', help='grep pattern to select datasets')
    add_arg('-t', '--tag', default=defaultTag, help='SusyNt tag to assign')
    add_arg('--athenaTag', default='17.3.1.1', help='athena tag')
    add_arg('--nickname', default=defaultNickname, help='grid nickname, for naming output DS')
    add_arg('--destSE', default='SLACXRD_SCRATCHDISK', help='replicate output dataset to specified site')
    add_arg('--met', default=defaultMet, help='MET flavor to use', choices=['Egamma10NoTau', 'Egamma10NoTau_STVF'])
    add_arg('--nGBPerJob', default='MAX', help='prun option')
    add_arg('--sys', type=bool, default=True, help='toggle systematics')
    add_arg('--noSubmit', action='store_true', help='test prun without submitting')
    args = parser.parse_args()

    # Standard options for data
    if args.job == 'data':
        input_file = 'dataSamples.txt'
        pattern = 'data'

    # Standard options for standard model mc
    elif args.job == 'mc':
        input_file = 'mcSamples.txt'
        pattern = 'mc'

    # Standard options for susy signals
    else:
        input_file = 'susySamples.txt'
        pattern = 'mc'

    # Override standards
    if args.input_file: input_file = args.input_file
    if args.pattern: pattern = args.pattern

    # Blacklisted sites
    with open('./blacklist.txt') as f:
        blacklist = f.read()
        blacklist = blacklist.replace('\n', '')

    # Print job
    print 'Submitting', args.job, args.tag
    print 'input file:', input_file
    print 'pattern:   ', pattern

    # Loop over inputs
    with open(input_file) as inputs:
        for line in inputs:

            info = line.split()
            if len(info) == 0: continue

            # Match pattern
            if re.search(pattern, line) == None: continue

            inDS = info[0]
            sumw, xsec = '1', '-1'
            if len(info) > 1: sumw = info[1]
            if len(info) > 2: xsec = info[2]

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

            # Grid command
            gridCommand = './gridScript.sh %IN --saveTau --metFlav ' + args.met
            gridCommand += ' -w ' + sumw + ' -x ' + xsec + ' -s ' + sample

            # AF2 sample option
            if re.search('_a[0-9]*_', inDS): gridCommand += ' --af2'

            print '\n______________________________________________________________________________________________'
            print 'Input   ', inDS
            print 'Output  ', outDS
            print 'Sample  ', sample
            print 'Command ', gridCommand
            print ''

            # The prun command
            prunCommand = 'prun --exec "' + gridCommand + '" --useRootCore --tmpDir /tmp --inTarBall=area.tar ' + \
                          '--extFile "*.so,*.root" --match "*root*" --outputs "susyNt.root,gridFileList.txt" ' + \
                          '--destSE=' + args.destSE + ' --nGBPerJob=' + args.nGBPerJob + \
                          ' --athenaTag=' + args.athenaTag + ' --excludedSite=' + blacklist + \
                          ' --inDS ' + inDS + ' --outDS ' + outDS

            # For testing
            if(args.noSubmit): prunCommand += ' --noSubmit'

            # Execute prun command
            print prunCommand
            subprocess.call(prunCommand, shell=True)

if __name__ == '__main__':
    main()
