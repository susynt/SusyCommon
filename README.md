#  SusyCommon
A RootCore package for common UCI SUSY analysis

To set up and area to write the ntuples, follow the instructions in gerbaudo/susynt-write
To set up and area to read the ntuples, follow the instructions in gerbaudo/susynt-read

After compiling, you can run with
```
cd SusyCommon/run
NtMaker -s [mc12_dummy,data_muons,data_egamma] -f xaod_filelist.txt
```

### Overview of the package
For an overview of the classes, see the doxygen documentation at [this
link](http://gerbaudo.github.io/SusyCommon/doxygen-html/)

## Old links and obsolete instructions
For useful information on using these packages, please consult the TWiki:
https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/UCISusyNtuples

Link to code in [SVN browser](https://svnweb.cern.ch/trac/atlasinst/browser/Institutes/UCIrvine/SUSYAnalysis/SusyCommon)

To checkout the package trunk
```
svn co svn+ssh://svn.cern.ch/reps/atlasinst/Institutes/UCIrvine/SUSYAnalysis/SusyCommon/trunk SusyCommon
```

For the currently recommended tag, please see the TWiki
```
svn co svn+ssh://svn.cern.ch/reps/atlasinst/Institutes/UCIrvine/SUSYAnalysis/SusyCommon/tags/SusyCommon-XX-YY-ZZ SusyCommon
```

To make a tag, use the script/makeTag.sh script (only for developers)
```
./scripts/makeTag.sh SusyCommon-XX-YY-ZZ "Some message"
```

To install, check out SusyNtuple, the MultiLep Common code, and run the MultiLep install script 
**NOTE** that this snippet will checkout the trunks of each code.  If you need to use specific tags, use the appropriate svn commands
It is recommended that you do this installation in a clean workarea, to avoid conflicts with other package tags.
```
svn co svn+ssh://svn.cern.ch/reps/atlasinst/Institutes/UCIrvine/SUSYAnalysis/SusyNtuple/trunk SusyNtuple
svn co svn+ssh://svn.cern.ch/reps/atlasphys/Physics/SUSY/Analyses/WeakProduction/MultiLep/trunk MultiLep
source MultiLep/installscripts/install_script.sh
```
Or, for convenience, these commands have been put in one script (kinit first):
`source SusyCommon/scripts/install.sh`

 for examples, here is my analysis code which depends on this package
`https://svnweb.cern.ch/trac/atlasinst/browser/Institutes/UCIrvine/sfarrell/SusyAna`

### Overview of the package


SusyD3PDInterface
  - TSelector for processing Susy D3PDs
  - Provides access to d3pd objects via 'd3pd' container

SusyD3PDAna
  - Inherits from SusyD3PDInterface
  - Holds all the functionality for doing the analysis out of d3pds.
    Object selection, event cleaning, reweighting, etc.

SusyNtMaker
  - Inherits from SusyD3PDAna
  - Applies cleaning cuts and event filtering
  - Writes a SusyNt tree


# Executables
SusyD3PDTest
  - Just dumps D3PD variables, for testing

NtMaker
  - Runs the SusyNtMaker

SusyNtTest
  - Example executable for processing SusyNt



#------------------------------------------------------------------------------------
# Examples on how to run
#------------------------------------------------------------------------------------

##################
# To Make SusyNt #
##################

1.) compile with Rootcore, or in standalone mode
2.) cd run/
3.) Create a list of d3pd files into a text file

command: NtMaker -s SAMPLE -f FILELIST.txt

The '-s' option sets the sample name and is used to control data/MC toggle.
If the sample name contains 'data' (case insensitive), the job will treat the input as data.
This is necessary when running on data so the job will not try to access MC specific branches of the D3PD.

Also, NtMaker -h will show you the current list of options.  

After you run NtMaker, it will spit out a susyNt.root which will contain the tree.

For more advanced grid production of SusyNt, please see the scripts under
[https://github.com/gerbaudo/susynt-submit](https://github.com/gerbaudo/susynt-submit)


