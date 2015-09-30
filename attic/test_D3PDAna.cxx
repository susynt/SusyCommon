

/**
   Test D3PDAna and check whether we can read a few events from NTUP_COMMON

   davide.gerbaudo@gmail.com
   June 2014
 */
#include <cstdlib>
#include <string>

#include "TChain.h"
#include "Cintex/Cintex.h"
#include "TSystem.h"

#include "SusyCommon/D3PDAna.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/ChainHelper.h"

using namespace std;
using namespace Susy;

//----------------------------------------------------------
void printHelp(const char *exeName)
{
  cout<<"Usage :"<<endl
      <<exeName<<endl
      <<"\t -f filelist.txt"<<endl
      <<"\t -s samplename"<<endl
      <<endl;
}
//----------------------------------------------------------
int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();
  string sampleName;
  string fileListName;
  bool verbose(false);

  int optind(1);
  while ((optind < argc)) {
    if(argv[optind][0]!='-'){optind++; continue;}
    string sw = argv[optind];
    if     (sw == "-h"){ printHelp(argv[0]); return 0; }
    else if(sw == "-f"){ optind++; fileListName = argv[optind]; }
    else if(sw == "-s"){ optind++; sampleName   = argv[optind]; }
    else if(sw == "-v"){ verbose = true; }
    else if(argv[optind][0]=='-') cout<<"Unknown switch "<<sw<<endl;
    optind++;
  } // end if(optind<argc)
  cout<<"Using the following options:"<<endl
      <<"fileListName   : "<<fileListName<<endl
      <<"sample         : "<<sampleName<<endl
      <<"verbose        : "<<verbose<<endl
      <<endl;

  TChain chain("physics");
  int fileErr = ChainHelper::addFileList(&chain, fileListName);
  if(!fileErr){
      Long64_t nEntries = chain.GetEntries();
      Long64_t nEntriesToProcess = 10;
      if(nEntries<nEntriesToProcess) nEntriesToProcess = nEntries;
      chain.ls();
      susy::D3PDAna analysis;
      if(verbose) analysis.setDebug(1);
      analysis.setSample(sampleName);
      cout<<"Total number of entries      : "<<nEntries<<endl;
      cout<<"Number of entries to process : "<<nEntriesToProcess<<endl;
      chain.Process(&analysis, sampleName.c_str(), nEntriesToProcess);
  } else {
      cout<<"error when adding files to the list"<<endl;
  }

  return 0;
}
