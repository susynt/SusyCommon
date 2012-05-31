
#include <cstdlib>
#include <string>

#include "TChain.h"
#include "Cintex/Cintex.h"

#include "SusyCommon/SusyD3PDInterface.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/ChainHelper.h"

using namespace std;

/*

    SusyD3PDTest - an executable for running a basic looper/dumper of SUSY D3PDs
    It is not intended to be very advanced but more to look quickly at a d3pd and test the code

*/

void help()
{
  cout << "  Options:"                          << endl;
  cout << "  -n number of events to process"    << endl;
  cout << "     defaults: -1 (all events)"      << endl;

  cout << "  -k number of events to skip"       << endl;
  cout << "     defaults: 0"                    << endl;

  cout << "  -d debug printout level"           << endl;
  cout << "     defaults: 0 (quiet) "           << endl;

  cout << "  -f name of input filelist"         << endl;
  cout << "     defaults: fileList.txt"         << endl;

  cout << "  -s sample name, for naming files"  << endl;
  cout << "     defaults: ntuple sample name"   << endl;

  cout << "  -h print this help"                << endl;
}


int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();

  int nEvt = -1;
  int nSkip = 0;
  int dbg = 0;
  string sample;
  string fileList = "fileList.txt";
  
  cout << "SusyD3PDTest" << endl;
  cout << endl;

  /** Read inputs to program */
  for(int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0)
      nEvt = atoi(argv[++i]);
    else if (strcmp(argv[i], "-k") == 0)
      nSkip = atoi(argv[++i]);
    else if (strcmp(argv[i], "-d") == 0)
      dbg = atoi(argv[++i]);
    else if (strcmp(argv[i], "-f") == 0)
      fileList = argv[++i];
    else if (strcmp(argv[i], "-s") == 0)
      sample = argv[++i];
    //if (strcmp(argv[i], "-h") == 0)
    else
    {
      help();
      return 0;
    }
  }

  cout << "flags:" << endl;
  cout << "  sample  " << sample   << endl;
  cout << "  nEvt    " << nEvt     << endl;
  cout << "  nSkip   " << nSkip    << endl;
  cout << "  dbg     " << dbg      << endl;
  cout << "  input   " << fileList << endl;
  cout << endl;

  // Build the input chain
  TChain* chain = new TChain("susy");
  int fileErr = ChainHelper::addFileList(chain, fileList);
  if(fileErr) return 1;
  Long64_t nEntries = chain->GetEntries();
  chain->ls();

  // Build the TSelector
  SusyD3PDInterface* susyAna = new SusyD3PDInterface();
  susyAna->setDebug(dbg);

  // Run the job
  if(nEvt<0) nEvt = nEntries;
  cout << endl;
  cout << "Total entries:   " << nEntries << endl;
  cout << "Process entries: " << nEvt << endl;
  chain->Process(susyAna, sample.c_str(), nEvt, nSkip);

  cout << endl;
  cout << "SusyD3PDTest job done" << endl;

  delete chain;
  return 0;
}
