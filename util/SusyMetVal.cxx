
#include <cstdlib>
#include <string>

#include "TChain.h"
#include "Cintex/Cintex.h"
#include "TSystem.h"

#include "SusyCommon/SusyMetValidation.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/ChainHelper.h"

using namespace std;

/*

    SusyMetValidation - an executable for running a basic looper/dumper of SUSY D3PDs
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

  cout << "  -o name of output hist file"       << endl;
  cout << "     defaults: ''"                   << endl;

  cout << "  -s sample name, sets isMC flag"    << endl;
  cout << "     use e.g. 'ttbar', 'DataG', etc" << endl;

  cout << "  -w set sum of mc weights for norm" << endl;
  cout << "     default: 1"                     << endl;

  cout << "  -x set cross section"              << endl;
  cout << "     default: -1 (use susy db)"      << endl;

  //cout << "  -l set lumi"                       << endl;
  //out << "     default: 5312/pb"               << endl;

  cout << "  --af2 specifies AF2 samples"       << endl;
  cout << "     default: off"                   << endl;

  cout << "  --metFlav set met flavor"          << endl;
  cout << "     default: Egamma10NoTau_STVF"    << endl;

  cout << "  -h print this help"                << endl;
}


int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();

  int nEvt        = -1;
  int nSkip       = 0;
  int dbg         = 0;
  string sample;
  string fileList = "fileList.txt";
  float xsec      = -1;
  float sumw      = 1;
  bool isAF2      = false;
  TString metFlav = "Egamma10NoTau_STVF";
  string histFileName = "";

  cout << "SusyMetValidation" << endl;
  cout << endl;

  // Read inputs to program
  for(int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-n") == 0)
      nEvt = atoi(argv[++i]);
    else if (strcmp(argv[i], "-k") == 0)
      nSkip = atoi(argv[++i]);
    else if (strcmp(argv[i], "-d") == 0)
      dbg = atoi(argv[++i]);
    else if (strcmp(argv[i], "-f") == 0)
      fileList = argv[++i];
    else if (strcmp(argv[i], "-o") == 0)
      histFileName = argv[++i];
    else if (strcmp(argv[i], "-s") == 0)
      sample = argv[++i];
    else if (strcmp(argv[i], "-w") == 0)
      sumw = atof(argv[++i]);
    else if (strcmp(argv[i], "-x") == 0)
      xsec = atof(argv[++i]);
    else if (strcmp(argv[i], "--af2") == 0)
      isAF2 = true;
    else if (strcmp(argv[i], "--metFlav") == 0)
      metFlav = argv[++i];
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
  cout << "  sumw    " << sumw     << endl;
  cout << "  xsec    " << xsec     << endl;
  cout << "  isAF2   " << isAF2    << endl;
  cout << "  metFlav " << metFlav  << endl;
  cout << endl;

  // Build the input chain
  TChain* chain = new TChain("susy");
  int fileErr = ChainHelper::addFileList(chain, fileList);
  if(fileErr) return 1;
  Long64_t nEntries = chain->GetEntries();
  chain->ls();

  // Try to catch cache errors that give empty files
  if(nEntries == 0 && nEvt != 0){
    cerr << "ERROR - zero input entries! Probably a problem. Aborting!" << endl;
    abort();
  }

  // Build the TSelector
  SusyMetValidation* susyAna = new SusyMetValidation();
  susyAna->setHistFileName(histFileName);
  susyAna->setDebug(dbg);
  susyAna->setSample(sample);
  susyAna->setSumw(sumw);
  susyAna->setAF2(isAF2);
  susyAna->setXsec(xsec);
  susyAna->setMetFlavor(metFlav);

  // GRL - now set in SusyD3PDAna
  //TString grl = gSystem->ExpandPathName("$ROOTCOREDIR/data/MultiLep/data12_8TeV.periodAllYear_DetStatus-v53-pro13-04_CoolRunQuery-00-04-08_All_Good_HCP.xml");
  //susyAna->setGRLFile(grl);

  // Run the job
  if(nEvt<0) nEvt = nEntries;
  cout << endl;
  cout << "Total entries:   " << nEntries << endl;
  cout << "Process entries: " << nEvt << endl;
  chain->Process(susyAna, sample.c_str(), nEvt, nSkip);

  cout << endl;
  cout << "SusyMetValidation job done" << endl;

  delete chain;
  return 0;
}
