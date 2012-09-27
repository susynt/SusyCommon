
#include <cstdlib>
#include <string>

#include "TChain.h"
#include "Cintex/Cintex.h"
#include "TSystem.h"

#include "SusyCommon/SusyNtMaker.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/ChainHelper.h"

using namespace std;
using namespace Susy;

/*

    SusyNtMaker - an executable for running a basic looper/dumper of SUSY D3PDs
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

  cout << "  -s sample name, sets isMC flag"    << endl;
  cout << "     use e.g. 'ttbar', 'DataG', etc" << endl;

  cout << "  -w set sum of mc weights for norm" << endl;
  cout << "     default: 1"                     << endl;

  cout << "  -x set cross section"              << endl;
  cout << "     default: -1 (use susy db)"      << endl;

  cout << "  -l set lumi"                       << endl;
  cout << "     default: 5312/pb"               << endl;

  cout << "  -m write output ntuple"            << endl;
  cout << "     default: 1 (true)"              << endl;

  cout << "  --sys will turn on systematic run" << endl;
  cout << "     default: off"                   << endl;

  cout << "  --savePh will save photons"        << endl;
  cout << "     default: off"                   << endl;

  cout << "  --saveTau will save taus"          << endl;
  cout << "     default: off"                   << endl;

  cout << "  --af2 specifies AF2 samples"       << endl;
  cout << "     default: off"                   << endl;

  cout << "  --d3pd1032 sets d3pd tag to p1032" << endl;
  cout << "     default: p1181"                 << endl;

  cout << "  --metFlav set met flavor"          << endl;
  cout << "     default: Egamma10NoTau"         << endl;

  cout << "  -h print this help"                << endl;
}


int main(int argc, char** argv)
{
  ROOT::Cintex::Cintex::Enable();

  int nEvt        = -1;
  int nSkip       = 0;
  int dbg         = 0;
  float lumi      = 5831;
  float xsec      = -1;
  float sumw      = 1;
  string sample;
  string fileList = "fileList.txt";
  bool sysOn      = false;
  bool savePh     = false;
  bool saveTau    = false;
  bool isAF2      = false;
  bool writeNt    = true;
  D3PDTag tag     = D3PD_p1181;
  TString metFlav = "Egamma10NoTau";

  cout << "SusyNtMaker" << endl;
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
    else if (strcmp(argv[i], "-s") == 0)
      sample = argv[++i];
    else if (strcmp(argv[i], "-w") == 0)
      sumw = atof(argv[++i]);
    else if (strcmp(argv[i], "-x") == 0)
      xsec = atof(argv[++i]);
    else if (strcmp(argv[i], "-l") == 0)
      lumi = atof(argv[++i]);
    else if (strcmp(argv[i], "-m") == 0)
      writeNt = atoi(argv[++i]);
    else if (strcmp(argv[i], "--sys") == 0)
      sysOn = true;
    else if (strcmp(argv[i], "--savePh") == 0)
      savePh = true;
    else if (strcmp(argv[i], "--saveTau") == 0)
      saveTau = true;
    else if (strcmp(argv[i], "--af2") == 0)
      isAF2 = true;
    else if (strcmp(argv[i], "--d3pd1032") == 0)
      tag = D3PD_p1032;
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
  cout << "  sys     " << sysOn    << endl;
  cout << "  savePh  " << savePh   << endl;
  cout << "  saveTau " << saveTau  << endl;
  cout << "  isAF2   " << isAF2    << endl;
  cout << "  d3pdtag " << tag      << endl;
  cout << "  metFlav " << metFlav  << endl;
  cout << "  lumi    " << lumi     << endl;
  cout << "  xsec    " << xsec     << endl;
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
  SusyNtMaker* susyAna = new SusyNtMaker();
  susyAna->setDebug(dbg);
  susyAna->setSample(sample);
  susyAna->setLumi(lumi);
  susyAna->setSumw(sumw);
  susyAna->setSys(sysOn);
  //susyAna->setSavePhotons(savePh);
  susyAna->setSelectPhotons(savePh);
  susyAna->setSelectTaus(saveTau);
  susyAna->setAF2(isAF2);
  susyAna->setXsec(xsec);
  susyAna->setFillNt(writeNt);
  susyAna->setD3PDTag(tag);
  susyAna->setMetFlavor(metFlav);

  // GRL
  TString grl = gSystem->ExpandPathName("$ROOTCOREDIR/data/MultiLep/data12_8TeV.periodAllYear_DetStatus-v53-pro13-04_CoolRunQuery-00-04-08_All_Good_HCP.xml");
  susyAna->setGRLFile(grl);

  // Run the job
  if(nEvt<0) nEvt = nEntries;
  cout << endl;
  cout << "Total entries:   " << nEntries << endl;
  cout << "Process entries: " << nEvt << endl;
  chain->Process(susyAna, sample.c_str(), nEvt, nSkip);

  cout << endl;
  cout << "SusyNtMaker job done" << endl;

  delete chain;
  return 0;
}
