
#include <cstdlib>
#include <string>

#include "TChain.h"
#include "Cintex/Cintex.h"
#include "TSystem.h"

#include "SusyCommon/SusyNtMaker.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/ChainHelper.h"

using namespace std;

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
  bool writeNt    = true;

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
  cout << "  lumi    " << lumi     << endl;
  cout << "  xsec    " << xsec     << endl;
  cout << endl;

  // Build the input chain
  TChain* chain = new TChain("susy");
  int fileErr = ChainHelper::addFileList(chain, fileList);
  if(fileErr) return 1;
  Long64_t nEntries = chain->GetEntries();
  chain->ls();

  // Build the TSelector
  SusyNtMaker* susyAna = new SusyNtMaker();
  susyAna->setDebug(dbg);
  susyAna->setSample(sample);
  susyAna->setLumi(lumi);
  susyAna->setSumw(sumw);
  susyAna->setSys(sysOn);
  susyAna->setSavePhotons(savePh);
  susyAna->setXsec(xsec);
  susyAna->setFillNt(writeNt);


  // GRL
  TString grl = gSystem->ExpandPathName("$ROOTCOREDIR/data/MultiLep/data12_8TeV.periodAllYear_DetStatus-v46-pro13_CoolRunQuery-00-04-08_Susy.xml");
  //TString grl = gSystem->ExpandPathName("$ROOTCOREDIR/data/MultiLep/data11_7TeV.periodAllYear_DetStatus-v36-pro10_CoolRunQuery-00-04-08_Susy.xml");
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
