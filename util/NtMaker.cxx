
#include <cstdlib>
#include <string>

#include "TChain.h"
#include "TSystem.h"

#include "SusyCommon/SusyNtMaker.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"


using namespace std;
using namespace Susy;
using Susy::SusyNtMaker;

/*

    SusyNtMaker - an executable for running a basic looper/dumper of SUSY D3PDs
    It is not intended to be very advanced but more to look quickly at a d3pd and test the code

*/

void help()
{
  cout<<"Options:"                                                    <<endl
      <<"  -i|--input       name of the input container"              <<endl
      <<"  -o|--output      name of the output container"             <<endl
      <<"  -t|--tag         SustNtuple production tag"                <<endl
      <<"  -n|--num-events  default: -1 (all events)"                 <<endl
      <<"  -k|--num-skip    default: 0"                               <<endl
      <<"  -d|--dbg-level   default: 0 (quiet) "                      <<endl
      <<"  -f|--filelist    default: fileList.txt"                    <<endl
      <<"  -s|--sample      default: ''"                              <<endl
      <<"  -h|--help        print this help"                          <<endl
      <<"  -l|--lumi        default: 5312/pb"                         <<endl /// \todo obsolete option?
      <<"  -m|--write-nt    default: 1 (true, write tuple)"           <<endl
      <<"  --outfile-name   default: 'susyNt.root'"                   <<endl
      <<"  --grl            default: XaodAnalysis::defaultGrlFile()"  <<endl
      <<"  --sys            default: off"                             <<endl
      <<"  --savePh         save photons"                             <<endl
      <<"  --saveTau        save taus"                                <<endl
      <<"  --saveContTau    save container taus (instead of selected)"<<endl
      <<"  --saveTruth      save truth info"                          <<endl
      <<"  --af2            toggle atlas-fast-2 sample"               <<endl
      <<"  --filterOff      save all events"                          <<endl
      <<"  --nLepFilter     default 0 (require N  light leptons)"     <<endl
      <<"  --nLepTauFilter  default 2 (require N lepton+tau)"         <<endl
      <<"  --mc15b          toggle mc15b (default: false)"            <<endl
      <<endl;
}


int main(int argc, char** argv)
{
  int nEvt        = -1;
  int nSkip       = 0;
  int dbg         = 0;
  float lumi      = 5831;
  string sample   = "";
  string fileList = "fileList.txt";
  string grl      = "";
  bool sysOn      = false;
  bool savePh     = false;
  bool saveTau    = false;
  bool saveContTau= false;
  bool saveTruth  = false;
  bool isAF2      = false;
  bool writeNt    = true;
  bool filter     = true;
  uint nLepFilter = 0;
  uint nLepTauFilter = 2;
  bool mc15b      = false;
  string inputContainer, outputContainer, ntTag;
  string outputFileName = "susyNt.root";

  cout<<"SusyNtMaker"<<endl;
  cout<<endl;

  int optind(1);
  while(optind < argc) {
      std::string sw = argv[optind];
      if      (sw=="-i" || sw=="--input"     ) { inputContainer = argv[++optind]; }
      else if (sw=="-o" || sw=="--output"    ) { outputContainer = argv[++optind]; }
      else if (sw=="-t" || sw=="--tag"       ) { ntTag = argv[++optind]; }
      else if (sw=="-n" || sw=="--num-events") { nEvt = atoi(argv[++optind]); }
      else if (sw=="-k" || sw=="--num-skip"  ) { nSkip = atoi(argv[++optind]); }
      else if (sw=="-d" || sw=="--dbg-level" ) { dbg = atoi(argv[++optind]); }
      else if (sw=="-f" || sw=="--filelist"  ) { fileList = argv[++optind]; }
      else if (sw=="-s" || sw=="--sample"    ) { sample = argv[++optind]; }
      else if (sw=="-h" || sw=="--help"      ) { help(); return 0; }
      else if (sw=="-l" || sw=="--lumi"      ) { lumi = atof(argv[++optind]); }
      else if (sw=="-m" || sw=="--write-nt"  ) { writeNt = atoi(argv[++optind]); }
      else if (sw=="--outfile-name"          ) { outputFileName = argv[++optind]; } 
      else if (sw=="--grl"          ) { grl = argv[++optind]; }
      else if (sw=="--sys"          ) { sysOn = true; }
      else if (sw=="--savePh"       ) { savePh = true; }
      else if (sw=="--saveTau"      ) { saveTau = true; }
      else if (sw=="--saveContTau"  ) { saveTau = saveContTau = true; }
      else if (sw=="--saveTruth"    ) { saveTruth = true; }
      else if (sw=="--af2"          ) { isAF2 = true; }
      else if (sw=="--filterOff"    ) { filter = false; }
      else if (sw=="--nLepFilter"   ) { nLepFilter = atoi(argv[++optind]); }
      else if (sw=="--nLepTauFilter") { nLepTauFilter = atoi(argv[++optind]); }
      else if (sw=="--mc15b"        ) { mc15b = true; }
      else {
          cout<<"Unknown switch '"<<sw<<"'"<<endl;
          help();
          return 1;
      }
      optind++;
  } // while(optind)

  cout<<"flags:"<<endl;
  cout<<"  sample         "<<sample  <<endl;
  cout<<"  mc15b          "<<mc15b   <<endl;
  cout<<"  nEvt           "<<nEvt    <<endl;
  cout<<"  nSkip          "<<nSkip   <<endl;
  cout<<"  dbg            "<<dbg     <<endl;
  cout<<"  fileList       "<<fileList<<endl;
  cout<<"  grl            "<<grl     <<endl;
  cout<<"  sys            "<<sysOn   <<endl;
  cout<<"  savePh         "<<savePh  <<endl;
  cout<<"  saveTau        "<<saveTau <<endl;
  cout<<"  saveContTau    "<<saveContTau<<endl;
  cout<<"  saveTru        "<<saveTruth<< endl;
  cout<<"  isAF2          "<<isAF2   <<endl;
  cout<<"  lumi           "<<lumi    <<endl;
  cout<<"  filter         "<<filter  <<endl;
  cout<<"  nLepFilter     "<<nLepFilter   <<endl;
  cout<<"  nLepTauFilter  "<<nLepTauFilter<<endl;
  cout<<"  input          "<<inputContainer <<endl;
  cout<<"  output         "<<outputContainer<<endl;
  cout<<"  ntTag          "<<ntTag          <<endl;
  cout<<"  outputFileName "<<outputFileName <<endl;
  cout<<endl;


  // Build the input chain
  TChain* chain = new TChain("CollectionTree");
  int fileErr = ChainHelper::addInput(chain, fileList, dbg);
  if(fileErr) return 1;
  Long64_t nEntries = chain->GetEntries();
  chain->ls();

  // Build the TSelector
  SusyNtMaker* susyAna = new SusyNtMaker();
  susyAna->setChain(chain);
  susyAna->setDebug(dbg);
  //susyAna->setSample(sample);
  susyAna->setMC15b(mc15b);
  susyAna->setLumi(lumi);
  susyAna->setSys(sysOn);
  susyAna->setSelectPhotons(savePh);
  susyAna->setSelectTaus(saveTau);
  susyAna->setSaveContTaus(saveContTau);
  susyAna->setAF2(isAF2);
  susyAna->setFillNt(writeNt);
  susyAna->setSelectTruthObjects(saveTruth);
  susyAna->setFilter(filter);
  susyAna->setNLepFilter(nLepFilter);
  susyAna->setNLepTauFilter(nLepTauFilter);
  susyAna->m_inputContainerName = inputContainer;
  susyAna->m_outputContainerName = outputContainer;
  susyAna->m_productionTag = ntTag;
  susyAna->m_outputFileName = outputFileName;
  susyAna->m_productionCommand = Susy::utils::commandLineArguments(argc, argv);
  // GRL - default is set in SusyD3PDAna::Begin, but now we can override it here
  susyAna->setGRLFile(grl);

  // Run the job
  if(nEvt<0) nEvt = nEntries;
  cout<<endl;
  cout<<"Total entries:   "<<nEntries<<endl;
  cout<<"Process entries: "<<nEvt<<endl;
  chain->Process(susyAna, sample.c_str(), nEvt, nSkip);

  cout<<endl;
  cout<<"SusyNtMaker job done"<<endl;

  delete chain;
  return 0;
}
