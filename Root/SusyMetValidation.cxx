#include "egammaAnalysisUtils/CaloIsoCorrection.h"
#include "MultiLep/MuonTools.h"
#include "MultiLep/ElectronTools.h"
#include "MultiLep/CutflowTools.h"
#include "SusyCommon/SusyMetValidation.h"

#define GeV 1000.

using namespace std;

const bool dumpMet = true;

/*--------------------------------------------------------------------------------*/
// SusyMetValidation Constructor
/*--------------------------------------------------------------------------------*/
SusyMetValidation::SusyMetValidation()
{
  n_base_ele=0;
  n_base_muo=0;
  n_base_tau=0;
  n_base_jet=0;
  n_sig_ele=0;
  n_sig_muo=0;
  n_sig_tau=0;
  n_sig_jet=0;
  n_evt_initial=0;
  n_evt_grl=0;
  n_evt_larErr=0;
  n_evt_larHole=0;
  n_evt_hotSpot=0;
  n_evt_badJet=0;
  n_evt_goodVtx=0;
  n_evt_badMu=0;
  n_evt_cosmic=0;
  n_evt_1Lep=0;
  n_evt_2Lep=0;
  n_evt_ee=0;
  n_evt_sfos=0;
  n_evt_zMass=0;
}
/*--------------------------------------------------------------------------------*/
// Destructor
/*--------------------------------------------------------------------------------*/
SusyMetValidation::~SusyMetValidation()
{
}
/*--------------------------------------------------------------------------------*/
// The Begin() function is called at the start of the query.
// When running with PROOF Begin() is only called on the client.
// The tree argument is deprecated (on PROOF 0 is passed).
/*--------------------------------------------------------------------------------*/
void SusyMetValidation::Begin(TTree* /*tree*/)
{
  SusyD3PDAna::Begin(0);
  if(m_dbg) cout << "SusyMetValidation::Begin" << endl;

  // Start the timer
  m_timer.Start();
}

/*--------------------------------------------------------------------------------*/
// Main process loop function - This is just an example for testing
/*--------------------------------------------------------------------------------*/
Bool_t SusyMetValidation::Process(Long64_t entry)
{
  // Communicate the entry number to the interface objects
  GetEntry(entry);

  static Long64_t chainEntry = -1;
  chainEntry++;
  if(m_dbg || chainEntry%5000==0)
  {
    cout << "**** Processing entry " << setw(6) << chainEntry
         << " run " << setw(6) << d3pd.evt.RunNumber()
         << " event " << setw(7) << d3pd.evt.EventNumber() << " ****" << endl;
  }

  if(selectEvent()){
    checkMet();
  }

  return kTRUE;
}

/*--------------------------------------------------------------------------------*/
// The Terminate() function is the last function to be called during
// a query. It always runs on the client, it can be used to present
// the results graphically or save the results to file.
/*--------------------------------------------------------------------------------*/
void SusyMetValidation::Terminate()
{
  // Stop the timer
  m_timer.Stop();

  SusyD3PDAna::Terminate();
  if(m_dbg) cout << "SusyMetValidation::Terminate" << endl;

  cout << endl;
  cout << "Object counter" << endl;
  cout << "  BaseEle  " << n_base_ele    << endl;
  cout << "  BaseMuo  " << n_base_muo    << endl;
  cout << "  BaseTau  " << n_base_tau    << endl;
  cout << "  BaseJet  " << n_base_jet    << endl;
  cout << "  SigEle   " << n_sig_ele     << endl;
  cout << "  SigMuo   " << n_sig_muo     << endl;
  cout << "  SigTau   " << n_sig_tau     << endl;
  cout << "  SigJet   " << n_sig_jet     << endl;
  cout << endl;
  cout << "Event counter" << endl;
  cout << "  Initial  " << n_evt_initial << endl;
  cout << "  GRL      " << n_evt_grl     << endl;
  cout << "  LarErr   " << n_evt_larErr  << endl;
  cout << "  HotSpot  " << n_evt_hotSpot << endl;
  cout << "  BadJet   " << n_evt_badJet  << endl;
  cout << "  GoodVtx  " << n_evt_goodVtx << endl;
  cout << "  BadMuon  " << n_evt_badMu   << endl;
  cout << "  Cosmic   " << n_evt_cosmic  << endl;
  cout << "  >=1 Lep  " << n_evt_1Lep    << endl;
  cout << "  >=2 Lep  " << n_evt_2Lep    << endl;
  cout << "  Is EE    " << n_evt_ee      << endl;
  cout << "  SFOS     " << n_evt_sfos    << endl;
  cout << "  Z Mass   " << n_evt_zMass   << endl;
  cout << endl;

  // Report timer
  double realTime = m_timer.RealTime();
  double cpuTime  = m_timer.CpuTime();
  int hours = int(realTime / 3600);
  realTime -= hours * 3600;
  int min   = int(realTime / 60);
  realTime -= min * 60;
  int sec   = int(realTime);

  float speed = n_evt_initial/m_timer.RealTime()/1000;

  printf("---------------------------------------------------\n");
  printf(" Number of events processed: %d \n",n_evt_initial);
  printf("\t Analysis time: Real %d:%02d:%02d, CPU %.3f      \n", hours, min, sec, cpuTime);
  printf("\t Analysis speed [kHz]: %2.3f                     \n",speed);
  printf("---------------------------------------------------\n\n");
}

/*--------------------------------------------------------------------------------*/
// Select event
/*--------------------------------------------------------------------------------*/
bool SusyMetValidation::selectEvent()
{
  if(m_dbg) cout << "selectEvent" << endl;

  n_evt_initial++;

  m_susyObj.Reset();
  clearObjects();
 
  // grl
  if(!passGRL()) return false;
  n_evt_grl++;

  // larErr
  if(!passLarErr()) return false;
  n_evt_larErr++;

  // Get Nominal Objects
  selectObjects();
  buildMet();
  //evtCheck();
  checkEventCleaning();
  checkObjectCleaning();

  // Tile hot spot
  //if((m_evtFlag & PASS_HotSpot) == 0) return false;
  if((m_cutFlags & ECut_HotSpot) == 0) return false;
  n_evt_hotSpot++;
  // Bad jet cut
  if((m_cutFlags & ECut_BadJet) == 0) return false;
  n_evt_badJet++;

  // primary vertex cut 
  if(!passGoodVtx()) return false; 
  n_evt_goodVtx++;

  // Bad muon veto
  if((m_cutFlags & ECut_BadMuon) == 0) return false;
  n_evt_badMu++;
  // Cosmic muon veto
  if((m_cutFlags & ECut_Cosmic) == 0) return false;
  n_evt_cosmic++;
  
  n_base_ele += m_baseElectrons.size();
  n_base_muo += m_baseMuons.size();
  n_base_tau += m_baseTaus.size();
  n_base_jet += m_baseJets.size();
  n_sig_ele += m_sigElectrons.size();
  n_sig_muo += m_sigMuons.size();
  n_sig_tau += m_sigTaus.size();
  n_sig_jet += m_sigJets.size();
  
  // Lepton multiplicity
  uint nBaseLep = m_baseElectrons.size() + m_baseMuons.size();
  uint nSigLep = m_sigElectrons.size() + m_sigMuons.size();
  //cout << "nSigLep " << nSigLep << endl;
  if(nSigLep < 1) return false;
  n_evt_1Lep++;
  if(nBaseLep != 2) return false;
  if(nSigLep != 2) return false;
  n_evt_2Lep++;

  // Get the signal leptons
  const LeptonInfo* lep1 = & m_sigLeptons[0];
  const LeptonInfo* lep2 = & m_sigLeptons[1];

  // For now, only looking at the ee channel
  if(!lep1->isElectron() || !lep2->isElectron()) return false;
  n_evt_ee++;

  if(!IsSFOS(lep1, lep2)) return false;
  n_evt_sfos++;

  float mll = GetMll(lep1, lep2)/GeV;
  if(fabs(mll-MZ) > 10) return false;
  n_evt_zMass++;

  return true;
}

/*--------------------------------------------------------------------------------*/
// Check the MET
/*--------------------------------------------------------------------------------*/
void SusyMetValidation::checkMet()
{
  if(dumpMet){
  
    cout << endl;
    dumpEvent();
    dumpBaselineObjects();
    //dumpSignalObjects();

    // Compare the total MET before and after SUSYTools
    cout << "Total MET" << endl;

    // The D3PD MET
    cout << "  D3PD: " << d3pd.met.RefFinal_et()/GeV << endl;
    cout << "  SUSY: " << m_met.Et()/GeV << endl;

  }
}


#undef GeV
