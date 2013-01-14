#include "egammaAnalysisUtils/CaloIsoCorrection.h"
#include "MultiLep/MuonTools.h"
#include "MultiLep/ElectronTools.h"
#include "MultiLep/CutflowTools.h"
#include "SusyCommon/SusyMetValidation.h"

#define GeV 1000.

using namespace std;

const bool dumpMet = false;

/*--------------------------------------------------------------------------------*/
// SusyMetValidation Constructor
/*--------------------------------------------------------------------------------*/
SusyMetValidation::SusyMetValidation()
{
  // Met names
  metNames[Met_susy_stvf] = "susy_stvf";
  metNames[Met_d3pd_stvf] = "d3pd_stvf";
  metNames[Met_d3pd_reff] = "d3pd_reff";

  // Initialize counters
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
  n_evt_hfor=0;
  n_evt_larErr=0;
  n_evt_larHole=0;
  n_evt_hotSpot=0;
  n_evt_badJet=0;
  n_evt_goodVtx=0;
  n_evt_badMu=0;
  n_evt_cosmic=0;
  n_evt_3Lep=0;
  n_evt_trig=0;
  n_evt_sfos=0;
  n_evt_zMass=0;
  n_evt_met=0;
  n_evt_bJet=0;
  n_evt_mt=0;
  n_evt_lepPt=0;
  n_evt_2mu=0;
  n_evt_pileup=0;
  n_evt_lepSF=0;
  n_evt_bTagSF=0;
  n_evt_lumi=0;
}
/*--------------------------------------------------------------------------------*/
// Destructor
/*--------------------------------------------------------------------------------*/
SusyMetValidation::~SusyMetValidation()
{
  delete m_triggerMatch;
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

  // Initialize trigger match tool
  m_triggerMatch = new TriggerMatchMultiLep();

  // Histograms
  bookHistos();

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

  // Event selection
  if(selectEvent()){
    d3pd.ReadAllActive();
    m_outputTree->Fill();
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
  cout << "  HFOR     " << n_evt_hfor    << endl;
  cout << "  LarErr   " << n_evt_larErr  << endl;
  cout << "  HotSpot  " << n_evt_hotSpot << endl;
  cout << "  BadJet   " << n_evt_badJet  << endl;
  cout << "  GoodVtx  " << n_evt_goodVtx << endl;
  cout << "  BadMuon  " << n_evt_badMu   << endl;
  cout << "  Cosmic   " << n_evt_cosmic  << endl;
  cout << "  ==3 Lep  " << n_evt_3Lep    << endl;
  cout << "  Trig     " << n_evt_trig    << endl;
  
  // Add here the rest of the event selection
  cout << "  SFOS     " << n_evt_sfos    << endl;
  cout << "  Z Mass   " << n_evt_zMass   << endl;
  cout << "  MET      " << n_evt_met     << endl;
  cout << "  BJet     " << n_evt_bJet    << endl;
  cout << "  MT       " << n_evt_mt      << endl;
  cout << "  LepPt    " << n_evt_lepPt   << endl;
  cout << "  2 Muons  " << n_evt_2mu     << endl;

  // Weighted events
  //cout << "  Pileup   " << n_evt_pileup  << endl;
  cout << "  LepSF    " << n_evt_lepSF   << endl;
  cout << "  BTagSF   " << n_evt_bTagSF  << endl;
  cout << "  A-E Lumi " << n_evt_lumi    << endl;

  cout << endl;

  // Report timer
  double realTime = m_timer.RealTime();
  double cpuTime  = m_timer.CpuTime();
  int hours = int(realTime / 3600);
  realTime -= hours * 3600;
  int min   = int(realTime / 60);
  realTime -= min * 60;
  int sec   = int(realTime);

  float speed = n_evt_initial/m_timer.RealTime();

  printf("---------------------------------------------------\n");
  printf(" Number of events processed: %d \n",n_evt_initial);
  printf("\t Analysis time: Real %d:%02d:%02d, CPU %.3f      \n", hours, min, sec, cpuTime);
  printf("\t Analysis speed [Hz]: %2.1f                      \n",speed);
  printf("---------------------------------------------------\n\n");

  // Save histograms
  saveHistos();
}

/*--------------------------------------------------------------------------------*/
// Book histograms
/*--------------------------------------------------------------------------------*/
void SusyMetValidation::bookHistos()
{
  system("mkdir -p rootFiles/MetVal");
  if(m_histFileName.empty()) m_histFileName = "rootFiles/MetVal/"+m_sample+".metVal.root";
  m_histFile = new TFile(m_histFileName.c_str(), "recreate");
  TH1::SetDefaultSumw2(true);

  // Preprocessor convenience
  #define NEWHIST(name, xLbl, nBin, min, max) \
    new TH1F(name, name ";" xLbl ";Events", nBin, min, max)

  // Met flavor histos
  for(uint i=0; i<Met_N; i++){
    string metName = metNames[i];

    #define NEWMETHIST(name, xLbl, nBin, min, max) \
      new TH1F((metName+"_"+name).c_str(), (metName+"_"+name+";"+xLbl+";Events").c_str(), nBin, min, max)
    #define NEWVARMETHIST(name, xLbl, nBin, bins) \
      new TH1F((metName+"_"+name).c_str(), (metName+"_"+name+";"+xLbl+";Events").c_str(), nBin, bins)

    // MET terms
    h_met[i] = NEWMETHIST("met", "MET [GeV]", 20, 0, 100);
    h_metEle[i] = NEWMETHIST("metEle", "MET Electron [GeV]", 20, 0, 100);
    h_metMuo[i] = NEWMETHIST("metMuo", "MET Muon [GeV]", 20, 0, 100);
    h_metJet[i] = NEWMETHIST("metJet", "MET Jet [GeV]", 20, 0, 100);
    h_metCell[i] = NEWMETHIST("metCell", "MET Cell [GeV]", 20, 0, 100);

    // MET leptons
    h_nMetLep[i] = NEWMETHIST("nMetLep", "Number of MET leptons", 10, -0.5, 9.5);
    h_nMetMu[i] = NEWMETHIST("nMetMu", "Number of MET muons", 10, -0.5, 9.5);
    h_nMetEl[i] = NEWMETHIST("nMetEl", "Number of MET electrons", 10, -0.5, 9.5);

    h_metLepPt[i] = NEWMETHIST("metLepPt", "MET lepton_{} P_{T} [GeV]", 25, 0, 500);
    h_metMuPt[i] = NEWMETHIST("metMuPt", "MET muon_{} P_{T} [GeV]", 25, 0, 500);
    h_metElPt[i] = NEWMETHIST("metElPt", "MET electron_{} P_{T} [GeV]", 25, 0, 500);

    h_metLepEta[i] = NEWMETHIST("metLepEta", "MET lepton_{} #eta", 20, -5.0, 5.0);
    //h_metMuEta[i] = NEWMETHIST("metMuEta", "MET muon_{} #eta", 10, -2.5, 2.5);
    h_metMuEta[i] = NEWMETHIST("metMuEta", "MET muon_{} #eta", 20, -4.8, 4.8);
    h_metElEta[i] = NEWMETHIST("metElEta", "MET electron_{} #eta", 20, -5.0, 5.0);

    // MET muon weights
    h_muWet[i] = NEWMETHIST("muWet", "MET muon E_{T} weight", 55, 0, 1.1);
    h_muWpx[i] = NEWMETHIST("muWpx", "MET muon P_{X} weight", 55, 0, 1.1);
    h_muWpy[i] = NEWMETHIST("muWpy", "MET muon P_{Y} weight", 55, 0, 1.1);
  }

  // Other histos
  //h_nMuRaw = NEWHIST("nMuRaw", "Number of D3PD muons", 10, -0.5, 9.5);
  h_nMu = NEWHIST("nMu", "Signal muon multiplicity", 6, -0.5, 5.5);

  h_rawMuPt = NEWHIST("rawMuPt", "D3PD raw muon p_{T} [GeV]", 25, 0, 500);
  h_rawMuEta = NEWHIST("rawMuEta", "D3PD raw muon #eta", 20, -4.8, 4.8);

  // Output d3pd, just to save the run and event numbers
  // Testing out the new skimmer!
  m_outputTree = new TTree("susy", "susy");
  d3pd.evt.SetActive(true, "RunNumber");
  d3pd.evt.SetActive(true, "EventNumber");
  if(m_isMC) d3pd.truth.SetActive(true, "mc_channel_number");
  d3pd.WriteTo(m_outputTree);
}

/*--------------------------------------------------------------------------------*/
// Save histograms
/*--------------------------------------------------------------------------------*/
void SusyMetValidation::saveHistos()
{
  m_histFile->Write();
  cout << endl << "Histograms written to " << m_histFile->GetName() << endl;
  m_histFile->Close();
}

/*--------------------------------------------------------------------------------*/
// Select event
// The selection is currently setup to evaluate the 3L met bump
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

  // HF overlap removal
  if(m_isMC && getHFORDecision()==4) return false;
  n_evt_hfor++;

  // larErr
  if(!passLarErr()) return false;
  n_evt_larErr++;

  // Get Nominal Objects
  selectObjects();
  buildMet();
  checkEventCleaning();
  checkObjectCleaning();

  // Tile hot spot
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
  //uint nBaseLep = m_baseElectrons.size() + m_baseMuons.size();
  //uint nSigEle = m_sigElectrons.size();
  uint nSigMuo = m_sigMuons.size();
  //uint nSigLep = nSigEle + nSigMuo;
  //if(nBaseLep != 3) return false;
  //if(nSigLep != 3) return false;
  n_evt_3Lep++;

  // Trigger cut
  //if(!passTrigger()) return false;
  n_evt_trig++;

  // SFOS
  //if(!HasSFOS(m_sigLeptons)) return false;
  n_evt_sfos++;

  // Z mass
  static vector<float> msfos;
  msfos.clear();
  msfos = MassesOfSFOSPairs(m_susyObj, &d3pd.muo, m_sigMuons, &d3pd.ele, m_sigElectrons);
  bool hasZ = false;
  for(uint i=0; i<msfos.size(); i++){
    if(fabs(msfos[i]-91.2*GeV) < 10*GeV){
      hasZ = true;
      break;
    }
  }
  //if(!hasZ) return false;
  n_evt_zMass++;

  // MET
  //if(m_met.Et() < 30*GeV || m_met.Et() > 75*GeV) return false;
  n_evt_met++;

  // Bjet veto
  //if(IsBJetEvent(m_susyObj, &d3pd.jet, m_sigJets, SUSYBTagger::MV1, make_pair("0_122", 0.122))) return false;
  n_evt_bJet++;

  // Mt
  //float mt = getMt(m_susyObj, &d3pd.muo, m_sigMuons, &d3pd.ele, m_sigElectrons, 
                   //m_met.Vect().XYvector(), m_met.Et());
  //if(mt < 50*GeV || mt > 110*GeV) return false;
  n_evt_mt++;

  // Lepton pt
  //if(m_sigLeptons[2].lv()->Pt() < 20*GeV) return false;
  n_evt_lepPt++;

  // Select emm and mmm channels by requiring at least 2 muons
  //if(m_sigMuons.size() < 2) return false;
  n_evt_2mu++;

  // Get the event weight
  float w = getEventWeight();
  // Lepton eff SF
  float lepSF = getLepSF(m_sigLeptons);
  n_evt_lepSF += lepSF;
  // B-tag eff SF
  float bTagSF = getBTagSF(m_sigJets);
  n_evt_bTagSF += bTagSF;
  // Full weight
  float wTotal = w * lepSF * bTagSF;
  n_evt_lumi += wTotal;

  //
  // Calculate all of the met variables we will want to look at here
  //

  // The MET terms
  TVector2 met[Met_N];
  TVector2 metEle[Met_N];
  TVector2 metMuo[Met_N];
  TVector2 metJet[Met_N];
  TVector2 metSoftJet[Met_N];
  TVector2 metCell[Met_N];

  // SUSYTools Egamma10NoTau_STVF_RefFinal
  met[Met_susy_stvf]     = m_met.Vect().XYvector();
  metEle[Met_susy_stvf]  = m_susyObj.computeMETComponent(METUtil::RefEle);
  metMuo[Met_susy_stvf]  = m_susyObj.computeMETComponent(METUtil::MuonTotal);
  metJet[Met_susy_stvf]  = m_susyObj.computeMETComponent(METUtil::RefJet);
  metCell[Met_susy_stvf] = m_susyObj.computeMETComponent(METUtil::CellOutEflow);
  metSoftJet[Met_susy_stvf] = m_susyObj.computeMETComponent(METUtil::SoftJets);

  // D3PD Egamma10NoTau_STVF_RefFinal
  met[Met_d3pd_stvf]     = TVector2(d3pd.met.Egamma10NoTau_RefFinal_STVF_etx(), 
                                    d3pd.met.Egamma10NoTau_RefFinal_STVF_ety());
  metEle[Met_d3pd_stvf]  = TVector2(d3pd.met.Egamma10NoTau_RefEle_etx(), 
                                    d3pd.met.Egamma10NoTau_RefEle_ety());
  metMuo[Met_d3pd_stvf]  = TVector2(d3pd.met.Egamma10NoTau_Muon_Total_Staco_etx(), 
                                    d3pd.met.Egamma10NoTau_Muon_Total_Staco_ety());
  metJet[Met_d3pd_stvf]  = TVector2(d3pd.met.Egamma10NoTau_RefJet_etx(), 
                                    d3pd.met.Egamma10NoTau_RefJet_ety());
  metCell[Met_d3pd_stvf] = TVector2(d3pd.met.Egamma10NoTau_CellOut_Eflow_STVF_etx(), 
                                    d3pd.met.Egamma10NoTau_CellOut_Eflow_STVF_ety());
  metSoftJet[Met_d3pd_stvf] = TVector2(d3pd.met.Egamma10NoTau_SoftJets_etx(), 
                                       d3pd.met.Egamma10NoTau_SoftJets_ety());

  // D3PD RefFinal
  met[Met_d3pd_reff]     = TVector2(d3pd.met.Egamma10NoTau_RefFinal_etx(), 
                                    d3pd.met.Egamma10NoTau_RefFinal_ety());
  metEle[Met_d3pd_reff]  = TVector2(d3pd.met.Egamma10NoTau_RefEle_etx(), 
                                    d3pd.met.Egamma10NoTau_RefEle_ety());
  metMuo[Met_d3pd_reff]  = TVector2(d3pd.met.Egamma10NoTau_Muon_Total_Staco_etx(), 
                                    d3pd.met.Egamma10NoTau_Muon_Total_Staco_ety());
  metJet[Met_d3pd_reff]  = TVector2(d3pd.met.Egamma10NoTau_RefJet_etx(), 
                                    d3pd.met.Egamma10NoTau_RefJet_ety());
  // Not sure if I should use the CellOut or the CellOut_Eflow.  One might be zero..
  metCell[Met_d3pd_reff] = TVector2(d3pd.met.Egamma10NoTau_CellOut_Eflow_etx(), 
                                    d3pd.met.Egamma10NoTau_CellOut_Eflow_ety());
  metSoftJet[Met_d3pd_reff] = TVector2(d3pd.met.Egamma10NoTau_SoftJets_etx(), 
                                       d3pd.met.Egamma10NoTau_SoftJets_ety());

  // Dump variables
  if(dumpMet || m_dbg>=2){
    cout << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    dumpEvent();
    dumpBaselineObjects();
    dumpSignalObjects();

    // Compare the total MET
    cout << "Total MET" << endl;
    cout << "  SUSY STVF:     " << met[Met_susy_stvf].Mod()/GeV << endl;
    cout << "  D3PD STVF:     " << met[Met_d3pd_stvf].Mod()/GeV << endl;
    cout << "  D3PD RefFinal: " << met[Met_d3pd_reff].Mod()/GeV  << endl;

    cout << "Ref Electron" << endl;
    cout << "  SUSY STVF:     " << metEle[Met_susy_stvf].Mod()/GeV << endl;
    cout << "  D3PD STVF:     " << metEle[Met_d3pd_stvf].Mod()/GeV << endl;
    cout << "  D3PD RefFinal: " << metEle[Met_d3pd_reff].Mod()/GeV  << endl;

    cout << "Total Muon" << endl;
    cout << "  SUSY STVF:     " << metMuo[Met_susy_stvf].Mod()/GeV << endl;
    cout << "  D3PD STVF:     " << metMuo[Met_d3pd_stvf].Mod()/GeV << endl;
    cout << "  D3PD RefFinal: " << metMuo[Met_d3pd_reff].Mod()/GeV  << endl;

    cout << "Ref Jet" << endl;
    cout << "  SUSY STVF:     " << metJet[Met_susy_stvf].Mod()/GeV << endl;
    cout << "  D3PD STVF:     " << metJet[Met_d3pd_stvf].Mod()/GeV << endl;
    cout << "  D3PD RefFinal: " << metJet[Met_d3pd_reff].Mod()/GeV  << endl;

    cout << "Cell Out" << endl;
    cout << "  SUSY STVF:     " << metCell[Met_susy_stvf].Mod()/GeV << endl;
    cout << "  D3PD STVF:     " << metCell[Met_d3pd_stvf].Mod()/GeV << endl;
    cout << "  D3PD RefFinal: " << metCell[Met_d3pd_reff].Mod()/GeV  << endl;

    cout << "Soft Jet" << endl;
    cout << "  SUSY STVF:     " << metSoftJet[Met_susy_stvf].Mod()/GeV << endl;
    cout << "  D3PD STVF:     " << metSoftJet[Met_d3pd_stvf].Mod()/GeV << endl;
    cout << "  D3PD RefFinal: " << metSoftJet[Met_d3pd_reff].Mod()/GeV  << endl;

    // TESTING, more explicit dump of D3PD vars
    cout << "Egamma10NoTau_RefFinal_STVF_etx        " << d3pd.met.Egamma10NoTau_RefFinal_STVF_etx()/GeV << endl;
    cout << "Egamma10NoTau_RefFinal_STVF_ety        " << d3pd.met.Egamma10NoTau_RefFinal_STVF_ety()/GeV << endl;
    cout << "Egamma10NoTau_RefFinal_STVF_sumet      " << d3pd.met.Egamma10NoTau_RefFinal_STVF_sumet()/GeV << endl;
    cout << "Egamma10NoTau_RefEle_etx               " << d3pd.met.Egamma10NoTau_RefEle_etx()/GeV << endl;
    cout << "Egamma10NoTau_RefEle_ety               " << d3pd.met.Egamma10NoTau_RefEle_ety()/GeV << endl;
    cout << "Egamma10NoTau_RefEle_sumet             " << d3pd.met.Egamma10NoTau_RefEle_sumet()/GeV << endl;
    cout << "Egamma10NoTau_Muon_Total_Staco_etx     " << d3pd.met.Egamma10NoTau_Muon_Total_Staco_etx()/GeV << endl;
    cout << "Egamma10NoTau_Muon_Total_Staco_ety     " << d3pd.met.Egamma10NoTau_Muon_Total_Staco_ety()/GeV << endl;
    cout << "Egamma10NoTau_Muon_Total_Staco_sumet   " << d3pd.met.Egamma10NoTau_Muon_Total_Staco_sumet()/GeV << endl;
    cout << "Egamma10NoTau_RefJet_etx               " << d3pd.met.Egamma10NoTau_RefJet_etx()/GeV << endl;
    cout << "Egamma10NoTau_RefJet_ety               " << d3pd.met.Egamma10NoTau_RefJet_ety()/GeV << endl;
    cout << "Egamma10NoTau_RefJet_sumet             " << d3pd.met.Egamma10NoTau_RefJet_sumet()/GeV << endl;
    cout << "Egamma10NoTau_CellOut_Eflow_STVF_etx   " << d3pd.met.Egamma10NoTau_CellOut_Eflow_STVF_etx()/GeV << endl;
    cout << "Egamma10NoTau_CellOut_Eflow_STVF_ety   " << d3pd.met.Egamma10NoTau_CellOut_Eflow_STVF_ety()/GeV << endl;
    cout << "Egamma10NoTau_CellOut_Eflow_STVF_sumet " << d3pd.met.Egamma10NoTau_CellOut_Eflow_STVF_sumet()/GeV << endl;
    cout << "Egamma10NoTau_SoftJets_etx             " << d3pd.met.Egamma10NoTau_SoftJets_etx()/GeV << endl;
    cout << "Egamma10NoTau_SoftJets_ety             " << d3pd.met.Egamma10NoTau_SoftJets_ety()/GeV << endl;
    cout << "Egamma10NoTau_SoftJets_sumet           " << d3pd.met.Egamma10NoTau_SoftJets_sumet()/GeV << endl;

    // Dump muon weights
    cout.precision(2);
    for(int iMu=0; iMu<d3pd.muo.n(); iMu++){
      const MuonElement* mu = & d3pd.muo[iMu];
      const TLorentzVector* lv = & m_susyObj.GetMuonTLV(iMu);
      cout << "  Mu : " << fixed
           << " q " << setw(2) << (int) mu->charge()
           << " pt " << setw(6) << lv->Pt()/GeV
           << " eta " << setw(5) << lv->Eta()
           << endl;
      // Dump the met weights
      //cout << "       stvf weights: " << mu->MET_Egamma10NoTau_wet().size() << endl;
      for(uint iW = 0; iW < mu->MET_Egamma10NoTau_wet().size(); iW++)
      {
        cout << "       stvf weights:" 
             << " wet " << mu->MET_Egamma10NoTau_wet()[iW] 
             << " wpx " << mu->MET_Egamma10NoTau_wpx()[iW] 
             << " wpy " << mu->MET_Egamma10NoTau_wpy()[iW]
             << endl;
      }
      //cout << "       reff weights: " << mu->MET_Egamma10NoTau_wet().size() << endl;
      for(uint iW = 0; iW < mu->MET_Egamma10NoTau_wet().size(); iW++)
      {
        cout << "       reff weights:"
             << " wet " << mu->MET_Egamma10NoTau_wet()[iW] 
             << " wpx " << mu->MET_Egamma10NoTau_wpx()[iW] 
             << " wpy " << mu->MET_Egamma10NoTau_wpy()[iW]
             << endl;
      }
    }
    cout.precision(6);
    cout.unsetf(ios_base::fixed);
  }

  // Fill histograms here

  // Fill MET histograms by looping over met indices
  for(uint i=0; i<Met_N; i++){
    h_met[i]->Fill(met[i].Mod()/GeV, w);
    h_metEle[i]->Fill(metEle[i].Mod()/GeV, w);
    h_metMuo[i]->Fill(metMuo[i].Mod()/GeV, w);
    h_metJet[i]->Fill(metJet[i].Mod()/GeV, w);
    h_metCell[i]->Fill(metCell[i].Mod()/GeV, w);
  }

  // Fill muon MET variables

  // susy_stvf met
  // Use the preMuons for the susy_stvf met, and all weights are 1
  for(uint i=0; i<m_preMuons.size(); i++){

    const TLorentzVector* muLV = & m_susyObj.GetMuonTLV(m_preMuons[i]);

    // muon met weights - I think these three are identical
    h_muWet[Met_susy_stvf]->Fill(1, w);
    h_muWpx[Met_susy_stvf]->Fill(1, w);
    h_muWpy[Met_susy_stvf]->Fill(1, w);

    // muon kinematics
    h_metMuPt[Met_susy_stvf]->Fill(muLV->Pt()/GeV, w);
    h_metLepPt[Met_susy_stvf]->Fill(muLV->Pt()/GeV, w);
    h_metMuEta[Met_susy_stvf]->Fill(muLV->Eta(), w);
    h_metLepEta[Met_susy_stvf]->Fill(muLV->Eta(), w);

  }

  h_nMetMu[Met_susy_stvf]->Fill(m_preMuons.size(), w);

  // Use d3pd muons for the d3pd_* met, with d3pd weights
  // Each weight is a vector, but perhaps only one entry?
  // For now, just use the first entry.  Assume 1 entry.
  uint nMuStvf = 0;
  uint nMuReff = 0;
  for(int i=0; i<d3pd.muo.n(); i++){

    const MuonElement* mu = & d3pd.muo[i];
    //const TLorentzVector* muLV = & m_susyObj.GetMuonTLV(i);

    // muon met weights
    h_muWet[Met_d3pd_stvf]->Fill(mu->MET_Egamma10NoTau_wet().at(0), w);
    h_muWpx[Met_d3pd_stvf]->Fill(mu->MET_Egamma10NoTau_wpx().at(0), w);
    h_muWpy[Met_d3pd_stvf]->Fill(mu->MET_Egamma10NoTau_wpy().at(0), w);
    h_muWet[Met_d3pd_reff]->Fill(mu->MET_Egamma10NoTau_wet().at(0), w);
    h_muWpx[Met_d3pd_reff]->Fill(mu->MET_Egamma10NoTau_wpx().at(0), w);
    h_muWpy[Met_d3pd_reff]->Fill(mu->MET_Egamma10NoTau_wpy().at(0), w);

    // muon kinematics
    if(mu->MET_Egamma10NoTau_wet().at(0) != 0){
      nMuStvf++;
      h_metMuPt[Met_d3pd_stvf]->Fill(mu->pt()/GeV, w);
      h_metLepPt[Met_d3pd_stvf]->Fill(mu->pt()/GeV, w);
      h_metMuEta[Met_d3pd_stvf]->Fill(mu->eta(), w);
      h_metLepEta[Met_d3pd_stvf]->Fill(mu->eta(), w);
    }
    if(mu->MET_Egamma10NoTau_wet().at(0) != 0){
      nMuReff++;
      h_metMuPt[Met_d3pd_reff]->Fill(mu->pt()/GeV, w);
      h_metLepPt[Met_d3pd_reff]->Fill(mu->pt()/GeV, w);
      h_metMuEta[Met_d3pd_reff]->Fill(mu->eta(), w);
      h_metLepEta[Met_d3pd_reff]->Fill(mu->eta(), w);
    }

    // raw muon histos
    h_rawMuPt->Fill(mu->pt()/GeV, w);
    h_rawMuEta->Fill(mu->eta(), w);

  }
  h_nMetMu[Met_d3pd_stvf]->Fill(nMuStvf, w);
  h_nMetMu[Met_d3pd_reff]->Fill(nMuReff, w);

  // Fill electron MET variables
  // TODO

  // Fill other histos
  h_nMu->Fill(nSigMuo, w);

  // Raw muons

  return true;
}

/*--------------------------------------------------------------------------------*/
// Trigger cut
/*--------------------------------------------------------------------------------*/
bool SusyMetValidation::passTrigger()
{
  // Vectors of the signal leptons
  static vector<TLorentzVector> elLVs;
  static vector<TLorentzVector> muLVs;
  elLVs.clear();
  muLVs.clear();
  for(uint iEl=0; iEl<m_sigElectrons.size(); iEl++){
    elLVs.push_back(m_susyObj.GetElecTLV(m_sigElectrons[iEl]));
  }
  for(uint iMu=0; iMu<m_sigMuons.size(); iMu++){
    muLVs.push_back(m_susyObj.GetMuonTLV(m_sigMuons[iMu]));
  }
  string stream = streamName(m_stream);
  if(stream == "Egamma") stream = "EGamma";
  return m_triggerMatch->passesMultiLepTrigger(&elLVs, &muLVs, m_sigElectrons, m_sigMuons,
                                               !m_isMC, d3pd.evt.RunNumber(), stream,
                                               d3pd.trig.trig_EF_el_n(), d3pd.trig.trig_EF_el_px(), 
                                               d3pd.trig.trig_EF_el_py(), 
                                               d3pd.trig.trig_EF_el_pz(), d3pd.trig.trig_EF_el_E(),
                                               d3pd.trig.trig_EF_el_EF_e24vh_medium1(), 
                                               d3pd.trig.trig_EF_el_EF_e24vh_medium1_e7_medium1(),
                                               d3pd.trig.trig_EF_el_EF_2e12Tvh_loose1(), 
                                               d3pd.trig.trig_EF_el_EF_e7T_medium1(),
                                               d3pd.trig.trig_EF_el_EF_e12Tvh_medium1_mu8(), 
                                               d3pd.trig.trig_EF_trigmuonef_n(),
                                               d3pd.trig.trig_EF_trigmuonef_track_CB_eta(),
                                               d3pd.trig.trig_EF_trigmuonef_track_CB_phi(),
                                               d3pd.trig.trig_EF_trigmuonef_track_CB_hasCB(),
                                               d3pd.trig.trig_EF_trigmuonef_EF_mu18_tight(),
                                               d3pd.trig.trig_EF_trigmuonef_EF_mu18_tight_mu8_EFFS(),
                                               d3pd.trig.trig_EF_trigmuonef_EF_2mu13(), 
                                               d3pd.trig.trig_EF_trigmuonef_EF_mu8(),
                                               d3pd.trig.trig_EF_trigmuonef_EF_mu18_tight_e7_medium1(),
                                               d3pd.trig.EF_e24vh_medium1_e7_medium1(), 
                                               d3pd.trig.EF_mu18_tight_mu8_EFFS(),
                                               d3pd.trig.EF_2mu13(), d3pd.trig.EF_2e12Tvh_loose1(),
                                               d3pd.trig.EF_e12Tvh_medium1_mu8(), 
                                               d3pd.trig.EF_mu18_tight_e7_medium1(), false);
}


#undef GeV
