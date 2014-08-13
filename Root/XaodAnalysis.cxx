#include <limits>

#include "TSystem.h"

#include "SusyCommon/XaodAnalysis.h"
#include "SusyCommon/get_object_functions.h"
#include "egammaAnalysisUtils/egammaTriggerMatching.h"
#include "D3PDReader/JetD3PDObject.h"

using namespace std;
using susy::XaodAnalysis;

#define GeV 1000.

/*--------------------------------------------------------------------------------*/
// XaodAnalysis Constructor
/*--------------------------------------------------------------------------------*/
XaodAnalysis::XaodAnalysis() :
        m_sample(""),
        m_stream(Stream_Unknown),
        m_isAF2(false),
        m_mcProd(MCProd_Unknown),
        m_d3pdTag(D3PD_p1328),
        m_selectPhotons(false),
        m_selectTaus(false),
        m_selectTruth(false),
        m_metFlavor(SUSYMet::Default),
        m_doMetMuCorr(false),
        m_doMetFix(false),
        m_lumi(LUMI_A_E),
        m_sumw(1),
        m_xsec(-1),
        m_errXsec(-1),
        m_mcRun(0),
        m_mcLB(0),
        m_sys(false),
        m_eleMediumSFTool(0),
        m_pileup(0),
        m_pileup_up(0),
        m_pileup_dn(0),
        m_susyXsec(0),
        m_hforTool(),
        m_tree(0),
        m_entry(0),
        m_dbg(0),
        m_isMC(false),
        m_flagsAreConsistent(false),
        m_flagsHaveBeenChecked(false)
{
}
/*--------------------------------------------------------------------------------*/
// Destructor
/*--------------------------------------------------------------------------------*/
XaodAnalysis::~XaodAnalysis()
{
    cout<<"~XaodAnalysis : todo"<<endl;
}
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::SlaveBegin(TTree *tree)
{
  if(m_dbg) cout << "XaodAnalysis::SlaveBegin" << endl;

  bool isData = m_sample.Contains("data", TString::kIgnoreCase);
  m_isMC = !isData;

  // Make sure MC production is specified
  if(m_isMC && m_mcProd==MCProd_Unknown){
    cout << "XaodAnalysis::SlaveBegin : ERROR : Sample is flagged as MC but "
         << "MCProduction is Unknown! Use command line argument to set it!"
         << endl;
    abort();
  }

  // Use sample name to set data stream
  if(m_isMC) m_stream = Stream_MC;
  else if(m_sample.Contains("muons", TString::kIgnoreCase)) m_stream = Stream_Muons;
  else if(m_sample.Contains("egamma", TString::kIgnoreCase)) m_stream = Stream_Egamma;
  else if(m_sample.Contains("jettauetmiss", TString::kIgnoreCase)) m_stream = Stream_JetTauEtmiss;
  else m_stream = Stream_Unknown;

  if(m_isMC) cout << "Processing as MC"   << endl;
  else       cout << "Processing as DATA" << endl;
  cout << "DataStream: " << streamName(m_stream) << endl;

  // Setup SUSYTools
  bool isMC12b = (m_mcProd == MCProd_MC12b);
  bool useLeptonTrigger = false;
  if(true /*m_dbg*/) objTool.msg().setLevel( MSG::DEBUG);
  objTool.setProperty("IsData",isData);
  objTool.setProperty("IsAtlfast", m_isAF2);   
  objTool.setProperty("IsMC12b", isMC12b);
  objTool.setProperty("UseLeptonTrigger",useLeptonTrigger);
  if( objTool.initialize() != StatusCode::SUCCESS){
      Error( APP_NAME, "Cannot intialize SUSYObjDef_xAOD..." );
      Error( APP_NAME, "Exiting... " );
      exit(-1);
  }else{
      Info( APP_NAME, "SUSYObjDef_xAOD initialized... " );
  }

#warning TElectronEfficiencyCorrectionTool not initialized
#warning fakemet_est tool not initialized
  if(isData){ initGrlTool(); }
  if(m_isMC){
#warning susy xsec tool not initialized
#warning pileup rew tool not initialized
  }
}


/*--------------------------------------------------------------------------------*/
// Main process loop function - This is just an example for testing
/*--------------------------------------------------------------------------------*/
Bool_t XaodAnalysis::Process(Long64_t entry)
{
  m_event.GetEntry(entry);

  static Long64_t chainEntry = -1;
  chainEntry++;
  if(m_dbg || chainEntry%10000==0)
  {
    cout << "**** Processing entry " << setw(6) << chainEntry
         << " run " << setw(6) << m_event.eventinfo.RunNumber()
         << " event " << setw(7) << m_event.eventinfo.EventNumber() << " ****" << endl;
  }

  // Object selection
  m_susyObj.Reset();
  clearObjects();
  selectObjects();
  buildMet();

  return kTRUE;
}

/*--------------------------------------------------------------------------------*/
// The Terminate() function is the last function to be called during
// a query. It always runs on the client, it can be used to present
// the results graphically or save the results to file.
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::Terminate()
{
  if(m_dbg) cout << "XaodAnalysis::Terminate" << endl;
  m_susyObj.finalize();

  if(m_isMC){
    delete m_susyXsec;
    delete m_pileup;
    delete m_pileup_up;
    delete m_pileup_dn;
  }
}
//----------------------------------------------------------
D3PDReader::MuonD3PDObject* XaodAnalysis::d3pdMuons()
{
    return &m_event.mu_staco;
}
//----------------------------------------------------------
D3PDReader::ElectronD3PDObject* XaodAnalysis::d3pdElectrons()
{
    return &m_event.el;
}
//----------------------------------------------------------
D3PDReader::TauD3PDObject* XaodAnalysis::d3pdTaus()
{
    return &m_event.tau;
}
//----------------------------------------------------------
D3PDReader::JetD3PDObject* XaodAnalysis::d3pdJets()
{
    return &m_event.jet_AntiKt4LCTopo;
}
/*--------------------------------------------------------------------------------*/
// Baseline object selection
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::selectBaselineObjects(SusyNtSys sys)
{
  if(m_dbg>=5) cout << "selectBaselineObjects" << endl;
  vector<int> goodJets;  // What the hell is this??

  // SUSYTools takes a flag for AF2. Now set via command line flag
  //bool isAF2 = false;

  // Handle Systematic
  // New syntax for SUSYTools in mc12
  SystErr::Syste susySys = SystErr::NONE;
  if(sys == NtSys_NOM);                                           // No need to check needlessly
  //else if(sys == NtSys_EES_UP) susySys = SystErr::EESUP;        // E scale up
  //else if(sys == NtSys_EES_DN) susySys = SystErr::EESDOWN;      // E scale down
  else if(sys == NtSys_EES_Z_UP  ) susySys = SystErr::EGZEEUP;    // E scale Zee up
  else if(sys == NtSys_EES_Z_DN  ) susySys = SystErr::EGZEEDOWN;  // E scale Zee dn
  else if(sys == NtSys_EES_MAT_UP) susySys = SystErr::EGMATUP;    // E scale material up
  else if(sys == NtSys_EES_MAT_DN) susySys = SystErr::EGMATDOWN;  // E scale material down
  else if(sys == NtSys_EES_PS_UP ) susySys = SystErr::EGPSUP;     // E scale presampler up
  else if(sys == NtSys_EES_PS_DN ) susySys = SystErr::EGPSDOWN;   // E scale presampler down
  else if(sys == NtSys_EES_LOW_UP) susySys = SystErr::EGLOWUP;    // E low pt up
  else if(sys == NtSys_EES_LOW_DN) susySys = SystErr::EGLOWDOWN;  // E low pt down
  else if(sys == NtSys_EER_UP    ) susySys = SystErr::EGRESUP;    // E smear up
  else if(sys == NtSys_EER_DN    ) susySys = SystErr::EGRESDOWN;  // E smear down
  else if(sys == NtSys_MS_UP     ) susySys = SystErr::MMSUP;      // MS scale up
  else if(sys == NtSys_MS_DN     ) susySys = SystErr::MMSLOW;     // MS scale down
  else if(sys == NtSys_ID_UP     ) susySys = SystErr::MIDUP;      // ID scale up
  else if(sys == NtSys_ID_DN     ) susySys = SystErr::MIDLOW;     // ID scale down
  else if(sys == NtSys_JES_UP    ) susySys = SystErr::JESUP;      // JES up
  else if(sys == NtSys_JES_DN    ) susySys = SystErr::JESDOWN;    // JES down
  else if(sys == NtSys_JER       ) susySys = SystErr::JER;        // JER (gaussian)

  else if(sys == NtSys_TES_UP    ) susySys = SystErr::TESUP;      // TES up
  else if(sys == NtSys_TES_DN    ) susySys = SystErr::TESDOWN;    // TES down

  D3PDReader::JetD3PDObject *jets = d3pdJets();
  // Container object selection
  if(m_selectTaus) m_contTaus = get_taus_baseline(d3pdTaus(), m_susyObj, 20.*GeV, 2.47,
                                                  SUSYTau::TauNone, SUSYTau::TauNone, SUSYTau::TauNone,
                                                  susySys, true);

  // Preselection
  m_preElectrons = get_electrons_baseline(d3pdElectrons(), &m_event.el_MET_Egamma10NoTau,
                                          !m_isMC, m_event.eventinfo.RunNumber(), m_susyObj,
                                          7.*GeV, 2.47, susySys);
  m_preMuons = get_muons_baseline(d3pdMuons(), !m_isMC, m_susyObj,
                                  6.*GeV, 2.5, susySys);
  // Removing eta cut for baseline jets. This is for the bad jet veto.
  m_preJets = get_jet_baseline(jets, &m_event.vxp, &m_event.eventinfo, &m_event.Eventshape, !m_isMC, m_susyObj,
                               20.*GeV, std::numeric_limits<float>::max(), susySys, false, goodJets);
  //m_preJets = get_jet_baseline(jets, &m_event.vxp, &m_event.eventinfo, !m_isMC, m_susyObj,
  //                             20.*GeV, 4.9, susySys, false, goodJets);

  // Selection for met muons
  // Diff with preMuons is pt selection
  m_metMuons = get_muons_baseline(d3pdMuons(), !m_isMC, m_susyObj,
                                  10.*GeV, 2.5, susySys);

  // Preselect taus
  if(m_selectTaus) m_preTaus = get_taus_baseline(d3pdTaus(), m_susyObj, 20.*GeV, 2.47,
                                                 SUSYTau::TauLoose, SUSYTau::TauLoose, SUSYTau::TauLoose,
                                                 susySys, true);

  performOverlapRemoval();

  // combine leptons
  m_preLeptons    = buildLeptonInfos(d3pdElectrons(), m_preElectrons, d3pdMuons(), m_preMuons, m_susyObj);
  m_baseLeptons   = buildLeptonInfos(d3pdElectrons(), m_baseElectrons, d3pdMuons(), m_baseMuons, m_susyObj);
}

/*--------------------------------------------------------------------------------*/
// perform overlap
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::performOverlapRemoval()
{
  D3PDReader::JetD3PDObject *jets = d3pdJets();
  // e-e overlap removal
  m_baseElectrons = overlap_removal(m_susyObj, d3pdElectrons(), m_preElectrons, d3pdElectrons(), m_preElectrons,
                                    0.05, true, true);
  // jet-e overlap removal
  m_baseJets      = overlap_removal(m_susyObj, jets, m_preJets, d3pdElectrons(), m_baseElectrons,
                                    0.2, false, false);

  if(m_selectTaus) {
    // tau-e overlap removal
    m_baseTaus    = overlap_removal(m_susyObj, d3pdTaus(), m_preTaus, d3pdElectrons(), m_baseElectrons, 0.2, false, false);
    // tau-mu overlap removal
    m_baseTaus    = overlap_removal(m_susyObj, d3pdTaus(), m_baseTaus, d3pdMuons(), m_preMuons, 0.2, false, false);
  }

  // e-jet overlap removal
  m_baseElectrons = overlap_removal(m_susyObj, d3pdElectrons(), m_baseElectrons, jets, m_baseJets,
                                    0.4, false, false);

  // m-jet overlap removal
  m_baseMuons     = overlap_removal(m_susyObj, d3pdMuons(), m_preMuons, jets, m_baseJets, 0.4, false, false);

  // e-m overlap removal
  vector<int> copyElectrons = m_baseElectrons;
  m_baseElectrons = overlap_removal(m_susyObj, d3pdElectrons(), m_baseElectrons, d3pdMuons(), m_baseMuons,
                                    0.01, false, false);
  m_baseMuons     = overlap_removal(m_susyObj, d3pdMuons(), m_baseMuons, d3pdElectrons(), copyElectrons, 0.01, false, false);

  // m-m overlap removal
  m_baseMuons     = overlap_removal(m_susyObj, d3pdMuons(), m_baseMuons, d3pdMuons(), m_baseMuons, 0.05, true, false);

  // jet-tau overlap removal
  m_baseJets      = overlap_removal(m_susyObj, jets, m_baseJets, d3pdTaus(), m_baseTaus, 0.2, false, false);

  // remove SFOS lepton pairs with Mll < 12 GeV
  m_baseElectrons = RemoveSFOSPair(m_susyObj, d3pdElectrons(), m_baseElectrons, 12.*GeV);
  m_baseMuons     = RemoveSFOSPair(m_susyObj, d3pdMuons(), m_baseMuons,     12.*GeV);
  //m_baseTaus      = RemoveSFOSPair(m_susyObj, d3pdTaus(), m_baseTaus,      12.*GeV);
}

/*--------------------------------------------------------------------------------*/
// Signal object selection - do baseline selection first!
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::selectSignalObjects()
{
  if(m_dbg>=5) cout << "selectSignalObjects" << endl;
  uint nVtx = getNumGoodVtx();
  D3PDReader::JetD3PDObject *jets =  d3pdJets();
  m_sigElectrons = get_electrons_signal(d3pdElectrons(), m_baseElectrons, d3pdMuons(), m_baseMuons,
                                        nVtx, !m_isMC, m_susyObj, 10.*GeV, 0.16, 0.18, 5., 0.4);
  m_sigMuons     = get_muons_signal(d3pdMuons(), m_baseMuons, d3pdElectrons(), m_baseElectrons,
                                    nVtx, !m_isMC, m_susyObj, 10.*GeV, .12, 3., 1.);
  m_sigJets      = get_jet_signal(jets, m_susyObj, m_baseJets, 20.*GeV, 2.5, 0.75);
  m_sigTaus      = get_taus_signal(d3pdTaus(), m_baseTaus, m_susyObj);

  // combine light leptons
  m_sigLeptons   = buildLeptonInfos(d3pdElectrons(), m_sigElectrons, d3pdMuons(), m_sigMuons, m_susyObj);

  // photon selection done in separate method, why?
  if(m_selectPhotons) selectSignalPhotons();
}

/*--------------------------------------------------------------------------------*/
// Build MissingEt
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::buildMet(SusyNtSys sys)
{
  if(m_dbg>=5) cout << "buildMet" << endl;

  // Need the proper jet systematic for building systematic
  SystErr::Syste susySys = SystErr::NONE;
  if(sys == NtSys_NOM);
  else if(sys == NtSys_JES_UP)      susySys = SystErr::JESUP;       // JES up
  else if(sys == NtSys_JES_DN)      susySys = SystErr::JESDOWN;     // JES down
  else if(sys == NtSys_JER)         susySys = SystErr::JER;         // JER (gaussian)
  else if(sys == NtSys_SCALEST_UP)  susySys = SystErr::SCALESTUP;   // Met scale sys up
  else if(sys == NtSys_SCALEST_DN)  susySys = SystErr::SCALESTDOWN; // Met scale sys down
  // Only one of these now?
  //else if(sys == NtSys_RESOST_UP)   susySys = SystErr::RESOSTUP;    // Met resolution sys up
  //else if(sys == NtSys_RESOST_DN)   susySys = SystErr::RESOSTDOWN;  // Met resolution sys down
  else if(sys == NtSys_RESOST)      susySys = SystErr::RESOST;      // Met resolution sys up

  // Need electrons with nonzero met weight in order to calculate the MET
  vector<int> metElectrons = get_electrons_met(&m_event.el_MET_Egamma10NoTau, m_susyObj);

  // Calculate the MET
  // We use the metMuons instead of preMuons so that we can have a lower pt cut on preMuons
  TVector2 metVector =  m_susyObj.GetMET(m_event.jet_AntiKt4LCTopo_MET_Egamma10NoTau.wet(), m_event.jet_AntiKt4LCTopo_MET_Egamma10NoTau.wpx(),
                                         m_event.jet_AntiKt4LCTopo_MET_Egamma10NoTau.wpy(), m_event.jet_AntiKt4LCTopo_MET_Egamma10NoTau.statusWord(),
                                         metElectrons,
                                         m_event.el_MET_Egamma10NoTau.wet(), m_event.el_MET_Egamma10NoTau.wpx(),
                                         m_event.el_MET_Egamma10NoTau.wpy(), m_event.el_MET_Egamma10NoTau.statusWord(),
                                         m_event.MET_CellOut_Egamma10NoTau.etx(),
                                         m_event.MET_CellOut_Egamma10NoTau.ety(),
                                         m_event.MET_CellOut_Egamma10NoTau.sumet(),
                                         m_event.MET_CellOut_Eflow_STVF_Egamma10NoTau.etx(),
                                         m_event.MET_CellOut_Eflow_STVF_Egamma10NoTau.ety(),
                                         m_event.MET_CellOut_Eflow_STVF_Egamma10NoTau.sumet(),
                                         m_event.MET_RefGamma_Egamma10NoTau.etx(),
                                         m_event.MET_RefGamma_Egamma10NoTau.ety(),
                                         m_event.MET_RefGamma_Egamma10NoTau.sumet(),
                                         m_metMuons,
                                         d3pdMuons()->ms_qoverp(),
                                         d3pdMuons()->ms_theta(),
                                         d3pdMuons()->ms_phi(),
                                         d3pdMuons()->charge(),
                                         d3pdMuons()->energyLossPar(),
                                         m_event.eventinfo.averageIntPerXing(),
                                         m_metFlavor, susySys);
  m_met.SetPxPyPzE(metVector.X(), metVector.Y(), 0, metVector.Mod());
}

/*--------------------------------------------------------------------------------*/
// Signal photons
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::selectSignalPhotons()
{
  if(m_dbg>=5) cout << "selectSignalPhotons" << endl;

  int phoQual = 2;      // Quality::Tight
  uint isoType = 1;     // Corresponds to PTED corrected isolation
  float etcone40CorrCut = 3*GeV;

  vector<int> base_photons = get_photons_baseline(&m_event.ph, m_susyObj,
                                                  20.*GeV, 2.47, SystErr::NONE, phoQual);

  // Latest and Greatest
  int nPV = getNumGoodVtx();
  m_sigPhotons = get_photons_signal(&m_event.ph, base_photons, m_susyObj, nPV,
                                    20.*GeV, etcone40CorrCut, isoType);
}
/*--------------------------------------------------------------------------------*/
// Truth object selection
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::selectTruthObjects()
{
  if(m_dbg>=5) cout << "selectTruthObjects" << endl;

  // ==>> First the truth particles
  // Done under SusyNtMaker::fillTruthParticleVars

  // ==>> Second the truth jets
  for(int index=0; index < m_event.AntiKt4Truth.n(); index++) {
      const D3PDReader::JetD3PDObjectElement &trueJet = m_event.AntiKt4Truth[index];
      if( trueJet.pt()/GeV > 15. && fabs(trueJet.eta()) < 4.5) m_truJets.push_back(index);
  }

  // ==>> Third and last the truth met
  m_truMet.SetPxPyPzE(m_event.MET_Truth.NonInt_etx(), m_event.MET_Truth.NonInt_ety(), 0, m_event.MET_Truth.NonInt_sumet());
}

/*--------------------------------------------------------------------------------*/
// Clear selected objects
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::clearObjects()
{
  m_preElectrons.clear();
  m_preMuons.clear();
  m_preJets.clear();
  m_preLeptons.clear();
  m_baseElectrons.clear();
  m_baseMuons.clear();
  m_baseLeptons.clear();
  m_baseJets.clear();
  m_sigElectrons.clear();
  m_sigMuons.clear();
  m_sigLeptons.clear();
  m_sigJets.clear();
  m_cutFlags = 0;

  m_sigPhotons.clear();
  m_truParticles.clear();
  m_truJets.clear();

  m_metMuons.clear();
}

/*--------------------------------------------------------------------------------*/
// Count number of good vertices
/*--------------------------------------------------------------------------------*/
uint XaodAnalysis::getNumGoodVtx()
{
  uint nVtx = 0;
  for(int i=0; i < m_event.vxp.n(); i++){
    if(m_event.vxp.nTracks()->at(i) >= 5) nVtx++;
  }
  return nVtx;
}

/*--------------------------------------------------------------------------------*/
// Match reco jet to a truth jet
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::matchTruthJet(int iJet)
{
  // Loop over truth jets looking for a match
  const TLorentzVector &jetLV = m_susyObj.GetJetTLV(iJet);
  for(int i=0; i<m_event.AntiKt4Truth.n(); i++){
    const D3PDReader::JetD3PDObjectElement &trueJet = m_event.AntiKt4Truth[i];
    TLorentzVector trueJetLV;
    trueJetLV.SetPtEtaPhiE(trueJet.pt(), trueJet.eta(), trueJet.phi(), trueJet.E());
    if(jetLV.DeltaR(trueJetLV) < 0.3) return true;
  }
  return false;
}

/*--------------------------------------------------------------------------------*/
// Event trigger flags
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::fillEventTriggers()
{
  if(m_dbg>=5) cout << "fillEventTriggers" << endl;

  m_evtTrigFlags = 0;
  // e7_medium1 not available at the moment, so use e7T for now
  //if(m_event.triggerbits.EF_e7_medium1())                 m_evtTrigFlags |= TRIG_e7_medium1;
  if(m_event.triggerbits.EF_e7T_medium1())                m_evtTrigFlags |= TRIG_e7_medium1;
  if(m_event.triggerbits.EF_e12Tvh_loose1())              m_evtTrigFlags |= TRIG_e12Tvh_loose1;
  if(m_event.triggerbits.EF_e12Tvh_medium1())             m_evtTrigFlags |= TRIG_e12Tvh_medium1;
  if(m_event.triggerbits.EF_e24vh_medium1())              m_evtTrigFlags |= TRIG_e24vh_medium1;
  if(m_event.triggerbits.EF_e24vhi_medium1())             m_evtTrigFlags |= TRIG_e24vhi_medium1;
  if(m_event.triggerbits.EF_2e12Tvh_loose1())             m_evtTrigFlags |= TRIG_2e12Tvh_loose1;
  if(m_event.triggerbits.EF_e24vh_medium1_e7_medium1())   m_evtTrigFlags |= TRIG_e24vh_medium1_e7_medium1;
  if(m_event.triggerbits.EF_mu8())                        m_evtTrigFlags |= TRIG_mu8;
  if(m_event.triggerbits.EF_mu13())                       m_evtTrigFlags |= TRIG_mu13;
  if(m_event.triggerbits.EF_mu18_tight())                 m_evtTrigFlags |= TRIG_mu18_tight;
  if(m_event.triggerbits.EF_mu24i_tight())                m_evtTrigFlags |= TRIG_mu24i_tight;
  if(m_event.triggerbits.EF_2mu13())                      m_evtTrigFlags |= TRIG_2mu13;
  if(m_event.triggerbits.EF_mu18_tight_mu8_EFFS())        m_evtTrigFlags |= TRIG_mu18_tight_mu8_EFFS;
  if(m_event.triggerbits.EF_e12Tvh_medium1_mu8())         m_evtTrigFlags |= TRIG_e12Tvh_medium1_mu8;
  if(m_event.triggerbits.EF_mu18_tight_e7_medium1())      m_evtTrigFlags |= TRIG_mu18_tight_e7_medium1;

  if(m_event.triggerbits.EF_tau20_medium1())                   m_evtTrigFlags |= TRIG_tau20_medium1;
  if(m_event.triggerbits.EF_tau20Ti_medium1())                 m_evtTrigFlags |= TRIG_tau20Ti_medium1;
  if(m_event.triggerbits.EF_tau29Ti_medium1())                 m_evtTrigFlags |= TRIG_tau29Ti_medium1;
  if(m_event.triggerbits.EF_tau29Ti_medium1_tau20Ti_medium1()) m_evtTrigFlags |= TRIG_tau29Ti_medium1_tau20Ti_medium1;
  if(m_event.triggerbits.EF_tau20Ti_medium1_e18vh_medium1())   m_evtTrigFlags |= TRIG_tau20Ti_medium1_e18vh_medium1;
  if(m_event.triggerbits.EF_tau20_medium1_mu15())              m_evtTrigFlags |= TRIG_tau20_medium1_mu15;

  if(m_event.triggerbits.EF_e18vh_medium1())              m_evtTrigFlags |= TRIG_e18vh_medium1;
  if(m_event.triggerbits.EF_mu15())                       m_evtTrigFlags |= TRIG_mu15;

  // EF_2mu8_EFxe40wMu_tclcw trigger only available for data, in periods > B
  if(!m_isMC && m_event.eventinfo.RunNumber()>=206248 && m_event.triggerbits.EF_2mu8_EFxe40wMu_tclcw())
    m_evtTrigFlags |= TRIG_2mu8_EFxe40wMu_tclcw;

  // Triggers requested fro the ISR analysis studies
  if(m_event.triggerbits.EF_mu6())                                m_evtTrigFlags |= TRIG_mu6;
  if(m_event.triggerbits.EF_2mu6())                               m_evtTrigFlags |= TRIG_2mu6;
  if(m_event.triggerbits.EF_e18vh_medium1_2e7T_medium1())         m_evtTrigFlags |= TRIG_e18vh_medium1_2e7T_medium1;
  if(m_event.triggerbits.EF_3mu6())                               m_evtTrigFlags |= TRIG_3mu6;
  if(m_event.triggerbits.EF_mu18_tight_2mu4_EFFS())               m_evtTrigFlags |= TRIG_mu18_tight_2mu4_EFFS;
  if(m_event.triggerbits.EF_2e7T_medium1_mu6())                   m_evtTrigFlags |= TRIG_2e7T_medium1_mu6;
  if(m_event.triggerbits.EF_e7T_medium1_2mu6())                   m_evtTrigFlags |= TRIG_e7T_medium1_2mu6;
  if(m_event.triggerbits.EF_xe80_tclcw_loose())                   m_evtTrigFlags |= TRIG_xe80_tclcw_loose;
  if(m_event.triggerbits.EF_j110_a4tchad_xe90_tclcw_loose())      m_evtTrigFlags |= TRIG_j110_a4tchad_xe90_tclcw_loose;
  if(m_event.triggerbits.EF_j80_a4tchad_xe100_tclcw_loose())      m_evtTrigFlags |= TRIG_j80_a4tchad_xe100_tclcw_loose;
  if(m_event.triggerbits.EF_j80_a4tchad_xe70_tclcw_dphi2j45xe10())m_evtTrigFlags |= TRIG_j80_a4tchad_xe70_tclcw_dphi2j45xe10;

  // Not sure about the availability of these, so just adding some protection
  if(m_event.triggerbits.EF_mu4T())                               m_evtTrigFlags |= TRIG_mu4T;
  if(m_event.triggerbits.EF_mu24())                               m_evtTrigFlags |= TRIG_mu24;
  if(m_event.triggerbits.EF_mu4T_j65_a4tchad_xe70_tclcw_veryloose()) m_evtTrigFlags |= TRIG_mu4T_j65_a4tchad_xe70_tclcw_veryloose;
  if(m_event.triggerbits.EF_2mu4T_xe60_tclcw())                   m_evtTrigFlags |= TRIG_2mu4T_xe60_tclcw;
  if(m_event.triggerbits.EF_2mu8_EFxe40_tclcw.IsAvailable() && m_event.triggerbits.EF_2mu8_EFxe40_tclcw())
    m_evtTrigFlags |= TRIG_2mu8_EFxe40_tclcw;
  if(m_event.triggerbits.EF_e24vh_medium1_EFxe35_tclcw())         m_evtTrigFlags |= TRIG_e24vh_medium1_EFxe35_tclcw;
  if(m_event.triggerbits.EF_mu24_j65_a4tchad_EFxe40_tclcw())      m_evtTrigFlags |= TRIG_mu24_j65_a4tchad_EFxe40_tclcw;
  if(m_event.triggerbits.EF_mu24_j65_a4tchad_EFxe40wMu_tclcw.IsAvailable() && m_event.triggerbits.EF_mu24_j65_a4tchad_EFxe40wMu_tclcw())
    m_evtTrigFlags |= TRIG_mu24_j65_a4tchad_EFxe40wMu_tclcw;
}

/*--------------------------------------------------------------------------------*/
// Electron trigger matching
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::matchElectronTriggers()
{
  if(m_dbg>=5) cout << "matchElectronTriggers" << endl;
  //int run = m_event.eventinfo.RunNumber();

  // loop over all pre electrons
  for(uint i=0; i<m_preElectrons.size(); i++){
    int iEl = m_preElectrons[i];
    const TLorentzVector &lv = m_susyObj.GetElecTLV(iEl);

    // trigger flags
    long long flags = 0;

    // 2012 triggers only

    // e7_medium1
    // NOTE: This feature is not currently available in d3pds!! Use e7T for now!
    //if( matchElectronTrigger(lv, m_event.triggerbits.trig_EF_el_EF_e7_medium1()) )
    if( matchElectronTrigger(lv, m_event.trig_EF_el.EF_e7T_medium1()) ){
      flags |= TRIG_e7_medium1;
    }
    // e12Tvh_loose1
    if( matchElectronTrigger(lv, m_event.trig_EF_el.EF_e12Tvh_loose1()) ){
      flags |= TRIG_e12Tvh_loose1;
    }
    // e12Tvh_medium1
    if( matchElectronTrigger(lv, m_event.trig_EF_el.EF_e12Tvh_medium1()) ){
      flags |= TRIG_e12Tvh_medium1;
    }
    // e24vh_medium1
    if( matchElectronTrigger(lv, m_event.trig_EF_el.EF_e24vh_medium1()) ){
      flags |= TRIG_e24vh_medium1;
    }
    // e24vhi_medium1
    if( matchElectronTrigger(lv, m_event.trig_EF_el.EF_e24vhi_medium1()) ){
      flags |= TRIG_e24vhi_medium1;
    }
    // 2e12Tvh_loose1
    if( matchElectronTrigger(lv, m_event.trig_EF_el.EF_2e12Tvh_loose1()) ){
      flags |= TRIG_2e12Tvh_loose1;
    }
    // e24vh_medium1_e7_medium1 - NOTE: you don't know which feature it matches to!!
    if( matchElectronTrigger(lv, m_event.trig_EF_el.EF_e24vh_medium1_e7_medium1()) ){
      flags |= TRIG_e24vh_medium1_e7_medium1;
    }
    // e12Tvh_medium1_mu8
    if( matchElectronTrigger(lv, m_event.trig_EF_el.EF_e12Tvh_medium1_mu8()) ){
      flags |= TRIG_e12Tvh_medium1_mu8;
    }
    // mu18_tight_e7_medium1 - NOTE: feature not available, so use e7_medium1 above!
    //if( matchElectronTrigger(lv, m_event.triggerbits.trig_EF_el_EF_mu18_tight_e7_medium1()) ){
      //flags |= TRIG_mu18_tight_e7_medium1;
    //}

    // e18vh_medium1
    if( matchElectronTrigger(lv, m_event.trig_EF_el.EF_e18vh_medium1()) ){
      flags |= TRIG_e18vh_medium1;
    }

    // e18vh_medium1_2e7T_medium1
    if( matchElectronTrigger(lv, m_event.trig_EF_el.EF_e18vh_medium1_2e7T_medium1()) ){
      flags |= TRIG_e18vh_medium1_2e7T_medium1;
    }
    // 2e7T_medium1_mu6
    if( matchElectronTrigger(lv, m_event.trig_EF_el.EF_2e7T_medium1_mu6()) ){
      flags |= TRIG_2e7T_medium1_mu6;
    }
    // e7T_medium1_2mu6
    if( matchElectronTrigger(lv, m_event.trig_EF_el.EF_e7T_medium1_2mu6()) ){
      flags |= TRIG_e7T_medium1_2mu6;
    }

    // e24vh_medium1_EFxe35_tclcw
    if( matchElectronTrigger(lv, m_event.trig_EF_el.EF_e24vh_medium1_EFxe35_tclcw()) ){
      flags |= TRIG_e24vh_medium1_EFxe35_tclcw;
    }

    // assign the trigger flags for this electron
    m_eleTrigFlags[iEl] = flags;
  }
}
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::matchElectronTrigger(const TLorentzVector &lv, vector<int>* trigBools)
{
  // matched trigger index - not used
  static int indexEF = -1;
  // Use function defined in egammaAnalysisUtils/egammaTriggerMatching.h
  return PassedTriggerEF(lv.Eta(), lv.Phi(), trigBools, indexEF, m_event.trig_EF_el.n(),
                         m_event.trig_EF_el.eta(), m_event.trig_EF_el.phi());
}

/*--------------------------------------------------------------------------------*/
// Muon trigger matching
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::matchMuonTriggers()
{
  if(m_dbg>=5) cout << "matchMuonTriggers" << endl;

  //int run = m_event.eventinfo.RunNumber();

  // loop over all pre muons
  for(uint i=0; i<m_preMuons.size(); i++){

    int iMu = m_preMuons[i];
    const TLorentzVector &lv = m_susyObj.GetMuonTLV(iMu);

    // trigger flags
    long long flags = 0;

    // 2012 triggers only

    // mu8
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu8()) ) {
      flags |= TRIG_mu8;
    }
    // mu13
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu13()) ){
      flags |= TRIG_mu13;
    }
    // mu18_tight
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu18_tight()) ) {
      flags |= TRIG_mu18_tight;
    }
    // mu24i_tight
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu24i_tight()) ) {
      flags |= TRIG_mu24i_tight;
    }
    // 2mu13
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_2mu13()) ) {
      flags |= TRIG_2mu13;
    }
    // mu18_tight_mu8_EFFS
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu18_tight_mu8_EFFS()) ) {
      flags |= TRIG_mu18_tight_mu8_EFFS;
    }
    // e12Tvh_medium1_mu8 - NOTE: muon feature not available, so use mu8
    //if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu8()) ) {
      //flags |= TRIG_e12Tvh_medium1_mu8;
    //}
    // mu18_tight_e7_medium1
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu18_tight_e7_medium1()) ) {
      flags |= TRIG_mu18_tight_e7_medium1;
    }

    // mu15
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu15()) ) {
      flags |= TRIG_mu15;
    }

    // 2mu8_EFxe40wMu_tclcw
    if(!m_isMC && m_event.eventinfo.RunNumber()>=206248 &&
       matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_2mu8_EFxe40wMu_tclcw())) {
      flags |= TRIG_2mu8_EFxe40wMu_tclcw;
    }

    // mu6
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu6()) ) {
      flags |= TRIG_mu6;
    }
    // 2mu6
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_2mu6()) ) {
      flags |= TRIG_2mu6;
    }
    // 3mu6
    //if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_3mu6()) ) {
      //flags |= TRIG_3mu6;
    //}
    // mu18_tight_2mu4_EFFS
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu18_tight_2mu4_EFFS()) ) {
      flags |= TRIG_mu18_tight_2mu4_EFFS;
    }

    // mu4T
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu4T()) ) {
      flags |= TRIG_mu4T;
    }
    // mu24
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu24()) ) {
      flags |= TRIG_mu24;
    }
    // mu4T_j65_a4tchad_xe70_tclcw_veryloose
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu4T_j65_a4tchad_xe70_tclcw_veryloose()) ) {
      flags |= TRIG_mu4T_j65_a4tchad_xe70_tclcw_veryloose;
    }
    // 2mu4T_xe60_tclcw
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_2mu4T_xe60_tclcw()) ) {
      flags |= TRIG_2mu4T_xe60_tclcw;
    }
    // 2mu8_EFxe40_tclcw
    if(m_event.trig_EF_trigmuonef.EF_2mu8_EFxe40_tclcw.IsAvailable() &&
       matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_2mu8_EFxe40_tclcw()) ) {
      flags |= TRIG_2mu8_EFxe40_tclcw;
    }
    // mu24_j65_a4tchad_EFxe40_tclcw
    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu24_j65_a4tchad_EFxe40_tclcw()) ) {
      flags |= TRIG_mu24_j65_a4tchad_EFxe40_tclcw;
    }
    // mu24_j65_a4tchad_EFxe40wMu_tclcw
    if(m_event.trig_EF_trigmuonef.EF_mu24_j65_a4tchad_EFxe40wMu_tclcw.IsAvailable() &&
       matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu24_j65_a4tchad_EFxe40wMu_tclcw()) ) {
      flags |= TRIG_mu24_j65_a4tchad_EFxe40wMu_tclcw;
    }

    // assign the trigger flags for this muon
    m_muoTrigFlags[iMu] = flags;
  }
}
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::matchMuonTrigger(const TLorentzVector &lv, vector<int>* passTrig)
{
  // loop over muon trigger features
  for(int iTrig=0; iTrig < m_event.trig_EF_trigmuonef.n(); iTrig++){

    // Check to see if this feature passed chain we want
    if(passTrig->at(iTrig)){

      // Loop over muon EF tracks
      TLorentzVector lvTrig;
      for(int iTrk=0; iTrk < m_event.trig_EF_trigmuonef.track_n()->at(iTrig); iTrk++){

        lvTrig.SetPtEtaPhiM( m_event.trig_EF_trigmuonef.track_CB_pt()->at(iTrig).at(iTrk),
                             m_event.trig_EF_trigmuonef.track_CB_eta()->at(iTrig).at(iTrk),
                             m_event.trig_EF_trigmuonef.track_CB_phi()->at(iTrig).at(iTrk),
                             0 );       // only eta and phi used to compute dR anyway
        // Require combined offline track...?
        if(!m_event.trig_EF_trigmuonef.track_CB_hasCB()->at(iTrig).at(iTrk)) continue;
        float dR = lv.DeltaR(lvTrig);
        if(dR < 0.15){
          return true;
        }

      } // loop over EF tracks
    } // trigger object passes chain?
  } // loop over trigger objects

  // matching failed
  return false;
}

/*--------------------------------------------------------------------------------*/
// Tau trigger matching
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::matchTauTriggers()
{
  if(m_dbg>=5) cout << "matchTauTriggers" << endl;

  //int run = m_event.eventinfo.RunNumber();

  // loop over all pre taus
  for(uint i=0; i<m_preTaus.size(); i++){

    int iTau = m_preTaus[i];
    const TLorentzVector &lv = m_susyObj.GetTauTLV(iTau);

    // trigger flags
    long long flags = 0;

    // tau20_medium1
    if( matchTauTrigger(lv, m_event.trig_EF_tau.EF_tau20_medium1()) ){
      flags |= TRIG_tau20_medium1;
    }
    // tau20Ti_medium1
    if( matchTauTrigger(lv, m_event.trig_EF_tau.EF_tau20Ti_medium1()) ){
      flags |= TRIG_tau20Ti_medium1;
    }
    // tau29Ti_medium1
    if( matchTauTrigger(lv, m_event.trig_EF_tau.EF_tau29Ti_medium1()) ){
      flags |= TRIG_tau29Ti_medium1;
    }
    // tau29Ti_medium1_tau20Ti_medium1
    if( matchTauTrigger(lv, m_event.trig_EF_tau.EF_tau29Ti_medium1_tau20Ti_medium1()) ){
      flags |= TRIG_tau29Ti_medium1_tau20Ti_medium1;
    }
    // tau20Ti_medium1_e18vh_medium1
    if( matchTauTrigger(lv, m_event.trig_EF_tau.EF_tau20Ti_medium1_e18vh_medium1()) ){
      flags |= TRIG_tau20Ti_medium1_e18vh_medium1;
    }
    // tau20_medium1_mu15
    if( matchTauTrigger(lv, m_event.trig_EF_tau.EF_tau20_medium1_mu15()) ){
      flags |= TRIG_tau20_medium1_mu15;
    }

    // assign the trigger flags for this tau
    m_tauTrigFlags[iTau] = flags;
  }
}
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::matchTauTrigger(const TLorentzVector &lv, vector<int>* passTrig)
{
  // loop over tau trigger features
  for(int iTrig=0; iTrig < m_event.trig_EF_tau.n(); iTrig++){
    // Check to see if this feature passed chain we want
    if(passTrig->at(iTrig)){
      // Now, try to match offline tau to this online tau
      static TLorentzVector trigLV;
      trigLV.SetPtEtaPhiM(m_event.trig_EF_tau.pt()->at(iTrig), m_event.trig_EF_tau.eta()->at(iTrig),
                          m_event.trig_EF_tau.phi()->at(iTrig), m_event.trig_EF_tau.m()->at(iTrig));
      float dR = lv.DeltaR(trigLV);
      if(dR < 0.15) return true;
    }
  }
  // matching failed
  return false;
}

/*--------------------------------------------------------------------------------*/
// Check event level cleaning cuts like GRL, LarError, etc.
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::checkEventCleaning()
{
  if(passGRL())      m_cutFlags |= ECut_GRL;
  if(passTTCVeto())  m_cutFlags |= ECut_TTC;
  if(passLarErr())   m_cutFlags |= ECut_LarErr;
  if(passTileErr())  m_cutFlags |= ECut_TileErr;
  if(passGoodVtx())  m_cutFlags |= ECut_GoodVtx;
  if(passTileTrip()) m_cutFlags |= ECut_TileTrip;
}

/*--------------------------------------------------------------------------------*/
// Check object level cleaning cuts like BadJet, BadMu, etc.
// SELECT OBJECTS BEFORE CALLING THIS!!
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::checkObjectCleaning()
{
  if(passTileHotSpot()) m_cutFlags |= ECut_HotSpot;
  if(passBadJet())      m_cutFlags |= ECut_BadJet;
  if(passBadMuon())     m_cutFlags |= ECut_BadMuon;
  if(passCosmic())      m_cutFlags |= ECut_Cosmic;
  if(passLarHoleVeto()) m_cutFlags |= ECut_SmartVeto;
}

/*--------------------------------------------------------------------------------*/
// Pass Lar hole veto
// Prior to calling this, need jet and MET selection
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::passLarHoleVeto()
{
  // LAr veto is not used anymore
  return true;

  //TVector2 metVector = m_met.Vect().XYvector();
  //vector<int> goodJets;
  //// Do I still need these jets with no eta cut?
  //// This only uses nominal jets...?
  //vector<int> jets = get_jet_baseline(&m_event.jet, &m_event.vxp, &m_event.eventinfo, !m_isMC, m_susyObj,
  //                                    20.*GeV, 9999999, SystErr::NONE, false, goodJets);
  //return !check_jet_larhole(&m_event.jet, jets, !m_isMC, m_susyObj, 180614, metVector, &m_fakeMetEst);
}
/*--------------------------------------------------------------------------------*/
// Pass tile hot spot veto
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::passTileHotSpot()
{
  D3PDReader::JetD3PDObject *jets =  d3pdJets();
  return !check_jet_tileHotSpot(jets, m_preJets, m_susyObj, !m_isMC, m_event.eventinfo.RunNumber());
}
/*--------------------------------------------------------------------------------*/
// Pass bad jet cut
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::passBadJet()
{
  D3PDReader::JetD3PDObject *jets =  d3pdJets();
  return !IsBadJetEvent(jets, m_baseJets, 20.*GeV, m_susyObj);
}
/*--------------------------------------------------------------------------------*/
// Pass good vertex
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::passGoodVtx()
{
  return PrimaryVertexCut(m_susyObj, &m_event.vxp);
}
/*--------------------------------------------------------------------------------*/
// Pass tile trip
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::passTileTrip()
{
  return !m_susyObj.IsTileTrip(m_event.eventinfo.RunNumber(), m_event.eventinfo.lbn(), m_event.eventinfo.EventNumber());
}
/*--------------------------------------------------------------------------------*/
// Pass bad muon veto
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::passBadMuon()
{
  return !IsBadMuonEvent(m_susyObj, d3pdMuons(), m_preMuons, 0.2);
}
/*--------------------------------------------------------------------------------*/
// Pass cosmic veto
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::passCosmic()
{
  return !IsCosmic(m_susyObj, d3pdMuons(), m_baseMuons, 1., 0.2);
}

/*--------------------------------------------------------------------------------*/
// Radiative b quark check for sherpa WW fix
/*--------------------------------------------------------------------------------*/
/*bool XaodAnalysis::hasRadiativeBQuark(const vector<int>* pdg, const vector<int>* status)
{
  if(!pdg || !status || pdg->size()!=status->size()) return false;
  const vector<int>& p = *pdg;
  const vector<int>& s = *status;
  const int pdgB(5), statRad(3);
  for(size_t i=0; i<p.size(); ++i) if(abs(p[i])==pdgB && s[i]==statRad) return true;
  return false;
}*/

/*--------------------------------------------------------------------------------*/
// Get event weight, combine gen, pileup, xsec, and lumi weights
// Default weight uses ICHEP dataset, A-B14 lumi
// You can supply a different luminosity, but the pileup weights will still correspond to A-B14
/*--------------------------------------------------------------------------------*/
float XaodAnalysis::getEventWeight(float lumi)
{
  if(!m_isMC) return 1;
  return m_event.eventinfo.mc_event_weight() * getXsecWeight() * getPileupWeight() * lumi / m_sumw;
}
/*--------------------------------------------------------------------------------*/
// Cross section and lumi scaling
/*--------------------------------------------------------------------------------*/
float XaodAnalysis::getXsecWeight()
{
  // Use user cross section if it has been set
  if(m_xsec > 0) return m_xsec;

  // Use SUSY cross section file
  int id = m_event.eventinfo.mc_channel_number();
  if(m_xsecMap.find(id) == m_xsecMap.end()) {
    m_xsecMap[id] = m_susyXsec->process(id);
  }
  return m_xsecMap[id].xsect() * m_xsecMap[id].kfactor() * m_xsecMap[id].efficiency();
}

/*--------------------------------------------------------------------------------*/
// Luminosity normalization
/*--------------------------------------------------------------------------------*/
float XaodAnalysis::getLumiWeight()
{ return m_lumi / m_sumw; }

/*--------------------------------------------------------------------------------*/
// Pileup reweighting
/*--------------------------------------------------------------------------------*/
float XaodAnalysis::getPileupWeight()
{
  return m_pileup->GetCombinedWeight(m_event.eventinfo.RunNumber(), m_event.eventinfo.mc_channel_number(), m_event.eventinfo.averageIntPerXing());
}
/*--------------------------------------------------------------------------------*/
float XaodAnalysis::getPileupWeightUp()
{
  return m_pileup_up->GetCombinedWeight(m_event.eventinfo.RunNumber(), m_event.eventinfo.mc_channel_number(), m_event.eventinfo.averageIntPerXing());
}
/*--------------------------------------------------------------------------------*/
float XaodAnalysis::getPileupWeightDown()
{
  return m_pileup_dn->GetCombinedWeight(m_event.eventinfo.RunNumber(), m_event.eventinfo.mc_channel_number(), m_event.eventinfo.averageIntPerXing());
}

/*--------------------------------------------------------------------------------*/
// PDF reweighting of 7TeV -> 8TeV
/*--------------------------------------------------------------------------------*/
float XaodAnalysis::getPDFWeight8TeV()
{
  #ifdef USEPDFTOOL
  float scale = m_event.mcevt.pdf_scale()->at(0);
  float x1 = m_event.mcevt.pdf_x1()->at(0);
  float x2 = m_event.mcevt.pdf_x2()->at(0);
  int id1 = m_event.mcevt.pdf_id1()->at(0);
  int id2 = m_event.mcevt.pdf_id2()->at(0);

  // MultLeip function... Not working?
  //return scaleBeamEnergy(*m_pdfTool, 21000, m_event.mcevt.pdf_scale()->at(0), m_event.mcevt.pdf_x1()->at(0),
                         //m_event.mcevt.pdf_x2()->at(0), m_event.mcevt.pdf_id1()->at(0), m_event.mcevt.pdf_id2()->at(0));
  // Simple scaling
  //return m_pdfTool->event_weight( pow(scale,2), x1, x2, id1, id2, 21000 );

  // For scaling to/from arbitrary beam energy
  m_pdfTool->setEventInfo( scale*scale, x1, x2, id1, id2 );
  //return m_pdfTool->scale((3.5+4.)/3.5);
  // possible typo correction?
  return m_pdfTool->scale(4./3.5);

  #else
  return 1;
  #endif
}

/*--------------------------------------------------------------------------------*/
// Lepton efficiency SF
/*--------------------------------------------------------------------------------*/
float XaodAnalysis::getLepSF(const vector<LeptonInfo>& leptons)
{
  // TODO: incorporate systematics
  float lepSF = 1;

  if(m_isMC){
    // Loop over leptons
    for(uint iLep=0; iLep<leptons.size(); iLep++){
      const LeptonInfo &lep = leptons[iLep];
      // Electrons
      if(lep.isElectron()){
          const D3PDReader::ElectronD3PDObjectElement* el = lep.getElectronElement();
        lepSF *= m_susyObj.GetSignalElecSF(el->cl_eta(), lep.lv()->Pt(), true, true, false);
      }
      // Muons
      else{
        lepSF *= m_susyObj.GetSignalMuonSF(lep.idx());
      }
    }
  }

  return lepSF;
}

/*--------------------------------------------------------------------------------*/
// BTag efficiency SF
// TODO: finish me!
/*--------------------------------------------------------------------------------*/
float XaodAnalysis::getBTagSF(const vector<int>& jets)
{
  return 1;
}

/*--------------------------------------------------------------------------------*/
// Calculate random MC run and lb numbers for cleaning cuts, etc.
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::calcRandomRunLB()
{
  if(m_pileup){
    m_mcRun = m_pileup->GetRandomRunNumber(m_event.eventinfo.RunNumber());
    m_mcLB = m_pileup->GetRandomLumiBlockNumber(m_mcRun);
  }
}

/*--------------------------------------------------------------------------------*/
// Get the heavy flavor overlap removal decision
/*--------------------------------------------------------------------------------*/
int XaodAnalysis::getHFORDecision()
{
  return m_hforTool.getDecision(m_event.eventinfo.mc_channel_number(),
                                m_event.mc.n(),
                                m_event.mc.pt(),
                                m_event.mc.eta(),
                                m_event.mc.phi(),
                                m_event.mc.m(),
                                m_event.mc.pdgId(),
                                m_event.mc.status(),
                                m_event.mc.vx_barcode(),
                                m_event.mc.parent_index(),
                                m_event.mc.child_index(),
                                HforToolD3PD::ALL); //HforToolD3PD::DEFAULT
}

/*--------------------------------------------------------------------------------*/
// Set MET flavor via a string. Only a couple of options available so far
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::setMetFlavor(string metFlav)
{
  if(metFlav=="STVF") m_metFlavor = SUSYMet::STVF;
  else if(metFlav=="STVF_JVF") m_metFlavor = SUSYMet::STVF_JVF;
  else if(metFlav=="Default") m_metFlavor = SUSYMet::Default;
  else{
    cout << "XaodAnalysis::setMetFlavor : ERROR : MET flavor " << metFlav
         << " is not supported!" << endl;
    abort();
  }
}

/*--------------------------------------------------------------------------------*/
// Print event info
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::dumpEvent()
{
  cout << "Run " << setw(6) << m_event.eventinfo.RunNumber()
       << " Event " << setw(7) << m_event.eventinfo.EventNumber()
       << " Stream " << streamName(m_stream);
  if(m_isMC){
    cout << " MCID " << setw(6) << m_event.eventinfo.mc_channel_number()
         << " weight " << m_event.eventinfo.mc_event_weight();
  }
  cout << endl;
}

/*--------------------------------------------------------------------------------*/
// Print baseline objects
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::dumpBaselineObjects()
{
  uint nEle = m_baseElectrons.size();
  uint nMu  = m_baseMuons.size();
  //uint nTau = m_baseTaus.size();
  uint nJet = m_baseJets.size();

  cout.precision(2);
  if(nEle){
    cout << "Baseline electrons" << endl;
    for(uint i=0; i < nEle; i++){
      int iEl = m_baseElectrons[i];
      const TLorentzVector &lv = m_susyObj.GetElecTLV(iEl);
      const D3PDReader::ElectronD3PDObjectElement &ele = (*d3pdElectrons())[iEl];
      cout << "  El : " << fixed
           << " q " << setw(2) << (int) ele.charge()
           << " pt " << setw(6) << lv.Pt()/GeV
           << " eta " << setw(5) << lv.Eta()
           << " phi " << setw(5) << lv.Phi();
      if(m_isMC) cout << " type " << setw(2) << ele.type() << " origin " << setw(2) << ele.origin();
      cout << endl;
    }
  }
  if(nMu){
    cout << "Baseline muons" << endl;
    for(uint i=0; i < nMu; i++){
      int iMu = m_baseMuons[i];
      const TLorentzVector &lv = m_susyObj.GetMuonTLV(iMu);
      const D3PDReader::MuonD3PDObjectElement &muo = (*d3pdMuons())[iMu];
      cout << "  Mu : " << fixed
           << " q " << setw(2) << (int) muo.charge()
           << " pt " << setw(6) << lv.Pt()/GeV
           << " eta " << setw(5) << lv.Eta()
           << " phi " << setw(5) << lv.Phi();
      if(m_isMC) cout << " type " << setw(2) << muo.type() << " origin " << setw(2) << muo.origin();
      cout << endl;
    }
  }
  if(nJet){
    cout << "Baseline jets" << endl;
    for(uint i=0; i < nJet; i++){
      int iJet = m_baseJets[i];
      const TLorentzVector &lv = m_susyObj.GetJetTLV(iJet);
      const D3PDReader::JetD3PDObjectElement &jet = (*d3pdJets())[iJet];
      cout << "  Jet : " << fixed
           << " pt " << setw(6) << lv.Pt()/GeV
           << " eta " << setw(5) << lv.Eta()
           << " phi " << setw(5) << lv.Phi()
           << " mv1 " << jet.flavor_weight_MV1();
      cout << endl;
    }
  }
  cout.precision(6);
  cout.unsetf(ios_base::fixed);
}

/*--------------------------------------------------------------------------------*/
// Print signal objects
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::dumpSignalObjects()
{
  uint nEle = m_sigElectrons.size();
  uint nMu  = m_sigMuons.size();
  //uint nTau = m_sigTaus.size();
  uint nJet = m_sigJets.size();

  cout.precision(2);
  if(nEle){
    cout << "Signal electrons" << endl;
    for(uint i=0; i < nEle; i++){
      int iEl = m_sigElectrons[i];
      const TLorentzVector &lv = m_susyObj.GetElecTLV(iEl);
      const D3PDReader::ElectronD3PDObjectElement &ele = (*d3pdElectrons())[iEl];
      cout << "  El : " << fixed
           << " q " << setw(2) << (int) ele.charge()
           << " pt " << setw(6) << lv.Pt()/GeV
           << " eta " << setw(5) << lv.Eta()
           << " phi " << setw(5) << lv.Phi();
      if(m_isMC) cout << " type " << setw(2) << ele.type() << " origin " << setw(2) << ele.origin();
      cout << endl;
    }
  }
  if(nMu){
    cout << "Signal muons" << endl;
    for(uint i=0; i < nMu; i++){
      int iMu = m_sigMuons[i];
      const TLorentzVector &lv = m_susyObj.GetMuonTLV(iMu);
      const D3PDReader::MuonD3PDObjectElement &muo = (*d3pdMuons())[iMu];
      cout << "  Mu : " << fixed
           << " q " << setw(2) << (int) muo.charge()
           << " pt " << setw(6) << lv.Pt()/GeV
           << " eta " << setw(5) << lv.Eta()
           << " phi " << setw(5) << lv.Phi();
      if(m_isMC) cout << " type " << setw(2) << muo.type() << " origin " << setw(2) << muo.origin();
      cout << endl;
    }
  }
  if(nJet){
    cout << "Signal jets" << endl;
    for(uint i=0; i < nJet; i++){
      int iJet = m_sigJets[i];
      const TLorentzVector &lv = m_susyObj.GetJetTLV(iJet);
      const D3PDReader::JetD3PDObjectElement &jet = (*d3pdJets())[iJet];
      cout << "  Jet : " << fixed
           << " pt " << setw(6) << lv.Pt()/GeV
           << " eta " << setw(5) << lv.Eta()
           << " phi " << setw(5) << lv.Phi()
           << " mv1 " << jet.flavor_weight_MV1();
      cout << endl;
    }
  }
  cout.precision(6);
  cout.unsetf(ios_base::fixed);
}
//----------------------------------------------------------
bool XaodAnalysis::runningOptionsAreValid()
{
    bool valid=true;
    bool isSimulation = m_event.eventinfo.isSimulation();
    bool isData = !isSimulation;
    bool isStreamEgamma = m_event.eventinfo.streamDecision_Egamma();
    bool isStreamJetEt  = m_event.eventinfo.streamDecision_JetTauEtmiss();
    bool isStreamMuons  = m_event.eventinfo.streamDecision_Muons();
    if(m_isMC != isSimulation) {
        valid=false;
        if(m_dbg)
            cout<<"XaodAnalysis::runningOptionsAreValid invalid isMc:"
                <<" (m_isMC:"<<m_isMC<<" != isSimulation:"<<isSimulation<<")"
                <<endl;
    }
    if(isData) {
        bool consistentStream = (isStreamMuons  ? m_stream==Stream_Muons :
                                 isStreamEgamma ? m_stream==Stream_Egamma :
                                 isStreamJetEt  ? m_stream==Stream_JetTauEtmiss :
                                 false);
        if(!consistentStream) {
            valid=false;
            if(m_dbg)
                cout<<"XaodAnalysis::runningOptionsAreValid: inconsistent stream"
                    <<" m_stream: "<<(m_stream==Stream_Muons        ? "Stream_Muons":
                                      m_stream==Stream_Egamma       ? "Stream_Egamma":
                                      m_stream==Stream_JetTauEtmiss ? "Stream_JetTauEtmiss":
                                      "unknown")
                    <<" eventinfo: "<<(isStreamMuons ? "Muons":
                                       isStreamEgamma ? "Egamma":
                                       isStreamJetEt ? "JetTauEtmiss":
                                       "unknown")
                    <<endl;

        }
    }
    if(m_dbg)
        cout<<"XaodAnalysis::runningOptionsAreValid(): "<<(valid?"true":"false")<<endl;
    return valid;
}
//----------------------------------------------------------
std::string XaodAnalysis::defauldGrlFile()
{
    return std::string( "$ROOTCOREBIN/../SUSYTools/data/GRL/Summer2013/"
                        "data12_8TeV.periodAllYear_DetStatus-v61-pro14-02"
                        "_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml");
}
//----------------------------------------------------------
bool XaodAnalysis::initGrlTool()
{
    bool success = false;
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    std::vector<std::string> grl_files;
    grl_files.push_back(XaodAnalysis::defauldGrlFile());
    m_grl->setProperty("GoodRunsListVec", grl_files);
    m_grl->setProperty("PassThrough", false);
    success = m_grl->initialize(); // DG any check we should do here? (file_exists?)
    return success;
}
//----------------------------------------------------------

#undef GeV
