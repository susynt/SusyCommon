#include "TSystem.h"
#include "SusyCommon/SusyD3PDAna.h"
#include "MultiLep/ElectronTools.h"
#include "MultiLep/MuonTools.h"
#include "MultiLep/JetTools.h"
#include "MultiLep/TauTools.h"
#include "MultiLep/CutflowTools.h"

#include "MultiLep/PhotonTools.h"

#ifdef USEPDFTOOL
#include "MultiLep/PDFErrorTools.h"
#endif

using namespace std;

#define GeV 1000.

/*--------------------------------------------------------------------------------*/
// SusyD3PDAna Constructor
/*--------------------------------------------------------------------------------*/
SusyD3PDAna::SusyD3PDAna() : 
        m_sample(""),
        m_isAF2(false),
        m_d3pdTag(D3PD_p1328),
        m_selectPhotons(true),
        m_selectTaus(false),
        m_selectTruth(false),
        m_metFlavor(SUSYMet::STVF),
        m_doMetFix(false),
	//m_useMetMuons(false),
        m_lumi(LUMI_A_E),
        m_sumw(1),
	m_xsec(-1),
	m_errXsec(-1),
	m_sys(false),
        m_eleMediumSFTool(0),
        m_pileup(0),
        m_pileupAB3(0),
        m_pileupAB(0),
        m_pileupAE(0),
        m_susyXsec(0),
        m_hforTool()
{
  // Now works!
  m_hforTool.setVerbosity(HforToolD3PD::ERROR);

  // Create the addition electron efficiency SF tool for medium SFs
  m_eleMediumSFTool = new Root::TElectronEfficiencyCorrectionTool;

  #ifdef USEPDFTOOL
  m_pdfTool = new PDFTool(3500000, 1, -1, 21000);
  //m_pdfTool = new PDFTool(3500000, 4./3.5);
  #endif
}
/*--------------------------------------------------------------------------------*/
// Destructor
/*--------------------------------------------------------------------------------*/
SusyD3PDAna::~SusyD3PDAna()
{
  if(m_eleMediumSFTool) delete m_eleMediumSFTool;
  #ifdef USEPDFTOOL
  if(m_pdfTool) delete m_pdfTool;
  #endif
}
/*--------------------------------------------------------------------------------*/
// The Begin() function is called at the start of the query.
// When running with PROOF Begin() is only called on the client.
// The tree argument is deprecated (on PROOF 0 is passed).
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::Begin(TTree* /*tree*/)
{
  if(m_dbg) cout << "SusyD3PDAna::Begin" << endl;

  // Use sample name to set MC flag
  if(m_sample.Contains("data", TString::kIgnoreCase)) {
    m_isMC = false;
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
  // TODO: FINISH ME
  m_susyObj.initialize(!m_isMC, m_isAF2, false,
		       gSystem->ExpandPathName("$ROOTCOREDIR/data/MuonMomentumCorrections/"),
		       gSystem->ExpandPathName("$ROOTCOREDIR/data/MuonEfficiencyCorrections/"));
                       //"STACO_CB_plus_ST",
                       //"efficiencySF.offline.RecoTrk.2012.8TeV.rel17p2.v02.root",
                       //"efficiencySF.offline.Tight.2012.8TeV.rel17p2.v02.root",
                       //"efficiencySF.e24vhi_medium1_e60_medium1.Tight.2012.8TeV.rel17p2.v02.root",
                       //gSystem->ExpandPathName("$ROOTCOREDIR/data/MultiLep"));
                       
  // Turn on jet calibration
  m_susyObj.SetJetCalib(true);

  // Initialize electron medium SF
  string eleMedFile = getenv("ROOTCOREDIR");
  eleMedFile += "/data/ElectronEfficiencyCorrection/efficiencySF.offline.Medium.2012.8TeV.rel17p2.v01.root";
  m_eleMediumSFTool->addFileName(eleMedFile.c_str());
  if(!m_eleMediumSFTool->initialize()){
    cout << "SusyD3PDAna::Begin : ERROR initializing TElectronEfficiencyCorrectionTool with file "
         << eleMedFile << endl;
    abort();
  }

  // Set the MissingEt flag for STVF
  // This is now done automatically when you call SUSYObjDef::GetMET
  //if(m_metFlavor.Contains("STVF")) m_susyObj.GetMETUtility()->configMissingET(true, true);

  m_fakeMetEst.initialize("$ROOTCOREDIR/data/MultiLep/fest_periodF_v1.root");

  // SUSY cross sections
  if(m_isMC){
    // Back to using the SUSYTools file
    //string xsecFileName  = gSystem->ExpandPathName("$ROOTCOREDIR/data/MultiLep/susy_crosssections_8TeV_mod.txt");
    string xsecFileName  = gSystem->ExpandPathName("$ROOTCOREDIR/data/SUSYTools/susy_crosssections_8TeV.txt");
    m_susyXsec = new SUSY::CrossSectionDB(xsecFileName);
  }

  // GRL
  if(!m_isMC){
    if(m_grlFileName.Length() == 0){
      m_grlFileName = gSystem->ExpandPathName("$ROOTCOREDIR/data/MultiLep/data12_8TeV.periodAllYear_DetStatus-v58-pro14-01_DQDefects-00-00-33_PHYS_StandardGRL_All_Good.xml");
    }
    
    Root::TGoodRunsListReader* grlReader = new Root::TGoodRunsListReader();
    grlReader->AddXMLFile(m_grlFileName);
    grlReader->Interpret();
    m_grl = grlReader->GetMergedGoodRunsList();
    delete grlReader;
  }

  // Pileup reweighting
  if(m_isMC){
    m_pileup = new Root::TPileupReweighting("PileupReweighting");
    m_pileup->SetDataScaleFactors(1/1.11);
    m_pileup->AddConfigFile("$ROOTCOREDIR/data/PileupReweighting/mc12a_defaults.prw.root");
    m_pileup->AddLumiCalcFile("$ROOTCOREDIR/data/MultiLep/ilumicalc_histograms_EF_2e12Tvh_loose1_200842-215643_v2.root");
    m_pileup->SetUnrepresentedDataAction(2);
    int pileupError = m_pileup->Initialize();
    if(pileupError){
      cout << "Problem in pileup initialization.  pileupError = " << pileupError << endl;
      abort();
    }

    // pileup reweighting for 2012 A-B3 only
    m_pileupAB3 = new Root::TPileupReweighting("PileupReweightingAB3");
    m_pileupAB3->SetDataScaleFactors(1/1.11);
    m_pileupAB3->AddConfigFile("$ROOTCOREDIR/data/PileupReweighting/mc12a_defaults.prw.root");
    m_pileupAB3->AddLumiCalcFile("$ROOTCOREDIR/data/MultiLep/ilumicalc_histograms_EF_2e12Tvh_loose1_200841-203195.root");
    m_pileupAB3->SetUnrepresentedDataAction(2);
    pileupError = m_pileupAB3->Initialize();
    if(pileupError){
      cout << "Problem in pileup initialization.  pileupError = " << pileupError << endl;
      abort();
    }
    // pileup reweighting for 2012 A-B only
    m_pileupAB = new Root::TPileupReweighting("PileupReweightingAB");
    m_pileupAB->SetDataScaleFactors(1/1.11);
    m_pileupAB->AddConfigFile("$ROOTCOREDIR/data/PileupReweighting/mc12a_defaults.prw.root");
    m_pileupAB->AddLumiCalcFile("$ROOTCOREDIR/data/MultiLep/ilumicalc_histograms_EF_e24vhi_medium1_200841-205113.root");
    m_pileupAB->SetUnrepresentedDataAction(2);
    pileupError = m_pileupAB->Initialize();
    if(pileupError){
      cout << "Problem in pileup initialization.  pileupError = " << pileupError << endl;
      abort();
    }
    // pileup reweighting for 2012 A-E only
    m_pileupAE = new Root::TPileupReweighting("PileupReweightingAE");
    m_pileupAE->SetDataScaleFactors(1/1.11);
    m_pileupAE->AddConfigFile("$ROOTCOREDIR/data/PileupReweighting/mc12a_defaults.prw.root");
    m_pileupAE->AddLumiCalcFile("$ROOTCOREDIR/data/MultiLep/ilumicalc_histograms_EF_2e12Tvh_loose1_200841-210308.root");
    m_pileupAE->SetUnrepresentedDataAction(2);
    pileupError = m_pileupAE->Initialize();
    if(pileupError){
      cout << "Problem in pileup initialization.  pileupError = " << pileupError << endl;
      abort();
    }
  }
}

/*--------------------------------------------------------------------------------*/
// Main process loop function - This is just an example for testing
/*--------------------------------------------------------------------------------*/
Bool_t SusyD3PDAna::Process(Long64_t entry)
{
  // Communicate the entry number to the interface objects
  GetEntry(entry);

  static Long64_t chainEntry = -1;
  chainEntry++;
  if(m_dbg || chainEntry%10000==0)
  {
    cout << "**** Processing entry " << setw(6) << chainEntry
         << " run " << setw(6) << d3pd.evt.RunNumber()
         << " event " << setw(7) << d3pd.evt.EventNumber() << " ****" << endl;
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
void SusyD3PDAna::Terminate()
{
  if(m_dbg) cout << "SusyD3PDAna::Terminate" << endl;
  m_susyObj.finalize();

  if(m_isMC){
    delete m_susyXsec;
    delete m_pileup;
    delete m_pileupAB3;
    delete m_pileupAB;
    delete m_pileupAE;
  }
}

/*--------------------------------------------------------------------------------*/
// Baseline object selection
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::selectBaselineObjects(SusyNtSys sys)
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

  // Preselection
  m_preElectrons = get_electrons_baseline(&d3pd.ele, !m_isMC, d3pd.evt.RunNumber(), m_susyObj, 
                                          10.*GeV, 2.47, susySys);
  m_preMuons     = get_muons_baseline(&d3pd.muo, !m_isMC, m_susyObj, 
                                      10.*GeV, 2.5, susySys);
  m_preJets      = get_jet_baseline(&d3pd.jet, &d3pd.vtx, &d3pd.evt, !m_isMC, m_susyObj, 
                                    20.*GeV, 4.9, susySys, false, goodJets);

  // Selection for met muons
  // Diff with preMuons is |eta| selection
  m_metMuons     = get_muons_baseline(&d3pd.muo, !m_isMC, m_susyObj, 
                                      10.*GeV, 2.5, susySys);

  // Preselect taus
  if(m_selectTaus){
    m_preTaus    = get_taus_baseline(&d3pd.tau, m_susyObj, 20.*GeV, 2.5, 
                                     SUSYTau::TauLoose, SUSYTau::TauLoose, SUSYTau::TauLoose, 
                                     susySys);
  }
  else{
    m_preTaus.clear();
  }
  
  performOverlapRemoval();

  // combine leptons
  m_preLeptons    = buildLeptonInfos(&d3pd.ele, m_preElectrons, &d3pd.muo, m_preMuons, m_susyObj);
  m_baseLeptons   = buildLeptonInfos(&d3pd.ele, m_baseElectrons, &d3pd.muo, m_baseMuons, m_susyObj);
}

/*--------------------------------------------------------------------------------*/
// perform overlap
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::performOverlapRemoval()
{
  // e-e overlap removal
  m_baseElectrons = overlap_removal(m_susyObj, &d3pd.ele, m_preElectrons, &d3pd.ele, m_preElectrons, 
                                    0.1, true, true);

  // jet-e overlap removal
  m_baseJets      = overlap_removal(m_susyObj, &d3pd.jet, m_preJets, &d3pd.ele, m_baseElectrons, 
                                    0.2, false, false);

  if(m_selectTaus) {
    // tau-e overlap removal
    m_baseTaus    = overlap_removal(m_susyObj, &d3pd.tau, m_preTaus, &d3pd.ele, m_baseElectrons, 0.2, false, false);
    // tau-mu overlap removal
    m_baseTaus    = overlap_removal(m_susyObj, &d3pd.tau, m_baseTaus, &d3pd.muo, m_preMuons, 0.2, false, false);
    // jet-tau overlap removal
    // TODO: Figure out how to deal with this!!
    // Oops, this actually screws up the bad jet veto.  I can't apply this here...
    //m_baseJets    = overlap_removal(m_susyObj, &d3pd.jet, m_baseJets, &d3pd.tau, m_baseTaus, 0.2, false, false);
  }

  // e-jet overlap removal
  m_baseElectrons = overlap_removal(m_susyObj, &d3pd.ele, m_baseElectrons, &d3pd.jet, m_baseJets, 
                                    0.4, false, false);

  // m-jet overlap removal
  m_baseMuons     = overlap_removal(m_susyObj, &d3pd.muo, m_preMuons, &d3pd.jet, m_baseJets, 0.4, false, false);

  // e-m overlap removal
  vector<int> copyElectrons = m_baseElectrons;
  m_baseElectrons = overlap_removal(m_susyObj, &d3pd.ele, m_baseElectrons, &d3pd.muo, m_baseMuons, 
                                    0.1, false, false);
  m_baseMuons     = overlap_removal(m_susyObj, &d3pd.muo, m_baseMuons, &d3pd.ele, copyElectrons, 0.1, false, false);

  // m-m overlap removal - need to validate this!!
  m_baseMuons     = overlap_removal(m_susyObj, &d3pd.muo, m_baseMuons, &d3pd.muo, m_baseMuons, 0.05, true, false);

  // remove SFOS lepton pairs with Mll < 12 GeV
  m_baseElectrons = RemoveSFOSPair(m_susyObj, &d3pd.ele, m_baseElectrons, 12.*GeV);
  m_baseMuons     = RemoveSFOSPair(m_susyObj, &d3pd.muo, m_baseMuons,     12.*GeV);
  m_baseTaus      = RemoveSFOSPair(m_susyObj, &d3pd.tau, m_baseTaus,      12.*GeV);
}

/*--------------------------------------------------------------------------------*/
// Signal object selection - do baseline selection first!
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::selectSignalObjects()
{
  if(m_dbg>=5) cout << "selectSignalObjects" << endl;
  uint nVtx = getNumGoodVtx();
  m_sigElectrons = get_electrons_signal(&d3pd.ele, m_baseElectrons, nVtx, !m_isMC, m_susyObj, 
                                        10.*GeV, 0.16, 0.18, 5., 0.4);
  m_sigMuons     = get_muons_signal(&d3pd.muo, m_baseMuons, nVtx, !m_isMC, m_susyObj, 10.*GeV, .12, 3., 1.);
  m_sigJets      = get_jet_signal(&d3pd.jet, m_susyObj, m_baseJets, 20.*GeV, 2.5, 0.75);
  m_sigTaus      = get_taus_signal(&d3pd.tau, m_baseTaus, m_susyObj);

  // combine light leptons
  m_sigLeptons   = buildLeptonInfos(&d3pd.ele, m_sigElectrons, &d3pd.muo, m_sigMuons, m_susyObj);

  // photon selection done in separate method, why?
  if(m_selectPhotons) selectSignalPhotons();
}

/*--------------------------------------------------------------------------------*/
// Build MissingEt
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::buildMet(SusyNtSys sys)
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
  vector<int> metElectrons = get_electrons_met(&d3pd.ele, m_susyObj);

  // Calculate the MET
  // We shouldn't need the useMetMuons flag anymore, since baseline muons will extend to eta=2.5
  TVector2 metVector = GetMetVector(m_susyObj, &d3pd.jet, &d3pd.muo, &d3pd.ele, &d3pd.met, 
                                    &d3pd.evt, m_preMuons, m_baseElectrons, metElectrons, 
                                    susySys, m_metFlavor, m_doMetFix);
  //TVector2 metVector = m_useMetMuons? 
  //                     GetMetVector(m_susyObj, &d3pd.jet, &d3pd.muo, &d3pd.ele, &d3pd.met, 
  //                                  &d3pd.evt, m_preMuons, m_baseElectrons, metElectrons, 
  //                                  susySys, m_metFlavor, m_doMetFix) :
  //                     GetMetVector(m_susyObj, &d3pd.jet, &d3pd.muo, &d3pd.ele, &d3pd.met, 
  //                                  &d3pd.evt, m_metMuons, m_baseElectrons, metElectrons, 
  //                                  susySys, m_metFlavor, m_doMetFix);
    
  m_met.SetPxPyPzE(metVector.X(), metVector.Y(), 0, metVector.Mod());
}

/*--------------------------------------------------------------------------------*/
// Signal photons
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::selectSignalPhotons()
{
  if(m_dbg>=5) cout << "selectSignalPhotons" << endl;

  int phoQual = 2;      // Quality::Tight
  uint isoType = 1;     // Corresponds to PTED corrected isolation 
  float etcone40CorrCut = 3*GeV; 

  vector<int> base_photons = get_photons_baseline(&d3pd.pho, !m_isMC, d3pd.evt.RunNumber(), m_susyObj, 
						  20.*GeV, 2.47, SystErr::NONE, phoQual);

  // Uncomment when PhotonTools Updated in MultiLep
  //vector<int> base_photons = get_photons_baseline(&d3pd.pho, !m_isMC, d3pd.evt.RunNumber(), m_susyObj, 
  //20.*GeV, 2.47, SystErr::NONE, phoQual, true);

  // Latest and Greatest
  int nPV = getNumGoodVtx();
  m_sigPhotons = get_photons_signal(&d3pd.pho, base_photons, m_susyObj, nPV, !m_isMC, 20.*GeV, etcone40CorrCut, isoType);
}
/*--------------------------------------------------------------------------------*/
// Truth object selection
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::selectTruthObjects()
{
  if(m_dbg>=5) cout << "selectTruthObjects" << endl;

  // ==>> First the truth particles
  // Done under SusyNtMaker::fillTruthParticleVars

  // ==>> Second the truth jets
  for(int index=0; index < d3pd.truthJet.n(); index++) {
    const TruthJetElement* trueJet = & d3pd.truthJet[index];
    if( trueJet->pt()/GeV > 15. && fabs(trueJet->eta()) < 4.5) m_truJets.push_back(index);
  }

  // ==>> Third and last the truth met
  m_truMet.SetPxPyPzE(d3pd.met.Truth_NonInt_etx(), d3pd.met.Truth_NonInt_ety(), 0, d3pd.met.Truth_NonInt_sumet());
}

/*--------------------------------------------------------------------------------*/
// Clear selected objects
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::clearObjects()
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
uint SusyD3PDAna::getNumGoodVtx()
{
  uint nVtx = 0;
  for(int i=0; i < d3pd.vtx.n(); i++){
    if(d3pd.vtx.nTracks()->at(i) >= 5) nVtx++;
  }
  return nVtx;
}

/*--------------------------------------------------------------------------------*/
// Match reco jet to a truth jet
/*--------------------------------------------------------------------------------*/
bool SusyD3PDAna::matchTruthJet(int iJet)
{
  // Loop over truth jets looking for a match
  const TLorentzVector* jetLV = & m_susyObj.GetJetTLV(iJet);
  for(int i=0; i<d3pd.truthJet.n(); i++){
    const TruthJetElement* trueJet = & d3pd.truthJet[i];
    TLorentzVector trueJetLV;
    trueJetLV.SetPtEtaPhiE(trueJet->pt(), trueJet->eta(), trueJet->phi(), trueJet->E());
    if(jetLV->DeltaR(trueJetLV) < 0.3) return true;
  }
  return false;
}

/*--------------------------------------------------------------------------------*/
// Event trigger flags
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::fillEventTriggers()
{
  if(m_dbg>=5) cout << "fillEventTriggers" << endl;
  
  m_evtTrigFlags = 0;
  // e7_medium1 not available at the moment, so use e7T for now
  //if(d3pd.trig.EF_e7_medium1())                 m_evtTrigFlags |= TRIG_e7_medium1;
  if(d3pd.trig.EF_e7T_medium1())                m_evtTrigFlags |= TRIG_e7_medium1;
  if(d3pd.trig.EF_e12Tvh_loose1())              m_evtTrigFlags |= TRIG_e12Tvh_loose1;
  if(d3pd.trig.EF_e12Tvh_medium1())             m_evtTrigFlags |= TRIG_e12Tvh_medium1;
  if(d3pd.trig.EF_e24vh_medium1())              m_evtTrigFlags |= TRIG_e24vh_medium1;
  if(d3pd.trig.EF_e24vhi_medium1())             m_evtTrigFlags |= TRIG_e24vhi_medium1;
  if(d3pd.trig.EF_2e12Tvh_loose1())             m_evtTrigFlags |= TRIG_2e12Tvh_loose1;
  if(d3pd.trig.EF_e24vh_medium1_e7_medium1())   m_evtTrigFlags |= TRIG_e24vh_medium1_e7_medium1;
  if(d3pd.trig.EF_mu8())                        m_evtTrigFlags |= TRIG_mu8;
  if(d3pd.trig.EF_mu13())                       m_evtTrigFlags |= TRIG_mu13;
  if(d3pd.trig.EF_mu18_tight())                 m_evtTrigFlags |= TRIG_mu18_tight;
  if(d3pd.trig.EF_mu24i_tight())                m_evtTrigFlags |= TRIG_mu24i_tight;
  if(d3pd.trig.EF_2mu13())                      m_evtTrigFlags |= TRIG_2mu13;
  if(d3pd.trig.EF_mu18_tight_mu8_EFFS())        m_evtTrigFlags |= TRIG_mu18_tight_mu8_EFFS;
  if(d3pd.trig.EF_e12Tvh_medium1_mu8())         m_evtTrigFlags |= TRIG_e12Tvh_medium1_mu8;
  if(d3pd.trig.EF_mu18_tight_e7_medium1())      m_evtTrigFlags |= TRIG_mu18_tight_e7_medium1;

  if(d3pd.trig.EF_tau20_medium1())                   m_evtTrigFlags |= TRIG_tau20_medium1;
  if(d3pd.trig.EF_tau20Ti_medium1())                 m_evtTrigFlags |= TRIG_tau20Ti_medium1;
  if(d3pd.trig.EF_tau29Ti_medium1())                 m_evtTrigFlags |= TRIG_tau29Ti_medium1;
  if(d3pd.trig.EF_tau29Ti_medium1_tau20Ti_medium1()) m_evtTrigFlags |= TRIG_tau29Ti_medium1_tau20Ti_medium1;
  if(d3pd.trig.EF_tau20Ti_medium1_e18vh_medium1())   m_evtTrigFlags |= TRIG_tau20Ti_medium1_e18vh_medium1;
  if(d3pd.trig.EF_tau20_medium1_mu15())              m_evtTrigFlags |= TRIG_tau20_medium1_mu15;

  if(d3pd.trig.EF_e18vh_medium1())              m_evtTrigFlags |= TRIG_e18vh_medium1;
  if(d3pd.trig.EF_mu15())                       m_evtTrigFlags |= TRIG_mu15;

  // EF_2mu8_EFxe40wMu_tclcw trigger only available for data, in periods > B
  if(!m_isMC && d3pd.evt.RunNumber()>=206248 && d3pd.trig.EF_2mu8_EFxe40wMu_tclcw()) 
    m_evtTrigFlags |= TRIG_2mu8_EFxe40wMu_tclcw;
}

/*--------------------------------------------------------------------------------*/
// Electron trigger matching
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::matchElectronTriggers()
{
  if(m_dbg>=5) cout << "matchElectronTriggers" << endl;
  //int run = d3pd.evt.RunNumber();

  // loop over all pre electrons
  for(uint i=0; i<m_preElectrons.size(); i++){
    int iEl = m_preElectrons[i];
    const TLorentzVector* lv = & m_susyObj.GetElecTLV(iEl);

    // trigger flags
    uint flags = 0;

    // 2012 triggers only

    // e7_medium1
    // NOTE: This feature is not currently available in d3pds!! Use e7T for now!
    //if( matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_e7_medium1()) )
    if( matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_e7T_medium1()) ){
      flags |= TRIG_e7_medium1;
    }
    // e12Tvh_loose1
    if( matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_e12Tvh_loose1()) ){
      flags |= TRIG_e12Tvh_loose1;
    }
    // e12Tvh_medium1
    if( matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_e12Tvh_medium1()) ){
      flags |= TRIG_e12Tvh_medium1;
    }
    // e24vh_medium1
    if( matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_e24vh_medium1()) ){
      flags |= TRIG_e24vh_medium1;
    }
    // e24vhi_medium1
    if( matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_e24vhi_medium1()) ){
      flags |= TRIG_e24vhi_medium1;
    }
    // 2e12Tvh_loose1
    if( matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_2e12Tvh_loose1()) ){
      flags |= TRIG_2e12Tvh_loose1;
    }
    // e24vh_medium1_e7_medium1 - NOTE: you don't know which feature it matches to!!
    if( matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_e24vh_medium1_e7_medium1()) ){
      flags |= TRIG_e24vh_medium1_e7_medium1;
    }
    // e12Tvh_medium1_mu8
    if( matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_e12Tvh_medium1_mu8()) ){
      flags |= TRIG_e12Tvh_medium1_mu8;
    }
    // mu18_tight_e7_medium1 - NOTE: feature not available, so use e7_medium1 above!
    //if( matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_mu18_tight_e7_medium1()) ){
      //flags |= TRIG_mu18_tight_e7_medium1;
    //}

    // e18vh_medium1
    if( matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_e18vh_medium1()) ){
      flags |= TRIG_e18vh_medium1;
    }

    // assign the trigger flags for this electron
    m_eleTrigFlags[iEl] = flags;
  }
}
/*--------------------------------------------------------------------------------*/
bool SusyD3PDAna::matchElectronTrigger(const TLorentzVector* lv, vector<int>* trigBools)
{
  // matched trigger index - not used
  static int indexEF = -1;
  // Use function defined in egammaAnalysisUtils/egammaTriggerMatching.h
  return PassedTriggerEF(lv->Eta(), lv->Phi(), trigBools, indexEF, d3pd.trig.trig_EF_el_n(), 
                         d3pd.trig.trig_EF_el_eta(), d3pd.trig.trig_EF_el_phi());
}

/*--------------------------------------------------------------------------------*/
// Muon trigger matching
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::matchMuonTriggers()
{
  if(m_dbg>=5) cout << "matchMuonTriggers" << endl;

  //int run = d3pd.evt.RunNumber();

  // loop over all pre muons
  for(uint i=0; i<m_preMuons.size(); i++){

    int iMu = m_preMuons[i];
    const TLorentzVector* lv = & m_susyObj.GetMuonTLV(iMu);
    
    // trigger flags
    uint flags = 0;

    // 2012 triggers only

    // mu8
    if( matchMuonTrigger(lv, d3pd.trig.trig_EF_trigmuonef_EF_mu8()) ) {
      flags |= TRIG_mu8;
    }
    // mu13
    if( matchMuonTrigger(lv, d3pd.trig.trig_EF_trigmuonef_EF_mu13()) ){
      flags |= TRIG_mu13;
    }
    // mu18_tight
    if( matchMuonTrigger(lv, d3pd.trig.trig_EF_trigmuonef_EF_mu18_tight()) ) {
      flags |= TRIG_mu18_tight;
    }
    // mu24i_tight
    if( matchMuonTrigger(lv, d3pd.trig.trig_EF_trigmuonef_EF_mu24i_tight()) ) {
      flags |= TRIG_mu24i_tight;
    }
    // 2mu13
    if( matchMuonTrigger(lv, d3pd.trig.trig_EF_trigmuonef_EF_2mu13()) ) {
      flags |= TRIG_2mu13;
    }
    // mu18_tight_mu8_EFFS
    if( matchMuonTrigger(lv, d3pd.trig.trig_EF_trigmuonef_EF_mu18_tight_mu8_EFFS()) ) {
      flags |= TRIG_mu18_tight_mu8_EFFS;
    }
    // e12Tvh_medium1_mu8 - NOTE: muon feature not available, so use mu8
    //if( matchMuonTrigger(lv, d3pd.trig.trig_EF_trigmuonef_EF_mu8()) ) {
      //flags |= TRIG_e12Tvh_medium1_mu8;
    //}
    // mu18_tight_e7_medium1
    if( matchMuonTrigger(lv, d3pd.trig.trig_EF_trigmuonef_EF_mu18_tight_e7_medium1()) ) {
      flags |= TRIG_mu18_tight_e7_medium1;
    }

    // mu15
    if( matchMuonTrigger(lv, d3pd.trig.trig_EF_trigmuonef_EF_mu15()) ) {
      flags |= TRIG_mu15;
    }

    // assign the trigger flags for this muon
    m_muoTrigFlags[iMu] = flags;
  }
}
/*--------------------------------------------------------------------------------*/
bool SusyD3PDAna::matchMuonTrigger(const TLorentzVector* lv, vector<int>* passTrig)
{
  // loop over muon trigger features
  for(int iTrig=0; iTrig < d3pd.trig.trig_EF_trigmuonef_n(); iTrig++){

    // Check to see if this feature passed chain we want
    if(passTrig->at(iTrig)){

      // Loop over muon EF tracks
      TLorentzVector lvTrig;
      for(int iTrk=0; iTrk < d3pd.trig.trig_EF_trigmuonef_track_n()->at(iTrig); iTrk++){

        lvTrig.SetPtEtaPhiM( d3pd.trig.trig_EF_trigmuonef_track_CB_pt()->at(iTrig).at(iTrk),
                             d3pd.trig.trig_EF_trigmuonef_track_CB_eta()->at(iTrig).at(iTrk),
                             d3pd.trig.trig_EF_trigmuonef_track_CB_phi()->at(iTrig).at(iTrk),
                             0 );       // only eta and phi used to compute dR anyway
        // Require combined offline track...?
        if(!d3pd.trig.trig_EF_trigmuonef_track_CB_hasCB()->at(iTrig).at(iTrk)) continue;
        float dR = lv->DeltaR(lvTrig);
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
void SusyD3PDAna::matchTauTriggers()
{
  if(m_dbg>=5) cout << "matchTauTriggers" << endl;

  //int run = d3pd.evt.RunNumber();

  // loop over all pre taus
  for(uint i=0; i<m_preTaus.size(); i++){

    int iTau = m_preTaus[i];
    //const TLorentzVector* lv = & m_tauLVs[iTau];
    const TLorentzVector* lv = & m_susyObj.GetTauTLV(iTau);
    
    // trigger flags
    uint flags = 0;

    // tau20_medium1
    if( matchTauTrigger(lv, d3pd.trig.trig_EF_tau_EF_tau20_medium1()) ){
      flags |= TRIG_tau20_medium1;
    }
    // tau20Ti_medium1
    if( matchTauTrigger(lv, d3pd.trig.trig_EF_tau_EF_tau20Ti_medium1()) ){
      flags |= TRIG_tau20Ti_medium1;
    }
    // tau29Ti_medium1
    if( matchTauTrigger(lv, d3pd.trig.trig_EF_tau_EF_tau29Ti_medium1()) ){
      flags |= TRIG_tau29Ti_medium1;
    }
    // tau29Ti_medium1_tau20Ti_medium1
    if( matchTauTrigger(lv, d3pd.trig.trig_EF_tau_EF_tau29Ti_medium1_tau20Ti_medium1()) ){
      flags |= TRIG_tau29Ti_medium1_tau20Ti_medium1;
    }
    // tau20Ti_medium1_e18vh_medium1
    if( matchTauTrigger(lv, d3pd.trig.trig_EF_tau_EF_tau20Ti_medium1_e18vh_medium1()) ){
      flags |= TRIG_tau20Ti_medium1_e18vh_medium1;
    }
    // tau20_medium1_mu15
    if( matchTauTrigger(lv, d3pd.trig.trig_EF_tau_EF_tau20_medium1_mu15()) ){
      flags |= TRIG_tau20_medium1_mu15;
    }

    // assign the trigger flags for this tau
    m_tauTrigFlags[iTau] = flags;
  }
}
/*--------------------------------------------------------------------------------*/
bool SusyD3PDAna::matchTauTrigger(const TLorentzVector* lv, vector<int>* passTrig)
{
  // loop over tau trigger features
  for(int iTrig=0; iTrig < d3pd.trig.trig_EF_tau_n(); iTrig++){
    // Check to see if this feature passed chain we want
    if(passTrig->at(iTrig)){
      // Now, try to match offline tau to this online tau
      static TLorentzVector trigLV;
      trigLV.SetPtEtaPhiM(d3pd.trig.trig_EF_tau_pt()->at(iTrig), d3pd.trig.trig_EF_tau_eta()->at(iTrig), 
                          d3pd.trig.trig_EF_tau_phi()->at(iTrig), d3pd.trig.trig_EF_tau_m()->at(iTrig));
      float dR = lv->DeltaR(trigLV);
      if(dR < 0.15) return true;
    }
  }
  // matching failed
  return false;
}

/*--------------------------------------------------------------------------------*/
// Check event level cleaning cuts like GRL, LarError, etc.
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::checkEventCleaning()
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
void SusyD3PDAna::checkObjectCleaning()
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
bool SusyD3PDAna::passLarHoleVeto()
{
  TVector2 metVector = m_met.Vect().XYvector();
  vector<int> goodJets;
  // Do I still need these jets with no eta cut?
  // This only uses nominal jets...?  TODO
  vector<int> jets = get_jet_baseline( &d3pd.jet, &d3pd.vtx, &d3pd.evt, !m_isMC, m_susyObj, 
                                       20.*GeV, 9999999, SystErr::NONE, false, goodJets );
  return !check_jet_larhole(&d3pd.jet, jets, !m_isMC, m_susyObj, 180614, metVector, &m_fakeMetEst);
}
/*--------------------------------------------------------------------------------*/
// Pass tile hot spot veto
/*--------------------------------------------------------------------------------*/
bool SusyD3PDAna::passTileHotSpot()
{
  return !check_jet_tileHotSpot( &d3pd.jet, m_preJets, m_susyObj, !m_isMC, d3pd.evt.RunNumber() );
}
/*--------------------------------------------------------------------------------*/
// Pass bad jet cut
/*--------------------------------------------------------------------------------*/
bool SusyD3PDAna::passBadJet()
{
  return !IsBadJetEvent(&d3pd.jet, m_baseJets, 20.*GeV, m_susyObj);
}
/*--------------------------------------------------------------------------------*/
// Pass good vertex
/*--------------------------------------------------------------------------------*/
bool SusyD3PDAna::passGoodVtx()
{
  return PrimaryVertexCut(m_susyObj, &d3pd.vtx);
}
/*--------------------------------------------------------------------------------*/
// Pass tile trip
/*--------------------------------------------------------------------------------*/
bool SusyD3PDAna::passTileTrip()
{
  return !m_susyObj.IsTileTrip(d3pd.evt.RunNumber(), d3pd.evt.lbn(), d3pd.evt.EventNumber());
}
/*--------------------------------------------------------------------------------*/
// Pass bad muon veto
/*--------------------------------------------------------------------------------*/
bool SusyD3PDAna::passBadMuon()
{
  return !IsBadMuonEvent(m_susyObj, &d3pd.muo, m_preMuons, 0.2);
}
/*--------------------------------------------------------------------------------*/
// Pass cosmic veto
/*--------------------------------------------------------------------------------*/
bool SusyD3PDAna::passCosmic()
{
  return !IsCosmic(m_susyObj, &d3pd.muo, m_baseMuons, 1., 0.2);
}
/*--------------------------------------------------------------------------------*/
// Radiative b quark check for sherpa WW fix
/*--------------------------------------------------------------------------------*/
bool SusyD3PDAna::hasRadiativeBQuark(const vector<int>* pdg, const vector<int>* status)
{
  if(!pdg || !status || pdg->size()!=status->size()) return false;
  const vector<int>& p = *pdg;
  const vector<int>& s = *status;
  const int pdgB(5), statRad(3);
  for(size_t i=0; i<p.size(); ++i) if(abs(p[i])==pdgB && s[i]==statRad) return true;
  return false;
}

/*--------------------------------------------------------------------------------*/
// Get event weight, combine gen, pileup, xsec, and lumi weights
// Default weight uses ICHEP dataset, A-B14 lumi
// You can supply a different luminosity, but the pileup weights will still correspond to A-B14
/*--------------------------------------------------------------------------------*/
float SusyD3PDAna::getEventWeight(float lumi)
{
  if(!m_isMC) return 1;
  return d3pd.truth.event_weight() * getXsecWeight() * getPileupWeight() * lumi / m_sumw;
}
/*--------------------------------------------------------------------------------*/
float SusyD3PDAna::getEventWeightAtoB3()
{
  if(!m_isMC) return 1;
  return d3pd.truth.event_weight() * getXsecWeight() * getPileupWeightAB3() * LUMI_A_B3 / m_sumw;
}
/*--------------------------------------------------------------------------------*/
float SusyD3PDAna::getEventWeightAtoB()
{
  if(!m_isMC) return 1;
  return d3pd.truth.event_weight() * getXsecWeight() * getPileupWeightAB() * LUMI_A_B14 / m_sumw;
}

/*--------------------------------------------------------------------------------*/
// Cross section and lumi scaling
/*--------------------------------------------------------------------------------*/
float SusyD3PDAna::getXsecWeight()
{
  // Use user cross section if it has been set
  if(m_xsec > 0) return m_xsec;

  // Use SUSY cross section file
  int id = d3pd.truth.channel_number();
  if(m_xsecMap.find(id) == m_xsecMap.end()) {
    m_xsecMap[id] = m_susyXsec->process(id);
  }
  return m_xsecMap[id].xsect() * m_xsecMap[id].kfactor() * m_xsecMap[id].efficiency();
}

/*--------------------------------------------------------------------------------*/
// Luminosity normalization
/*--------------------------------------------------------------------------------*/
float SusyD3PDAna::getLumiWeight()
{ return m_lumi / m_sumw; }

/*--------------------------------------------------------------------------------*/
// Pileup reweighting
/*--------------------------------------------------------------------------------*/
float SusyD3PDAna::getPileupWeight()
{
  return m_pileup->GetCombinedWeight(d3pd.evt.RunNumber(), d3pd.truth.channel_number(), d3pd.evt.averageIntPerXing());
}
/*--------------------------------------------------------------------------------*/
float SusyD3PDAna::getPileupWeightAB3()
{
  return m_pileupAB3->GetCombinedWeight(d3pd.evt.RunNumber(), d3pd.truth.channel_number(), d3pd.evt.averageIntPerXing());
}
/*--------------------------------------------------------------------------------*/
float SusyD3PDAna::getPileupWeightAB()
{
  return m_pileupAB->GetCombinedWeight(d3pd.evt.RunNumber(), d3pd.truth.channel_number(), d3pd.evt.averageIntPerXing());
}
/*--------------------------------------------------------------------------------*/
float SusyD3PDAna::getPileupWeightAE()
{
  return m_pileupAE->GetCombinedWeight(d3pd.evt.RunNumber(), d3pd.truth.channel_number(), d3pd.evt.averageIntPerXing());
}

/*--------------------------------------------------------------------------------*/
// PDF reweighting of 7TeV -> 8TeV
/*--------------------------------------------------------------------------------*/
float SusyD3PDAna::getPDFWeight8TeV()
{
  #ifdef USEPDFTOOL
  float scale = d3pd.gen.pdf_scale()->at(0);
  float x1 = d3pd.gen.pdf_x1()->at(0);
  float x2 = d3pd.gen.pdf_x2()->at(0);
  int id1 = d3pd.gen.pdf_id1()->at(0);
  int id2 = d3pd.gen.pdf_id2()->at(0);

  // MultiLep function... Not working?
  //return scaleBeamEnergy(*m_pdfTool, 21000, d3pd.gen.pdf_scale()->at(0), d3pd.gen.pdf_x1()->at(0),
                         //d3pd.gen.pdf_x2()->at(0), d3pd.gen.pdf_id1()->at(0), d3pd.gen.pdf_id2()->at(0));
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
float SusyD3PDAna::getLepSF(const vector<LeptonInfo>& leptons)
{
  // TODO: incorporate systematics
  float lepSF = 1;

  if(m_isMC){
    // Loop over leptons
    for(uint iLep=0; iLep<leptons.size(); iLep++){
      const LeptonInfo* lep = & leptons[iLep];
      // Electrons
      if(lep->isElectron()){
        const ElectronElement* el = lep->getElectronElement();
        lepSF *= m_susyObj.GetSignalElecSF(el->cl_eta(), lep->lv()->Pt(), true, true, false);
      }
      // Muons
      else{
        lepSF *= m_susyObj.GetSignalMuonSF(lep->idx());
      }
    }
  }

  return lepSF;
}

/*--------------------------------------------------------------------------------*/
// BTag efficiency SF
// TODO: finish me!
/*--------------------------------------------------------------------------------*/
float SusyD3PDAna::getBTagSF(const vector<int>& jets)
{
  return 1;
}

/*--------------------------------------------------------------------------------*/
// Get the heavy flavor overlap removal decision
/*--------------------------------------------------------------------------------*/
int SusyD3PDAna::getHFORDecision()
{
  return m_hforTool.getDecision(d3pd.truth.channel_number(),
                                d3pd.truth.n(),
                                d3pd.truth.pt(),
                                d3pd.truth.eta(),
                                d3pd.truth.phi(),
                                d3pd.truth.m(),
                                d3pd.truth.pdgId(),
                                d3pd.truth.status(),
                                d3pd.truth.vx_barcode(),
                                d3pd.truth.parent_index(),
                                d3pd.truth.child_index(),
                                HforToolD3PD::ALL); //HforToolD3PD::DEFAULT
}

/*--------------------------------------------------------------------------------*/
// Set MET flavor via a string. Only a couple of options available so far
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::setMetFlavor(string metFlav)
{
  if(metFlav=="STVF") m_metFlavor = SUSYMet::STVF;
  else if(metFlav=="STVF_JVF") m_metFlavor = SUSYMet::STVF_JVF;
  else if(metFlav=="Default") m_metFlavor = SUSYMet::Default;
  else{
    cout << "SusyD3PDAna::setMetFlavor : ERROR : MET flavor " << metFlav 
         << " is not supported!" << endl;
    abort();
  }
}

/*--------------------------------------------------------------------------------*/
// Print event info
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::dumpEvent()
{
  cout << "Run " << setw(6) << d3pd.evt.RunNumber()
       << " Event " << setw(7) << d3pd.evt.EventNumber()
       << " Stream " << streamName(m_stream);
  if(m_isMC){
    cout << " MCID " << setw(6) << d3pd.truth.channel_number()
         << " weight " << d3pd.truth.event_weight();
  }
  cout << endl;
}

/*--------------------------------------------------------------------------------*/
// Print baseline objects
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::dumpBaselineObjects()
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
      const TLorentzVector* lv = & m_susyObj.GetElecTLV(iEl);
      const ElectronElement* ele = & d3pd.ele[iEl];
      cout << "  El : " << fixed
           << " q " << setw(2) << (int) ele->charge()
           << " pt " << setw(6) << lv->Pt()/GeV
           << " eta " << setw(5) << lv->Eta()
           << " phi " << setw(5) << lv->Phi();
      if(m_isMC) cout << " type " << setw(2) << ele->type() << " origin " << setw(2) << ele->origin();
      cout << endl;
    }
  }
  if(nMu){
    cout << "Baseline muons" << endl;
    for(uint i=0; i < nMu; i++){
      int iMu = m_baseMuons[i];
      const TLorentzVector* lv = & m_susyObj.GetMuonTLV(iMu);
      const MuonElement* muo = & d3pd.muo[iMu];
      cout << "  Mu : " << fixed
           << " q " << setw(2) << (int) muo->charge()
           << " pt " << setw(6) << lv->Pt()/GeV
           << " eta " << setw(5) << lv->Eta()
           << " phi " << setw(5) << lv->Phi();
      if(m_isMC) cout << " type " << setw(2) << muo->type() << " origin " << setw(2) << muo->origin();
      cout << endl;
    }
  }
  if(nJet){
    cout << "Baseline jets" << endl;
    for(uint i=0; i < nJet; i++){
      int iJet = m_baseJets[i];
      const TLorentzVector* lv = & m_susyObj.GetJetTLV(iJet);
      const JetElement* jet = & d3pd.jet[iJet];
      cout << "  Jet : " << fixed
           << " pt " << setw(6) << lv->Pt()/GeV
           << " eta " << setw(5) << lv->Eta()
           << " phi " << setw(5) << lv->Phi()
           << " mv1 " << jet->flavor_weight_MV1();
      cout << endl;
    }
  }
  cout.precision(6);
  cout.unsetf(ios_base::fixed);
}

/*--------------------------------------------------------------------------------*/
// Print signal objects
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::dumpSignalObjects()
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
      const TLorentzVector* lv = & m_susyObj.GetElecTLV(iEl);
      const ElectronElement* ele = & d3pd.ele[iEl];
      cout << "  El : " << fixed
           << " q " << setw(2) << (int) ele->charge()
           << " pt " << setw(6) << lv->Pt()/GeV
           << " eta " << setw(5) << lv->Eta()
           << " phi " << setw(5) << lv->Phi();
      if(m_isMC) cout << " type " << setw(2) << ele->type() << " origin " << setw(2) << ele->origin();
      cout << endl;
    }
  }
  if(nMu){
    cout << "Signal muons" << endl;
    for(uint i=0; i < nMu; i++){
      int iMu = m_sigMuons[i];
      const TLorentzVector* lv = & m_susyObj.GetMuonTLV(iMu);
      const MuonElement* muo = & d3pd.muo[iMu];
      cout << "  Mu : " << fixed
           << " q " << setw(2) << (int) muo->charge()
           << " pt " << setw(6) << lv->Pt()/GeV
           << " eta " << setw(5) << lv->Eta()
           << " phi " << setw(5) << lv->Phi();
      if(m_isMC) cout << " type " << setw(2) << muo->type() << " origin " << setw(2) << muo->origin();
      cout << endl;
    }
  }
  if(nJet){
    cout << "Signal jets" << endl;
    for(uint i=0; i < nJet; i++){
      int iJet = m_sigJets[i];
      const TLorentzVector* lv = & m_susyObj.GetJetTLV(iJet);
      const JetElement* jet = & d3pd.jet[iJet];
      cout << "  Jet : " << fixed
           << " pt " << setw(6) << lv->Pt()/GeV
           << " eta " << setw(5) << lv->Eta()
           << " phi " << setw(5) << lv->Phi()
           << " mv1 " << jet->flavor_weight_MV1();
      cout << endl;
    }
  }
  cout.precision(6);
  cout.unsetf(ios_base::fixed);
}

#undef GeV
