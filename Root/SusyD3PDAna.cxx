#include "TSystem.h"
#include "SusyCommon/SusyD3PDAna.h"
#include "MultiLep/ElectronTools.h"
#include "MultiLep/MuonTools.h"
#include "MultiLep/JetTools.h"
#include "MultiLep/CutflowTools.h"

using namespace std;

/*--------------------------------------------------------------------------------*/
// SusyD3PDAna Constructor
/*--------------------------------------------------------------------------------*/
SusyD3PDAna::SusyD3PDAna() : 
        m_sample(""),
        m_lumi(4700),
        m_sumw(1),
	m_xsec(-1),
	m_sys(false),
        m_pileup(0),
        m_susyXsec(0)
{
}
/*--------------------------------------------------------------------------------*/
// Destructor
/*--------------------------------------------------------------------------------*/
SusyD3PDAna::~SusyD3PDAna()
{
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
  else if(m_sample.Contains("muons", TString::kIgnoreCase))  m_stream = Stream_Muons;
  else if(m_sample.Contains("egamma", TString::kIgnoreCase)) m_stream = Stream_Egamma;
  else m_stream = Stream_Unknown;

  if(m_isMC) cout << "Processing as MC"   << endl;
  else       cout << "Processing as DATA" << endl;

  cout << "DataStream: " << streamName(m_stream) << endl;

  // Setup SUSYTools
  m_susyObj.initialize(!m_isMC);
  m_fakeMetEst.initialize("$ROOTCOREDIR/data/MultiLep/fest_periodF_v1.root");

  // SUSY cross sections
  if(m_isMC){
    string xsecFileName  = gSystem->ExpandPathName("$ROOTCOREDIR/data/SUSYTools/susy_crosssections.txt");
    m_susyXsec = new SUSY::CrossSectionDB(xsecFileName);
  }

  // GRL
  if(!m_isMC){
    Root::TGoodRunsListReader* grlReader = new Root::TGoodRunsListReader();
    grlReader->AddXMLFile(m_grlFileName);
    grlReader->Interpret();
    m_grl = grlReader->GetMergedGoodRunsList();
    delete grlReader;
  }

  // Pileup reweighting
  if(m_isMC){
    m_pileup = new Root::TPileupReweighting("PileupReweighting");
    m_pileup->AddLumiCalcFile("$ROOTCOREDIR/data/MultiLep/ilumicalc_histograms_None_178044-191933.root");
    m_pileup->AddConfigFile("$ROOTCOREDIR/data/MultiLep/mc11b_defaults.root");
    m_pileup->SetUnrepresentedDataAction(2);
    int pileupError = m_pileup->Initialize();

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
  }
}

/*--------------------------------------------------------------------------------*/
// Baseline object selection
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::selectBaselineObjects(SYSTEMATIC sys)
{
  if(m_dbg) cout << "selectBaselineObjects" << endl;
  vector<int> goodJets;  // What the hell is this??
  float mu = d3pd.evt.averageIntPerXing();

  // Handle Systematic
  int ees = 0, eer = 0;
  string musys = "";
  JetErr::Syste jetsys = JetErr::NONE;
  if(sys == NOM);                                  // No need to check needlessly
  else if(sys == EES_UP) ees = 1;                  // E scale up
  else if(sys == EES_DN) ees = 2;                  // E scale down
  else if(sys == EER_UP) eer = 1;                  // E smear up
  else if(sys == EER_DN) eer = 2;                  // E smear down
  else if(sys == MS_UP ) musys = "MSUP";           // MS scale up
  else if(sys == MS_DN ) musys = "MSLOW";          // MS scale down
  else if(sys == ID_UP ) musys = "IDUP";           // ID scale up
  else if(sys == ID_DN ) musys = "IDLOW";          // ID scale down
  else if(sys == JES_UP) jetsys = JetErr::JESUP;   // JES up
  else if(sys == JES_DN) jetsys = JetErr::JESDOWN; // JES down
  else if(sys == JER)    jetsys = JetErr::JER;     // JER (gaussian)

  // Preselection
  m_preElectrons = get_electrons_baseline( &d3pd.ele, !m_isMC, d3pd.evt.RunNumber(), m_susyObj, 10.*GeV, 2.47, ees, eer, false );
  m_preMuons     = get_muons_baseline( &d3pd.muo, !m_isMC, m_susyObj, 10.*GeV, 2.4, musys );
  m_preJets      = get_jet_baseline( &d3pd.jet, !m_isMC, m_susyObj, 20.*GeV, 4.9, jetsys, false, goodJets, mu, &d3pd.vtx );
  
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
  m_baseElectrons = overlap_removal(m_susyObj, &d3pd.ele, m_preElectrons, &d3pd.ele, m_preElectrons, 0.1, 1);

  // jet-e overlap removal
  m_baseJets      = overlap_removal(m_susyObj, &d3pd.jet, m_preJets, &d3pd.ele, m_baseElectrons, 0.2, 0);

  // e-jet overlap removal
  m_baseElectrons = overlap_removal(m_susyObj, &d3pd.ele, m_baseElectrons, &d3pd.jet, m_baseJets, 0.4, 0);

  // m-jet overlap removal
  m_baseMuons     = overlap_removal(m_susyObj, &d3pd.muo, m_preMuons, &d3pd.jet, m_baseJets, 0.4, 0);

  // e-m overlap removal
  vector<int> copyElectrons = m_baseElectrons;
  m_baseElectrons = overlap_removal(m_susyObj, &d3pd.ele, m_baseElectrons, &d3pd.muo, m_baseMuons, 0.1, 0);
  //m_baseMuons     = overlap_removal(m_susyObj, &d3pd.muo, m_baseMuons, &d3pd.ele, m_baseElectrons, 0.1, 0);
  m_baseMuons     = overlap_removal(m_susyObj, &d3pd.muo, m_baseMuons, &d3pd.ele, copyElectrons, 0.1, 0);

  // Moved to SusyNtTools to be called during ana
  // remove SFOS lepton pairs with Mll < 20 GeV
  //m_baseElectrons = RemoveSFOSPair(m_susyObj, &d3pd.ele, m_baseElectrons, 20.*GeV);
  //m_baseMuons     = RemoveSFOSPair(m_susyObj, &d3pd.muo, m_baseMuons,     20.*GeV);


}

/*--------------------------------------------------------------------------------*/
// Signal object selection - do baseline selection first!
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::selectSignalObjects()
{
  if(m_dbg) cout << "selectSignalObjects" << endl;
  // TODO: make these functions more symmetric
  m_sigElectrons = get_electrons_signal(&d3pd.ele, m_baseElectrons, m_susyObj, 10.*GeV);
  m_sigMuons     = get_muons_signal(&d3pd.muo, m_susyObj, m_baseMuons, 10.*GeV, 1.8*GeV);
  m_sigJets      = get_jet_signal(&d3pd.jet, m_susyObj, m_baseJets, 20.*GeV, 2.5, 0.75);

  // combine leptons
  m_sigLeptons   = buildLeptonInfos(&d3pd.ele, m_sigElectrons, &d3pd.muo, m_sigMuons, m_susyObj);
}

/*--------------------------------------------------------------------------------*/
// Build MissingEt
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::buildMet(SYSTEMATIC sys)
{
  if(m_dbg) cout << "buildMet" << endl;
 
  // Need the proper jet systematic for building systematic
  JetErr::Syste jetsys = JetErr::NONE;     // Nominal
  if(sys == NOM);
  else if(sys == JES_UP) jetsys = JetErr::JESUP;   // JES up
  else if(sys == JES_DN) jetsys = JetErr::JESDOWN; // JES down
  else if(sys == JER)    jetsys = JetErr::JER;     // JER (gaussian)
  

  // Need ALL electrons in order to calculate the MET
  // Actually, I see common code uses all electrons that have lv.Pt() != 0
  // That's fine though because SUSYObjDef specifically fills for electrons that
  // should enter the RefEle term
  vector<int> allElectrons = get_electrons_all(&d3pd.ele, m_susyObj);
  // MET uses muons before overlap removal
  //TVector2 metVector = GetMetVector(&d3pd.jet, m_susyObj, &d3pd.muo, &d3pd.ele, &d3pd.met, m_baseMuons, m_baseElectrons, allElectrons, JetErr::NONE);
  TVector2 metVector = GetMetVector(&d3pd.jet, m_susyObj, &d3pd.muo, &d3pd.ele, &d3pd.met, m_preMuons, m_baseElectrons, allElectrons, jetsys);
  m_met.SetPxPyPzE(metVector.X(), metVector.Y(), 0, metVector.Mod());
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
  m_evtFlag = 0;
}

/*--------------------------------------------------------------------------------*/
// Electron trigger matching
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::matchElectronTriggers(bool disregardPt)
{
  if(m_dbg) cout << "matchElectronTriggers" << endl;
  int run = d3pd.evt.RunNumber();

  // loop over all pre electrons
  for(uint i=0; i<m_preElectrons.size(); i++){
    int iEl = m_preElectrons[i];
    const TLorentzVector* lv = & m_susyObj.GetElecTLV(iEl);
    
    // trigger flags
    uint flags = 0;

    if(disregardPt || lv->Pt() > 25.*GeV){
      // e20_medium
      if( m_isMC || (run<186873 && matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_e20_medium())) ){
        flags |= TRIG_e20_medium;
      }
      // e22_medium
      if( m_isMC || (run<188902 && matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_e22_medium())) ){
        flags |= TRIG_e22_medium;
      }
      // e22vh_medium1
      if( m_isMC || (run>188901 && matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_e22vh_medium1())) ){
        flags |= TRIG_e22vh_medium1;
      }
    }
    if(disregardPt || lv->Pt() > 17.*GeV){
      // 2e12_medium
      if( m_isMC || (run<186873 && matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_2e12_medium())) ){
        flags |= TRIG_2e12_medium;
      }
      // 2e12T_medium
      if( m_isMC || (run>186873 && run<188902 && matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_2e12T_medium())) ){
        flags |= TRIG_2e12T_medium;
      }
      // 2e12Tvh_medium
      if( m_isMC || (run>188901 && matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_2e12Tvh_medium())) ){
        flags |= TRIG_2e12Tvh_medium;
      }
    }
    if(disregardPt || lv->Pt() > 15.*GeV){
      // e10_medium_mu6
      if( m_isMC || (run>185353 && matchElectronTrigger(lv, d3pd.trig.trig_EF_el_EF_e10_medium_mu6())) ){
        flags |= TRIG_e10_medium_mu6;
      }
    }

    // assign the flags in the map
    m_eleTrigFlags[iEl] = flags;
  }
}
/*--------------------------------------------------------------------------------*/
//bool SusyD3PDAna::matchElectronTrigger(float eta, float phi, vector<int>* trigBools)
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
void SusyD3PDAna::matchMuonTriggers(bool disregardPt)
{
  if(m_dbg) cout << "matchMuonTriggers" << endl;

  // New prescription!
  int run = d3pd.evt.RunNumber();
  // loop over all pre muons
  for(uint i=0; i<m_preMuons.size(); i++){

    int iMu = m_preMuons[i];
    const TLorentzVector* lv = & m_susyObj.GetMuonTLV(iMu);
    
    // trigger flags
    uint flags = 0;

    if(disregardPt || lv->Pt() > 20.*GeV){
      // mu18
      if( m_isMC || (run<186516 && matchMuonTrigger(lv, d3pd.trig.trig_EF_trigmuonef_EF_mu18()))) {
        flags |= TRIG_mu18;
      }
      // mu18_medium
      if( m_isMC || (run>=186516 && matchMuonTrigger(lv, d3pd.trig.trig_EF_trigmuonef_EF_mu18_medium()))) {
        flags |= TRIG_mu18_medium;
      }
    }
    if(disregardPt || lv->Pt() > 10.*GeV){
      // 2mu10_loose
      if( m_isMC || matchMuonTrigger(lv, d3pd.trig.trig_EF_trigmuonef_EF_2mu10_loose())) {
        flags |= TRIG_2mu10_loose;
      }
    }
    if(disregardPt || lv->Pt() > 8.*GeV){
      // e10_medium_mu6
      if( m_isMC || matchMuonTrigger(lv, d3pd.trig.trig_EF_trigmuonef_EF_mu6()) ) {
        flags |= TRIG_e10_medium_mu6;
      }
    }

    // assign the flags for this muon
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
// Check event level cuts, like LArHole veto, badJet, etc.
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::evtCheck(){

  // Lar Hole Veto
  if(passLarHoleVeto())
    m_evtFlag |= PASS_LAr;

  // Bad Jet
  if(passBadJet())
    m_evtFlag |= PASS_BadJet;
  
  // Bad Muon
  if(passBadMuon())
    m_evtFlag |= PASS_BadMuon;
  
  // Cosmic muon check
  if(passCosmic())
    m_evtFlag |= PASS_Cosmic;

  // Now store the pass all
  if( (m_evtFlag & PASS_LAr) &&
      (m_evtFlag & PASS_BadJet) &&
      (m_evtFlag & PASS_BadMuon) &&
      (m_evtFlag & PASS_Cosmic) )
    m_evtFlag |= PASS_Event;

}

/*--------------------------------------------------------------------------------*/
// Pass Lar hole veto
// Prior to calling this, need jet and MET selection
/*--------------------------------------------------------------------------------*/
bool SusyD3PDAna::passLarHoleVeto()
{
  TVector2 metVector = m_met.Vect().XYvector();
  vector<int> goodJets;
  float mu = d3pd.evt.averageIntPerXing();
  // Do I still need these jets with no eta cut?
  vector<int> jets = get_jet_baseline( &d3pd.jet, !m_isMC, m_susyObj, 20.*GeV, 9999999, JetErr::NONE, false, goodJets,mu,&d3pd.vtx );
  return !check_jet_larhole(&d3pd.jet, jets, !m_isMC, m_susyObj, 180614, metVector, &m_fakeMetEst);
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
// Method for quick debuggin'
/*--------------------------------------------------------------------------------*/
void SusyD3PDAna::dump(){

  // Right now I need to debug the jets, so that is what I will dump
  for(uint i=0; i<m_preJets.size(); ++i){
    int idx = m_preJets.at(i);
    TLorentzVector tlv = m_susyObj.GetJetTLV( idx );
    cout<<"Jet index: "<<idx<<" Pt: "<<tlv.Pt()/GeV<<" Eta: "<<tlv.Eta()<<" Phi: "<<tlv.Phi()<<endl;
  }

}
