#include "egammaAnalysisUtils/CaloIsoCorrection.h"

//#include "TauCorrections/TauCorrections.h"
#include "TauCorrUncert/TauSF.h"
//#include "SUSYTools/MV1.h"

#include "SusyCommon/SusyNtMaker.h"
#include "SusyCommon/TruthTools.h"
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/WhTruthExtractor.h"
#include "SusyNtuple/mc_truth_utils.h"
#include "SusyNtuple/RecoTruthClassification.h"

#include "ElectronEfficiencyCorrection/TElectronEfficiencyCorrectionTool.h"

#include "xAODPrimitives/IsolationType.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODEgamma/EgammaxAODHelpers.h"

// Amg include
#include "EventPrimitives/EventPrimitivesHelpers.h"


#include <algorithm> // max_element
#include <iomanip> // setw
#include <sstream> // std::ostringstream

using namespace std;
namespace smc =susy::mc;

using susy::SusyNtMaker;

#define GeV 1000.
const double MeV2GeV=1.0e-3;

//----------------------------------------------------------
SusyNtMaker::SusyNtMaker() :
    m_outTreeFile(NULL),
    m_outTree(NULL),
    m_susyNt(0),
    m_fillNt(true),
    m_filter(true),
    m_nLepFilter(0),
    m_nLepTauFilter(2),
    m_filterTrigger(false),
    m_saveContTaus(false),
    m_isWhSample(false),
    m_hDecay(0),
    m_hasSusyProp(false),
    h_rawCutFlow(NULL),
    h_genCutFlow(NULL),
    m_cutstageCounters(SusyNtMaker::cutflowLabels().size(), 0)
{
  n_base_ele=0;
  n_base_muo=0;
  n_base_tau=0;
  n_base_jet=0;
  n_sig_ele=0;
  n_sig_muo=0;
  n_sig_tau=0;
  n_sig_jet=0;

  //AT-2014-11-05: Correct place to put this ?
  CP::SystematicCode::enableFailure();
}
//----------------------------------------------------------
SusyNtMaker::~SusyNtMaker()
{
}
//----------------------------------------------------------
void SusyNtMaker::SlaveBegin(TTree* tree)
{
  XaodAnalysis::SlaveBegin(tree);
  if(m_dbg)
      cout<<"SusyNtMaker::SlaveBegin"<<endl;
  if(m_fillNt)
      initializeOuputTree();
  m_isWhSample = guessWhetherIsWhSample(m_sample);
  initializeCutflowHistograms();

  m_timer.Start();
}
//----------------------------------------------------------
const std::vector< std::string > SusyNtMaker::cutflowLabels()
{
    vector<string> labels;
    labels.push_back("Initial"        );
    labels.push_back("SusyProp Veto"  );
    labels.push_back("GRL"            );
    labels.push_back("LAr Error"      );
    labels.push_back("Tile Error"     );
    labels.push_back("TTC Veto"       );
    labels.push_back("Good Vertex"    );
    labels.push_back("Buggy WWSherpa" );
    labels.push_back("Hot Spot"       );
    labels.push_back("Bad Jet"        );
    labels.push_back("Bad Muon"       );
    labels.push_back("Cosmic"         );
    labels.push_back(">=1 lep"        );
    labels.push_back(">=2 lep"        );
    labels.push_back("==3 lep"        );
    return labels;
}
//----------------------------------------------------------
TH1F* SusyNtMaker::makeCutFlow(const char* name, const char* title)
{
    vector<string> labels = SusyNtMaker::cutflowLabels();
    int nCuts = labels.size();
    TH1F* h = new TH1F(name, title, nCuts, 0., static_cast<float>(nCuts));
    for(int iCut=0; iCut<nCuts; ++iCut)
        h->GetXaxis()->SetBinLabel(iCut+1, labels[iCut].c_str());
    return h;
}
//----------------------------------------------------------
TH1F* SusyNtMaker::getProcCutFlow(int signalProcess)
{
  // Look for it on the map
  map<int,TH1F*>::const_iterator it = m_procCutFlows.find(signalProcess);
  // New process
  if(it == m_procCutFlows.end()){
    stringstream stream;
    stream << signalProcess;
    string name = "procCutFlow" + stream.str();
    return m_procCutFlows[signalProcess] = makeCutFlow(name.c_str(),
                                                       (name+";Cuts;Events").c_str());
  }
  // Already saved process
  else{
    return it->second;
  }
}
//----------------------------------------------------------
Bool_t SusyNtMaker::Process(Long64_t entry)
{
  static Long64_t chainEntry = -1;
  chainEntry++;
  m_event.getEntry(chainEntry); // DG 2014-09-19 TEvent wants the chain entry, not the tree entry (?)
  clearOutputObjects();
  retrieveCollections();
  
  const xAOD::EventInfo* eventinfo = XaodAnalysis::xaodEventInfo();

  if(!m_flagsHaveBeenChecked) {
      m_flagsAreConsistent = runningOptionsAreValid();
      m_flagsHaveBeenChecked=true;
      if(!m_flagsAreConsistent) {
          cout<<"ERROR: Inconsistent options. Stopping here."<<endl;
          abort();
      }
  }

  if(m_dbg || chainEntry%5000==0)
  {
    cout << "***********************************************************" << endl;
    cout << "**** Processing entry " << setw(6) << chainEntry
         << " run " << setw(6) << eventinfo->runNumber()
         << " event " << setw(7) << eventinfo->eventNumber() << " ****" << endl;
    cout << "***********************************************************" << endl;
  }

  if(selectEvent() && m_fillNt){
      fillNtVars();
      if(m_isMC && m_sys) doSystematic();
      int bytes = m_outTree->Fill();
      if(bytes==-1){
          cout << "SusyNtMaker ERROR filling tree!  Abort!" << endl;
          abort();
      }
  }
  deleteShallowCopies();
  clearContainerPointers();
  return kTRUE;
}
//----------------------------------------------------------
void SusyNtMaker::Terminate()
{
  XaodAnalysis::Terminate();
  m_timer.Stop();
  if(m_dbg) cout<<"SusyNtMaker::Terminate"<<endl;
  cout<<counterSummary()<<endl;
  cout<<timerSummary()<<endl;
  if(m_fillNt) saveOutputTree();
}
//----------------------------------------------------------
bool isSimplifiedModel(const TString &sampleName)
{
    return sampleName.Contains("simplifiedModel");
}
//----------------------------------------------------------
bool SusyNtMaker::selectEvent()
{
  if(m_dbg>=5) cout << "selectEvent" << endl;
  clearOutputObjects();
  m_susyNt.clear();
  bool pass_event_and_object_sel = (passEventlevelSelection() && passObjectlevelSelection());
  bool keep_all_events(!m_filter);
  return (pass_event_and_object_sel || keep_all_events);
}
//----------------------------------------------------------
void SusyNtMaker::fillNtVars()
{
  fillEventVars();
  fillElectronVars();
  fillMuonVars();
  fillTauVars();
  fillJetVars();
  fillMetVars();
  fillPhotonVars();
  if(m_isMC && getSelectTruthObjects() ) {
    fillTruthParticleVars();
    fillTruthJetVars();
    fillTruthMetVars();
  }
}
//----------------------------------------------------------
void SusyNtMaker::fillEventVars()
{
  if(m_dbg>=5) cout << "fillEventVars" << endl;
  Susy::Event* evt = m_susyNt.evt();
  const xAOD::EventInfo* eventinfo = XaodAnalysis::xaodEventInfo();

  evt->run              = eventinfo->runNumber();
  evt->event            = eventinfo->eventNumber();
  evt->lb               = eventinfo->lumiBlock();
  evt->stream           = m_stream;

  evt->isMC             = m_isMC;
  evt->mcChannel        = m_isMC? eventinfo->mcChannelNumber() : 0;
  evt->w                = m_isMC? eventinfo->mcEventWeight()   : 1; // \todo DG this now has an arg, is 0 the right one?
  evt->nVtx             = getNumGoodVtx();
  evt->avgMu            = eventinfo->averageInteractionsPerCrossing();

  evt->hfor             = m_isMC? getHFORDecision() : -1;

  // SUSY final state
  evt->susyFinalState   = m_susyFinalState;
  // \todo for later (DG could become obsolete?)
  // evt->susySpartId1     = m_event.SUSY.Spart1_pdgId.IsAvailable()? m_event.SUSY.Spart1_pdgId() : 0;
  // evt->susySpartId2     = m_event.SUSY.Spart2_pdgId.IsAvailable()? m_event.SUSY.Spart2_pdgId() : 0;

  float mZ = -1.0, mZtruthMax = 40.0;
  if(m_isMC){
    // int dsid = m_event.eventinfo.mc_channel_number();
    // if(IsAlpgenLowMass(dsid) || IsAlpgenPythiaZll(dsid)) mZ = MllForAlpgen(&m_event.mc);
    // else if(IsSherpaZll(dsid)) mZ = MllForSherpa(&m_event.mc);
  }
  evt->mllMcTruth       = mZ;
  evt->passMllForAlpgen = m_isMC ? (mZ < mZtruthMax) : true;
  evt->hDecay           = m_hDecay;
  evt->eventWithSusyProp= m_hasSusyProp;
  evt->trigFlags        = m_evtTrigFlags;

  evt->wPileup          = m_isMC? getPileupWeight(eventinfo) : 1;
  evt->wPileup_up       = m_isMC? getPileupWeightUp() : 1;
  evt->wPileup_dn       = m_isMC? getPileupWeightDown() : 1;
  evt->xsec             = m_isMC? getXsecWeight() : 1;
  evt->errXsec          = m_isMC? m_errXsec : 1;
  evt->sumw             = m_isMC? m_sumw : 1;

  if(m_isMC){
      xAOD::TruthEventContainer::const_iterator truthE_itr = xaodTruthEvent()->begin();
      // ( *truthE_itr )->pdfInfoParameter(evt->pdf_id1   , xAOD::TruthEvent::PDGID1); // not available for some samples
      // ( *truthE_itr )->pdfInfoParameter(evt->pdf_id2   , xAOD::TruthEvent::PDGID2);
      // ( *truthE_itr )->pdfInfoParameter(evt->pdf_x1    , xAOD::TruthEvent::X1);
      // ( *truthE_itr )->pdfInfoParameter(evt->pdf_x2    , xAOD::TruthEvent::X2);
      // ( *truthE_itr )->pdfInfoParameter(evt->pdf_scale , xAOD::TruthEvent::SCALE);
      // DG what are these two?
      //( *truthE_itr )->pdfInfoParameter(evt->pdf_x1   , xAOD::TruthEvent::x1);
      //( *truthE_itr )->pdfInfoParameter(evt->pdf_x2   , xAOD::TruthEvent::x2);
  }
  evt->pdfSF            = m_isMC? getPDFWeight8TeV() : 1;
  m_susyNt.evt()->cutFlags[NtSys::NOM] = m_cutFlags;
}
//----------------------------------------------------------
void SusyNtMaker::fillElectronVars()
{
    if(m_dbg>=5) cout<<"fillElectronVars"<<endl;
    xAOD::ElectronContainer* electrons = XaodAnalysis::xaodElectrons();
    for(auto &i : m_preElectrons){
        storeElectron(*(electrons->at(i)));
    }
}
//----------------------------------------------------------
void SusyNtMaker::fillMuonVars()
{
    if(m_dbg>=5) cout<<"fillMuonVars"<<endl;
    xAOD::MuonContainer* muons = XaodAnalysis::xaodMuons();
    for(auto &i : m_preMuons){
        storeMuon(*(muons->at(i)));
    }
}
//----------------------------------------------------------
void SusyNtMaker::fillJetVars()
{
    if(m_dbg>=5) cout<<"fillJetVars"<<endl;
    xAOD::JetContainer* jets = XaodAnalysis::xaodJets();
    for(auto &i : m_preJets){
        storeJet(*(jets->at(i)));
    }
}
//----------------------------------------------------------
void SusyNtMaker::fillTauVars()
{
    if(m_dbg>=5) cout<<"fillTauVars"<<endl;
    xAOD::TauJetContainer* taus =  XaodAnalysis::xaodTaus();
    // vector<int>& saveTaus = m_saveContTaus? m_contTaus : m_preTaus; // container taus are meant to be used only for background estimates?
    vector<int>& saveTaus = m_preTaus;
    for(auto &i : saveTaus){
        storeTau(*(taus->at(i)));
    }
}
//----------------------------------------------------------
void SusyNtMaker::fillPhotonVars()
{
    if(m_dbg>=5) cout<<"fillPhotonVars"<<endl;
    const xAOD::PhotonContainer* photons = XaodAnalysis::xaodPhotons();
    for(auto &i : m_sigPhotons){
        storePhoton(*(photons->at(i)));
    }
}
//----------------------------------------------------------
bool isMcAtNloTtbar(const int &channel) { return channel==105200; }
//----------------------------------------------------------
void SusyNtMaker::fillTruthParticleVars()
{
    if(m_dbg>=5) cout<<"fillTruthParticleVars"<<endl;
    // DG-2014-08-29 todo
  // // Retrieve indicies -- should go elsewhere
  // m_truParticles        = m_recoTruthMatch.LepFromHS_McIdx();
  // vector<int> truthTaus = m_recoTruthMatch.TauFromHS_McIdx();
  // m_truParticles.insert( m_truParticles.end(), truthTaus.begin(), truthTaus.end() );
  // if(m_isMC && isMcAtNloTtbar(m_event.eventinfo.mc_channel_number())){
  //   vector<int> ttbarPart(WhTruthExtractor::ttbarMcAtNloParticles(m_event.mc.pdgId(),
  //                                                                 m_event.mc.child_index()));
  //   m_truParticles.insert(m_truParticles.end(), ttbarPart.begin(), ttbarPart.end());
  // }
    const xAOD::TruthParticleContainer* particles = xaodTruthParticles();
    for(auto &i : m_truParticles){
        storeTruthParticle(*(particles->at(i)));
    }
}

//----------------------------------------------------------
void SusyNtMaker::storeElectron(const xAOD::Electron &in)
{
    Susy::Electron out;
    double pt(in.pt()*MeV2GeV), eta(in.eta()), phi(in.phi()), m(in.m()*MeV2GeV);
    out.SetPtEtaPhiM(pt, eta, phi, m);
    out.pt  = pt;
    out.eta = eta;
    out.phi = phi;
    out.m   = m;
    out.isBaseline = in.auxdata< char >("baseline");
    out.isSignal = in.auxdata< char >("signal");
    out.q   = in.charge();
    bool all_available=true;
    
    // IsEM quality flags - no need to recalculate them
    all_available &= in.passSelection(out.mediumPP,"Medium");
    all_available &= in.passSelection(out.tightPP,"Tight"); 
    all_available &= in.passSelection(out.looseLLH,"LooseLLH");
    all_available &= in.passSelection(out.mediumLLH,"MediumLLH");
    all_available &= in.passSelection(out.veryTightLLH,"VeryTightLLH");

    //Isolations
    all_available &= in.isolationValue(out.etcone20, xAOD::Iso::etcone20); 
    all_available &= in.isolationValue(out.topoEtcone30Corr, xAOD::Iso::topoetcone30); 
    all_available &= in.isolationValue(out.ptcone20, xAOD::Iso::ptcone20);
    all_available &= in.isolationValue(out.ptcone30, xAOD::Iso::ptcone30);
    out.etcone20 *= MeV2GeV;
    out.topoEtcone30Corr *= MeV2GeV;
    out.ptcone20 *= MeV2GeV;
    out.ptcone30 *= MeV2GeV;

    if(m_isMC){
      //Store the SF of the tightest ID
      bool recoSF=true;
      bool idSF=true;
      bool trigSF=false;
      //AT 2014-10-29: To be updated once SusyTools function return also the error.
      if(eleIsOfType(in, Tight))
	out.effSF = m_susyObj[Tight]->GetSignalElecSF(in, recoSF, idSF, trigSF);
      else if(eleIsOfType(in, Medium))
	out.effSF = m_susyObj[Medium]->GetSignalElecSF(in, recoSF, idSF, trigSF);
      
      if(eleIsOfType(in, VeryTightLLH))
	out.effSF_LLH = m_susyObj[VeryTightLLH]->GetSignalElecSF(in, recoSF, idSF, trigSF);
      else if(eleIsOfType(in, MediumLLH))
	out.effSF_LLH = m_susyObj[MediumLLH]->GetSignalElecSF(in, recoSF, idSF, trigSF);	 
      else if(eleIsOfType(in, LooseLLH))
	out.effSF_LLH = m_susyObj[LooseLLH]->GetSignalElecSF(in, recoSF, idSF, trigSF);

      if(m_dbg) cout << "AT: susyTool electron SF " << out.effSF << " LLH " << out.effSF_LLH << endl;
      /*
      const Root::TResult &result =  m_electronEfficiencySFTool->calculate(in);
      out.effSF    = result.getScaleFactor();
      out.errEffSF = result.getTotalUncertainty();
      if(m_dbg) cout << "AT: electron SF " << out.effSF << " " << out.errEffSF << endl;
      */
    
      out.mcType   = xAOD::EgammaHelpers::getParticleTruthType(&in);
      out.mcOrigin = xAOD::EgammaHelpers::getParticleTruthOrigin(&in);    
      const xAOD::TruthParticle* truthEle = xAOD::EgammaHelpers::getTruthParticle(&in); //AT 10/12/14: Always false ????
      out.matched2TruthLepton   = truthEle ? true : false;
      int matchedPdgId = truthEle ? truthEle->pdgId() : -999;
      out.truthType  = isFakeLepton(out.mcOrigin, out.mcType, matchedPdgId); 
      
      //AT 12/09/14: Need to get this from eGamma - at some point
      out.isChargeFlip          = m_isMC? m_recoTruthMatch.isChargeFlip(out, out.q) : false;
    }

    if(const xAOD::CaloCluster* c = in.caloCluster()) {
        out.clusE   = c->e()*MeV2GeV;
        out.clusEta = c->eta();
        out.clusPhi = c->phi();
    } else {
        all_available = false;
    }
    if(const xAOD::TrackParticle* t = in.trackParticle()){
        out.trackPt = t->pt()*MeV2GeV;
	out.d0      = t->d0();//AT:: wrt to PV ???

	const xAOD::Vertex* PV = getPV();
	double  primvertex_z = (PV) ? PV->z() : -999;
	out.z0 = t->z0() + t->vz() - primvertex_z;

	out.errD0         = Amg::error(t->definingParametersCovMatrix(),0);
	out.errZ0         = Amg::error(t->definingParametersCovMatrix(),1);
    } else {
        all_available = false;
    }
    // DG-2014-08-29 mc info not available yet
    // 
    // // Trigger flags
    // eleOut->trigFlags     = m_eleTrigFlags[ lepIn->idx() ];

 
    if(m_dbg && !all_available) cout<<"missing some electron variables"<<endl;
    m_susyNt.ele()->push_back(out);
}
//----------------------------------------------------------
// // match muon to muon truth, returns element object if success, else NULL
// D3PDReader::TruthMuonD3PDObjectElement* getMuonTruth(D3PDReader::MuonD3PDObject* muons, int muIdx, D3PDReader::TruthMuonD3PDObject* truthMuons)
// {
//     D3PDReader::TruthMuonD3PDObjectElement* result = NULL;
//     int bc = muons->truth_barcode()->at(muIdx);
//     if(bc==0){ // if barcode is zero then matching has already failed
//         return result;
//     }
//     // loop over truth muons, comparing barcode
//     for(int matchIdx=0; matchIdx < truthMuons->n(); matchIdx++){
//         if(bc == truthMuons->barcode()->at(matchIdx)){
//             result = & (*truthMuons)[matchIdx];
//             break;
//         }
//     }
//     return result;
// }
//----------------------------------------------------------
void SusyNtMaker::storeMuon(const xAOD::Muon &in)
{
    Susy::Muon out;
    double pt(in.pt()*MeV2GeV), eta(in.eta()), phi(in.phi()), m(in.m()*MeV2GeV);
    out.SetPtEtaPhiM(pt, eta, phi, m);
    out.pt  = pt;
    out.eta = eta;
    out.phi = phi;
    out.m   = m;
    out.q   = in.charge();
    out.isBaseline = in.auxdata< char >("baseline");
    out.isSignal   = in.auxdata< char >("signal");
    out.isCombined = in.muonType()==xAOD::Muon::Combined;
    out.isCosmic   = in.auxdata< char >("cosmic");
    out.isBadMuon  = m_susyObj[m_eleIDDefault]->IsBadMuon(in); // Uses default qoverpcut of 0.2

    bool all_available=true;

    // Isolation
    all_available &= in.isolation(out.etcone20, xAOD::Iso::etcone20); out.etcone20 *= MeV2GeV;
    all_available &= in.isolation(out.ptcone20, xAOD::Iso::ptcone20); out.ptcone20 *= MeV2GeV;
    all_available &= in.isolation(out.etcone30, xAOD::Iso::etcone30); out.etcone30 *= MeV2GeV;
    all_available &= in.isolation(out.ptcone30, xAOD::Iso::ptcone30); out.ptcone30 *= MeV2GeV;

    // ASM-2014-12-11 
    if(const xAOD::TrackParticle* t = in.primaryTrackParticle()){
      const xAOD::Vertex* PV = getPV();
      double  primvertex_z = (PV) ? PV->z() : 0.;
      out.d0             = t->d0();
      out.errD0          = Amg::error(t->definingParametersCovMatrix(),0); 
      out.z0             = t->z0() + t->vz() - primvertex_z;
      out.errZ0          = Amg::error(t->definingParametersCovMatrix(),1); 
    }
    // Inner Detector Track - if exists
    if(const xAOD::TrackParticle* idtrack = in.trackParticle( xAOD::Muon::InnerDetectorTrackParticle )){
      out.idTrackPt      = idtrack->pt()*MeV2GeV;  
      out.idTrackEta     = idtrack->eta();  
      out.idTrackPhi     = idtrack->phi(); 
      out.idTrackQ       = idtrack->qOverP() < 0 ? -1 : 1;
      out.id_qoverp      = idtrack->qOverP()*MeV2GeV;
      out.id_theta       = idtrack->theta();
      out.id_phi         = idtrack->phi();
    }
    // Muon Spectrometer Track - if exists
    if(const xAOD::TrackParticle* mstrack = in.trackParticle( xAOD::Muon::MuonSpectrometerTrackParticle )){
      out.msTrackPt      = mstrack->pt()*MeV2GeV;
      out.msTrackEta     = mstrack->eta();
      out.msTrackPhi     = mstrack->phi();
      out.msTrackQ       = mstrack->qOverP() < 0 ? -1 : 1;
      out.ms_qoverp      = mstrack->qOverP()*MeV2GeV;
      out.ms_theta       = mstrack->theta();
      out.ms_phi         = mstrack->phi();
    }

    // Truth Flags 
    if(m_isMC) {
      //AT 09/12/14 added Type/Origin
      const xAOD::TrackParticle* trackParticle = in.primaryTrackParticle();
      if(trackParticle){
	static SG::AuxElement::Accessor<int> acc_truthType("truthType");
	static SG::AuxElement::Accessor<int> acc_truthOrigin("truthOrigin");
	if (acc_truthType.isAvailable(*trackParticle)  ) out.mcType    = acc_truthType(*trackParticle);
	if (acc_truthOrigin.isAvailable(*trackParticle)) out.mcOrigin  = acc_truthOrigin(*trackParticle);

	const xAOD::TruthParticle* truthMu = xAOD::EgammaHelpers::getTruthParticle(trackParticle);
	out.matched2TruthLepton = truthMu ? true : false;
	int matchedPdgId = truthMu ? truthMu->pdgId() : -999;
	out.truthType  = isFakeLepton(out.mcOrigin, out.mcType, matchedPdgId); 

      }
      //// Old method tried to loop over all truth particles and do the matching by hand if above two were zero.
      //// ASM-2014-11-02, as as "AT 2014-10-29: Do not work... Need to know about all the truth particles in the event."
      //out.matched2TruthLepton = m_recoTruthMatch.Matched2TruthLepton(out); 
    }

    // Trigger Flags 
    // ASM-2014-11-02 :: Trigger information in DC14 samples are problematic
    // muOut->trigFlags      = m_muoTrigFlags[ lepIn->idx() ];

    // Scale Factors
    // ASM-2014-11-02 :: How to get the uncertatinty?
    if(m_isMC) {
      float value = 0.;
      CP::CorrectionCode result = m_muonEfficiencySFTool->getEfficiencyScaleFactor( in, value );

      if( result == CP::CorrectionCode::OutOfValidityRange ) {
        cout << "ASM :: getEfficiencyScaleFactor out of validity range " << endl;
      }
      else {
        out.effSF    = value; 
        out.errEffSF = 0.;  // ASM-2014-11-02 0. for the time being
      }
    }
    else {
      out.effSF    = 1.; 
      out.errEffSF = 0.; 
    }

    // ASM-2014-11-02 :: Store to be true at the moment
    all_available =  false;
    if(m_dbg && !all_available) cout<<"missing some muon variables"<<endl;
    m_susyNt.muo()->push_back(out);
}
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::storeJet(const xAOD::Jet &in)
{
    Susy::Jet out;
    double pt(in.pt()*MeV2GeV), eta(in.eta()), phi(in.phi()), m(in.m()*MeV2GeV);
    out.SetPtEtaPhiM(pt, eta, phi, m);
    out.pt  = pt;
    out.eta = eta;
    out.phi = phi;
    out.m   = m;
    bool all_available=true;

    // JVF 
    // ASM-2014-11-04 :: Remember JVT is gonna replace JVF in Run-II but not yet available
    vector<float> jetJVF;
    in.getAttribute(xAOD::JetAttribute::JVF,jetJVF); // JVF returns a vector that holds jvf per vertex
    const xAOD::Vertex* PV = getPV();                // Need to know the PV
    out.jvf = (PV) ? jetJVF.at(PV->index()) : 0.;    // Upon discussion w/ TJ (2014-12-11)   

    // Truth Label/Matching 
    if (m_isMC) in.getAttribute(xAOD::JetAttribute::JetLabel, out.truthLabel); 
    // jetOut->matchTruth    = m_isMC? matchTruthJet(jetIdx) : false;

    // B-tagging 
    out.mv1           = (in.btagging())->MV1_discriminant();
    out.sv1plusip3d   = (in.btagging())->SV1plusIP3D_discriminant();
    // Most of these are not available in DC14 samples, some obselete (ASM)
  // jetOut->sv0           = element->flavor_weight_SV0();
  // jetOut->combNN        = element->flavor_weight_JetFitterCOMBNN();
  // jetOut->jfit_mass     = element->flavor_component_jfit_mass();
  // jetOut->sv0p_mass     = element->flavor_component_sv0p_mass();
  // jetOut->svp_mass      = element->flavor_component_svp_mass();

    // Misc
    out.detEta = (in.jetP4(xAOD::JetConstitScaleMomentum)).eta();
    in.getAttribute(xAOD::JetAttribute::EMFrac,out.emfrac);
    in.getAttribute(xAOD::JetAttribute::BchCorrJet,out.bch_corr_jet);
    in.getAttribute(xAOD::JetAttribute::BchCorrCell,out.bch_corr_cell);

    // isBadJet 
    out.isBadVeryLoose = false; // DG tmp-2014-11-02 in.isAvailable("bad") ? in.auxdata<char>("bad") : false;

    // Hot Tile
    float fracSamplingMax, samplingMax;
    in.getAttribute(xAOD::JetAttribute::SamplingMax, samplingMax);
    in.getAttribute(xAOD::JetAttribute::FracSamplingMax, fracSamplingMax);
    const xAOD::EventInfo* eventinfo = XaodAnalysis::xaodEventInfo();
    out.isHotTile = m_susyObj[m_eleIDDefault]->isHotTile(eventinfo->runNumber(),
                                                         fracSamplingMax,
                                                         samplingMax,
                                                         eta, phi); 

  // // BCH cleaning flags - ASM-2014-11-04 :: Obsolete???
  // uint bchRun = m_isMC? m_mcRun : m_event.eventinfo.RunNumber();
  // uint bchLB = m_isMC? m_mcLB : m_event.eventinfo.lbn();
  // #define BCH_ARGS bchRun, bchLB, jetOut->detEta, jetOut->phi, jetOut->bch_corr_cell, jetOut->emfrac, jetOut->pt*1000.
  // jetOut->isBadMediumBCH = !m_susyObj[m_eleIDDefault]->passBCHCleaningMedium(BCH_ARGS, 0);
  // jetOut->isBadMediumBCH_up = !m_susyObj[m_eleIDDefault]->passBCHCleaningMedium(BCH_ARGS, 1);
  // jetOut->isBadMediumBCH_dn = !m_susyObj[m_eleIDDefault]->passBCHCleaningMedium(BCH_ARGS, -1);
  // jetOut->isBadTightBCH = !m_susyObj[m_eleIDDefault]->passBCHCleaningTight(BCH_ARGS);
  // #undef BCH_ARGS

  // // Save the met weights for the jets
  // // by checking status word similar to
  // // what is done in met utility
  // const D3PDReader::MissingETCompositionD3PDObjectElement &jetMetEgamma10NoTau = m_event.jet_AntiKt4LCTopo_MET_Egamma10NoTau[jetIdx];
  // // 0th element is what we care about
  // int sWord = jetMetEgamma10NoTau.statusWord().at(0);
  // bool passSWord = (MissingETTags::DEFAULT == sWord);       // Note assuming default met..
  // jetOut->met_wpx = 0; passSWord ? jetMetEgamma10NoTau.wpx().at(0) : 0;
  // jetOut->met_wpy = 0; passSWord ? jetMetEgamma10NoTau.wpy().at(0) : 0;
    if(m_dbg && !all_available) cout<<"missing some jet variables"<<endl;
    m_susyNt.jet()->push_back(out);
}
//----------------------------------------------------------
void SusyNtMaker::storePhoton(const xAOD::Photon &in)
{
    Susy::Photon out;
    double pt(in.pt()*MeV2GeV), eta(in.eta()), phi(in.phi()), m(in.m()*MeV2GeV);
    out.SetPtEtaPhiM(pt, eta, phi, m);
    out.pt  = pt;
    out.eta = eta;
    out.phi = phi;
    out.m   = m;
    out.isConv = xAOD::EgammaHelpers::isConvertedPhoton(&in);
    bool all_available=true;

    all_available &= in.passSelection(out.tight,"Tight");

    if(const xAOD::CaloCluster* c = in.caloCluster()) {
      out.clusE   = c->e()*MeV2GeV;
      out.clusEta = c->eta();
      out.clusPhi = c->phi();
    } else {
      all_available = false;
    }
    out.OQ = in.isGoodOQ(xAOD::EgammaParameters::BADCLUSPHOTON);
    in.isolationValue(out.topoEtcone40,xAOD::Iso::topoetcone40);

    if(m_dbg) cout << "AT: storePhoton: " << out.pt << " " << out.tight << " " << out.isConv << endl;
    // // Miscellaneous
    // phoOut->idx    = phIdx;
    // if(m_dbg>=5) cout << "fillPhotonVar" << endl;
    if(m_dbg && !all_available) cout<<"missing some photon variables"<<endl;
    m_susyNt.pho()->push_back(out);
}
//----------------------------------------------------------
void SusyNtMaker::storeTau(const xAOD::TauJet &tau)
{
    Susy::Tau out;
    double pt(tau.pt()*MeV2GeV), eta(tau.eta()), phi(tau.phi()), m(tau.m()*MeV2GeV);
    out.SetPtEtaPhiM(pt, eta, phi, m);
    out.pt  = pt;
    out.eta = eta;
    out.phi = phi;
    out.m   = m;
    bool all_available=true;
    out.q = tau.charge();
    
    // tauOut->author                = element->author(); // suneet: there is no author flag anymore?
    out.nTrack = tau.nTracks();
    out.eleBDT = tau.discriminant(xAOD::TauJetParameters::BDTEleScore);
    out.jetBDT = tau.discriminant(xAOD::TauJetParameters::BDTJetScore);

    out.jetBDTSigLoose = tau.isTau(xAOD::TauJetParameters::JetBDTSigLoose);
    out.jetBDTSigMedium = tau.isTau(xAOD::TauJetParameters::JetBDTSigMedium);
    out.jetBDTSigTight = tau.isTau(xAOD::TauJetParameters::JetBDTSigTight);

    out.eleBDTLoose = tau.isTau(xAOD::TauJetParameters::EleBDTLoose);
    out.eleBDTMedium = tau.isTau(xAOD::TauJetParameters::EleBDTMedium);
    out.eleBDTTight = tau.isTau(xAOD::TauJetParameters::EleBDTTight);

    out.muonVeto = tau.isTau(xAOD::TauJetParameters::MuonVeto);

    // tauOut->trueTau               = m_isMC? element->trueTauAssoc_matched() : false;
    
    // tauOut->matched2TruthLepton   = m_isMC? m_recoTruthMatch.Matched2TruthLepton(*tauLV, true) : false;
    // tauOut->detailedTruthType     = m_isMC? m_recoTruthMatch.TauDetailedFakeType(*tauLV) : -1;
    // tauOut->truthType             = m_isMC? m_recoTruthMatch.TauFakeType(tauOut->detailedTruthType) : -1;
    
    // // ID efficiency scale factors
    // if(m_isMC){
    //   #define TAU_ARGS TauCorrUncert::BDTLOOSE, tauLV->Eta(), element->numTrack()
    //   //TauCorrections* tauSF       = m_susyObj[m_eleIDDefault]->GetTauCorrectionsProvider();
    //   TauCorrUncert::TauSF* tauSF = m_susyObj[m_eleIDDefault]->GetSFTool();
    //   //tauOut->looseEffSF        = tauSF->GetIDSF(TauCorrUncert::BDTLOOSE, tauLV->Eta(), element->numTrack());
    //   //tauOut->mediumEffSF       = tauSF->GetIDSF(TauCorrUncert::BDTMEDIUM, tauLV->Eta(), element->numTrack());
    //   //tauOut->tightEffSF        = tauSF->GetIDSF(TauCorrUncert::BDTTIGHT, tauLV->Eta(), element->numTrack());
    //   //tauOut->errLooseEffSF     = tauSF->GetIDSFUnc(TauCorrUncert::BDTLOOSE, tauLV->Eta(), element->numTrack());
    //   //tauOut->errMediumEffSF    = tauSF->GetIDSFUnc(TauCorrUncert::BDTMEDIUM, tauLV->Eta(), element->numTrack());
    //   //tauOut->errTightEffSF     = tauSF->GetIDSFUnc(TauCorrUncert::BDTTIGHT, tauLV->Eta(), element->numTrack());
    //   tauOut->looseEffSF          = tauSF->GetIDSF(TAU_ARGS);
    //   tauOut->mediumEffSF         = tauSF->GetIDSF(TAU_ARGS);
    //   tauOut->tightEffSF          = tauSF->GetIDSF(TAU_ARGS);
    //   tauOut->errLooseEffSF       = sqrt(pow(tauSF->GetIDStatUnc(TAU_ARGS), 2) + pow(tauSF->GetIDSysUnc(TAU_ARGS), 2));
    //   tauOut->errMediumEffSF      = sqrt(pow(tauSF->GetIDStatUnc(TAU_ARGS), 2) + pow(tauSF->GetIDSysUnc(TAU_ARGS), 2));
    //   tauOut->errTightEffSF       = sqrt(pow(tauSF->GetIDStatUnc(TAU_ARGS), 2) + pow(tauSF->GetIDSysUnc(TAU_ARGS), 2));
    //   #undef TAU_ARGS
    
    //   if(element->numTrack()==1){
    //     float eta = element->leadTrack_eta();
    //     tauOut->looseEVetoSF      = tauSF->GetEVetoSF(eta, TauCorrUncert::BDTLOOSE, TauCorrUncert::LOOSE, TauCorrUncert::MEDIUMPP);
    //     tauOut->mediumEVetoSF     = tauSF->GetEVetoSF(eta, TauCorrUncert::BDTMEDIUM, TauCorrUncert::MEDIUM, TauCorrUncert::MEDIUMPP);
    //     // Doesn't currently work. Not sure why. Maybe they don't provide SFs for this combo
    //     //tauOut->tightEVetoSF      = tauSF->GetEVetoSF(eta, TauCorrUncert::BDTTIGHT, TauCorrUncert::TIGHT, TauCorrUncert::MEDIUMPP);
    //     tauOut->errLooseEVetoSF   = tauSF->GetEVetoSFUnc(eta, TauCorrUncert::BDTLOOSE, TauCorrUncert::LOOSE, TauCorrUncert::MEDIUMPP, 1);
    //     tauOut->errMediumEVetoSF  = tauSF->GetEVetoSFUnc(eta, TauCorrUncert::BDTMEDIUM, TauCorrUncert::MEDIUM, TauCorrUncert::MEDIUMPP, 1);
    //     //tauOut->errTightEVetoSF   = tauSF->GetEVetoSFUnc(eta, TauCorrUncert::BDTTIGHT, TauCorrUncert::TIGHT, TauCorrUncert::MEDIUMPP, 1);
    //   }
    // }
    
    // tauOut->trigFlags             = m_tauTrigFlags[tauIdx];
    
    // tauOut->idx   = tauIdx;
    if(m_dbg && !all_available) cout<<"missing some tau variables"<<endl;
    m_susyNt.tau()->push_back(out);
}

/*--------------------------------------------------------------------------------*/
// Fill MET variables
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillMetVars(SusyNtSys sys)
{

  xAOD::MissingETContainer::const_iterator met_it = m_metContainer->find("Final");
  
  if (met_it == m_metContainer->end()) {
    cout<<"No RefFinal inside MET container"<<endl;
    return;
  }
  
  m_susyNt.met()->push_back( Susy::Met() );
  Susy::Met* metOut = & m_susyNt.met()->back();

  metOut->Et = (*met_it)->met()*MeV2GeV;// m_met.Et();
  metOut->phi = (*met_it)->phi()*MeV2GeV;// m_met.Phi();
  metOut->sys = sys;
  metOut->sumet = (*met_it)->sumet()*MeV2GeV;

  if(m_dbg) cout << " AT:fillMetVars " << metOut->Et << " " << metOut->phi << " " << metOut->lv().Pt() << endl;

  // RefEle
  xAOD::MissingETContainer::const_iterator met_find = m_metContainer->find("RefEle");
  if (met_find == m_metContainer->end()) {
      // cout << "No RefEle inside MET container" << endl;
  }
  else {
      metOut->refEle = (*met_find)->met()*MeV2GeV;
      metOut->refEle_etx = (*met_find)->mpx()*MeV2GeV;
      metOut->refEle_ety = (*met_find)->mpy()*MeV2GeV;
      metOut->refEle_sumet = (*met_find)->sumet()*MeV2GeV;
  }
  
  // RefGamma
  met_find = m_metContainer->find("RefGamma");
  if (met_find == m_metContainer->end()) {
      // cout << "No RefGamma inside MET container" << endl;
  }
  else {
      metOut->refGamma = (*met_find)->met()*MeV2GeV;
      metOut->refGamma_etx = (*met_find)->mpx()*MeV2GeV;
      metOut->refGamma_ety = (*met_find)->mpy()*MeV2GeV;
      metOut->refGamma_sumet = (*met_find)->sumet()*MeV2GeV;
  }
  
  // RefTau
  met_find = m_metContainer->find("RefTau");
  if (met_find == m_metContainer->end()) {
      // cout << "No RefTau inside MET container" << endl;
  }
  else {
      // cout << "Found RefTau inside MET container, not stored." << endl;
  }

  // Muons
  met_find = m_metContainer->find("Muons");
  if (met_find == m_metContainer->end()) {
      // cout << "No Muons inside MET container" << endl;
  }
  else {
      metOut->refMuo = (*met_find)->met()*MeV2GeV;
      metOut->refMuo_etx = (*met_find)->mpx()*MeV2GeV;
      metOut->refMuo_ety = (*met_find)->mpy()*MeV2GeV;
      metOut->refMuo_sumet = (*met_find)->sumet()*MeV2GeV;
  }

  // RefJet
  met_find = m_metContainer->find("RefJet");
  if (met_find == m_metContainer->end()) {
      // cout << "No RefJet inside MET container" << endl;
  }
  else {
      metOut->refJet = (*met_find)->met()*MeV2GeV;
      metOut->refJet_etx = (*met_find)->mpx()*MeV2GeV;
      metOut->refJet_ety = (*met_find)->mpy()*MeV2GeV;
      metOut->refJet_sumet = (*met_find)->sumet()*MeV2GeV;
  }

  // SoftClus
  met_find = m_metContainer->find("SoftClus");
  if (met_find == m_metContainer->end()) {
      // cout << "No SoftClus (softTerm) inside MET container" << endl;
  }
  else {
      metOut->softTerm = (*met_find)->met()*MeV2GeV;
      metOut->softTerm_etx = (*met_find)->mpx()*MeV2GeV;
      metOut->softTerm_ety = (*met_find)->mpy()*MeV2GeV;
      metOut->softTerm_sumet = (*met_find)->sumet()*MeV2GeV;
  }

  // cout << "Done looking for MET terms!" << endl;

#warning fillMetVars not implemented
  // if(m_dbg>=5) cout << "fillMetVars: sys " << sys << endl;

  // // Just fill the lv for now
  // double Et  = m_met.Et()/GeV;
  // double phi = m_met.Phi();

  // //double px = m_met.Px()/GeV;
  // //double py = m_met.Py()/GeV;
  // //double pz = m_met.Pz()/GeV;
  // //double E  = m_met.E()/GeV;

  // // Need to get the metUtility in order to
  // // get all the sumet terms.  In the future,
  // // we could use the metUtility to get all the
  // // comonents instead of the SUSYTools method
  // // computeMetComponent, but that is up to Steve,
  // // Lord of the Ntuples.
  // METUtility* metUtil = m_susyObj[m_eleIDDefault]->GetMETUtility();

  // m_susyNt.met()->push_back( Susy::Met() );
  // Susy::Met* metOut = & m_susyNt.met()->back();
  // metOut->Et    = Et;
  // metOut->phi   = phi;
  // metOut->sys   = sys;
  // metOut->sumet = metUtil->getMissingET(METUtil::RefFinal, METUtil::None).sumet()/GeV;

  // // MET comp terms
  // // Need to save these for the MET systematics as well.
  // // Use the sys enum to determine which argument to pass to SUSYTools
  // METUtil::Systematics metSys = METUtil::None;

  // // I guess these are the only ones we need to specify, the ones specified in SUSYTools...
  // // All the rest should be automatic (I think), e.g. JES
  // if(sys == NtSys_SCALEST_UP) metSys = METUtil::ScaleSoftTermsUp;
  // else if(sys == NtSys_SCALEST_DN) metSys = METUtil::ScaleSoftTermsDown;
  // else if(sys == NtSys_RESOST) metSys = METUtil::ResoSoftTermsUp;

  // // Save the MET terms
  // TVector2 refEleV   = m_susyObj[m_eleIDDefault]->computeMETComponent(METUtil::RefEle, metSys);
  // TVector2 refMuoV   = m_susyObj[m_eleIDDefault]->computeMETComponent(METUtil::MuonTotal, metSys);
  // TVector2 refJetV   = m_susyObj[m_eleIDDefault]->computeMETComponent(METUtil::RefJet, metSys);
  // TVector2 refGammaV = m_susyObj[m_eleIDDefault]->computeMETComponent(METUtil::RefGamma, metSys);
  // //TVector2 softJetV  = m_susyObj[m_eleIDDefault]->computeMETComponent(METUtil::SoftJets, metSys);
  // //TVector2 refCellV  = m_susyObj[m_eleIDDefault]->computeMETComponent(METUtil::CellOutEflow, metSys);
  // TVector2 softTermV = m_susyObj[m_eleIDDefault]->computeMETComponent(METUtil::SoftTerms, metSys);
  // //float sumet = m_susyObj[m_eleIDDefault]->_metUtility->getMissingET(METUtil::SoftTerms).sumet();

  //migrated// metOut->refEle     = refEleV.Mod()/GeV;
  //migrated// metOut->refEle_etx = refEleV.Px()/GeV;
  //migrated// metOut->refEle_ety = refEleV.Py()/GeV;
  //migrated// metOut->refEle_sumet = metUtil->getMissingET(METUtil::RefEle, metSys).sumet()/GeV;

  //migrated// metOut->refMuo     = refMuoV.Mod()/GeV;
  //migrated// metOut->refMuo_etx = refMuoV.Px()/GeV;
  //migrated// metOut->refMuo_ety = refMuoV.Py()/GeV;
  //migrated// metOut->refMuo_sumet = metUtil->getMissingET(METUtil::MuonTotal, metSys).sumet()/GeV;

  //migrated// metOut->refJet     = refJetV.Mod()/GeV;
  //migrated// metOut->refJet_etx = refJetV.Px()/GeV;
  //migrated// metOut->refJet_ety = refJetV.Py()/GeV;
  //migrated// metOut->refJet_sumet = metUtil->getMissingET(METUtil::RefJet, metSys).sumet()/GeV;

  //migrated// metOut->refGamma     = refGammaV.Mod()/GeV;
  //migrated// metOut->refGamma_etx = refGammaV.Px()/GeV;
  //migrated// metOut->refGamma_ety = refGammaV.Py()/GeV;
  //migrated// metOut->refGamma_sumet = metUtil->getMissingET(METUtil::RefGamma, metSys).sumet()/GeV;

  // //metOut->softJet     = softJetV.Mod()/GeV;
  // //metOut->softJet_etx = softJetV.Px()/GeV;
  // //metOut->softJet_ety = softJetV.Py()/GeV;

  // //metOut->refCell     = refCellV.Mod()/GeV;
  // //metOut->refCell_etx = refCellV.Px()/GeV;
  // //metOut->refCell_ety = refCellV.Py()/GeV;

  //migrated// metOut->softTerm     = softTermV.Mod()/GeV;
  //migrated// metOut->softTerm_etx = softTermV.Px()/GeV;
  //migrated// metOut->softTerm_ety = softTermV.Py()/GeV;
  //migrated// metOut->softTerm_sumet = metUtil->getMissingET(METUtil::SoftTerms, metSys).sumet()/GeV;

  // //metOut->refEle        = m_susyObj[m_eleIDDefault]->computeMETComponent(METUtil::RefEle, metSys).Mod()/GeV;
  // //metOut->refMuo        = m_susyObj[m_eleIDDefault]->computeMETComponent(METUtil::MuonTotal, metSys).Mod()/GeV;
  // //metOut->refJet        = m_susyObj[m_eleIDDefault]->computeMETComponent(METUtil::RefJet, metSys).Mod()/GeV;
  // //metOut->refGamma      = m_susyObj[m_eleIDDefault]->computeMETComponent(METUtil::RefGamma, metSys).Mod()/GeV;
  // //metOut->softJet       = m_susyObj[m_eleIDDefault]->computeMETComponent(METUtil::SoftJets, metSys).Mod()/GeV;
  // //metOut->refCell       = m_susyObj[m_eleIDDefault]->computeMETComponent(METUtil::CellOutEflow, metSys).Mod()/GeV;
}
//----------------------------------------------------------
void SusyNtMaker::storeTruthParticle(const xAOD::TruthParticle &in)
{
    Susy::TruthParticle out;
    double pt(in.pt()*MeV2GeV), eta(in.eta()), phi(in.phi()), m(in.m()*MeV2GeV);
    out.SetPtEtaPhiM(pt, eta, phi, m);
    out.pt  = pt;
    out.eta = eta;
    out.phi = phi;
    out.m   = m;
    bool all_available=true;
    // out.charge = in.charge(); // DG 2014-08-29 discards const ??
    out.pdgId = in.pdgId();
    out.status = in.status();
  //   tprOut->motherPdgId = smc::determineParentPdg(m_event.mc.pdgId(),
  //                                                 m_event.mc.parent_index(),
  //                                                 truParIdx);
}
/*--------------------------------------------------------------------------------*/
// Fill Truth Jet variables
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillTruthJetVars()
{
#warning fillTruthJetVars not implemented
  // if(m_dbg>=5) cout << "fillTruthJetVars" << endl;

  // for(uint iTruJet=0; iTruJet<m_truJets.size(); iTruJet++){
  //   int truJetIdx = m_truJets[iTruJet];

  //   m_susyNt.tjt()->push_back( Susy::TruthJet() );
  //   Susy::TruthJet* truJetOut = & m_susyNt.tjt()->back();
  //   const D3PDReader::JetD3PDObjectElement* element = & m_event.AntiKt4Truth[truJetIdx];

  //   // Set TLV
  //   float pt  = element->pt() / GeV;
  //   float eta = element->eta();
  //   float phi = element->phi();
  //   float m   = element->m()  / GeV;

  //   truJetOut->SetPtEtaPhiM(pt, eta, phi, m);
  //   truJetOut->pt     = pt;
  //   truJetOut->eta    = eta;
  //   truJetOut->phi    = phi;
  //   truJetOut->m      = m;

  //   truJetOut->flavor = element->flavor_truth_label();
  // }
}
/*--------------------------------------------------------------------------------*/
// Fill Truth Met variables
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillTruthMetVars()
{
#warning fillTruthMetVars not implemented
  // if(m_dbg>=5) cout << "fillTruthMetVars" << endl;

  // // Just fill the lv for now
  // double Et  = m_truMet.Et()/GeV;
  // double phi = m_truMet.Phi();

  // m_susyNt.tmt()->push_back( Susy::TruthMet() );
  // Susy::TruthMet* truMetOut = & m_susyNt.tmt()->back();
  // truMetOut->Et  = Et;
  // truMetOut->phi = phi;
}



/*--------------------------------------------------------------------------------*/
// Handle Systematic
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::doSystematic()
{
  if(m_dbg>=5) cout<< "doSystematic " << sysList.size() << endl;
  std::vector<CP::SystematicSet>::iterator sysListItr;
  for (sysListItr = sysList.begin(); sysListItr != sysList.end(); ++sysListItr){
    if((*sysListItr).name()=="") continue; // skip Nominal
    NtSys::SusyNtSys sys = NtSys::CPsys2sys((*sysListItr).name());
    if( sys == NtSys::SYSUNKNOWN ) continue;
    if(!NtSys::isObjectSystematic(sys)) continue;
    
    cout << "Found syst in global registry: " << (*sysListItr).name() << endl;
    cout << " Our systematic " << NtSys::SusyNtSysNames[sys] << endl;
    
    if ( m_susyObj[m_eleIDDefault]->applySystematicVariation(*sysListItr) != CP::SystematicCode::Ok){
      cout << "SusyNtMaker::doSystematic - cannot configure SUSYTools for " << (*sysListItr).name() << endl;
      continue;
    }

    
    /*
      Do our stuff here
    */



     m_susyObj[m_eleIDDefault]->resetSystematics();
  }

  // Loop over the systematics: start at 1, nominal already saved
  /*
  for(int i = 1; i < NtSys_N; i++){
      SusyNtSys sys = static_cast<SusyNtSys>(i);
    if(m_dbg>=5) cout << "Doing sys " << SusyNtSystNames[sys] << endl;
    clearOutputObjects();
    selectObjects(sys);
    buildMet(sys);
    assignEventCleaningFlags();
    assignObjectCleaningFlags();
    if     (isElecSys(sys)) saveElectronSF(sys); // Lepton Specific sys
    else if(isMuonSys(sys)) saveMuonSF(sys);
    else if(isJetSys(sys))  saveJetSF(sys);
    else if(isTauSys(sys))  saveTauSF(sys);
    fillMetVars(sys);
    m_susyNt.evt()->cutFlags[sys] = m_cutFlags;
  }
  */
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::saveElectronSF(SusyNtSys sys)
{
#warning saveElectronSF not implemented
  // // Loop over preselected leptons and fill the systematic shifts
  // for(uint iLep=0; iLep < m_preLeptons.size(); iLep++){
  //   const LeptonInfo* lep = & m_preLeptons[iLep];
  //   if(!lep->isElectron()) continue;

  //   // Systematic shifted energy
  //   float E_sys = lep->lv()->E() / GeV;

  //   // Try to find this electron in the list of SusyNt electrons
  //   Susy::Electron* eleOut = 0;
  //   for(uint iEl=0; iEl<m_susyNt.ele()->size(); iEl++){
  //     Susy::Electron* ele = & m_susyNt.ele()->at(iEl);
  //     if(ele->idx == lep->idx()){
  //       eleOut = ele;
  //       break;
  //     }
  //   }

  //   // If electron not found, then we need to add it
  //   if(eleOut == 0){
  //     addMissingElectron(lep, sys);
  //     eleOut = & m_susyNt.ele()->back();
  //   }

  //   // Calculate systematic scale factor
  //   float sf = E_sys / eleOut->E();
  //   if(sys == NtSys_EES_Z_UP)        eleOut->ees_z_up = sf;
  //   else if(sys == NtSys_EES_Z_DN)   eleOut->ees_z_dn = sf;
  //   else if(sys == NtSys_EES_MAT_UP) eleOut->ees_mat_up = sf;
  //   else if(sys == NtSys_EES_MAT_DN) eleOut->ees_mat_dn = sf;
  //   else if(sys == NtSys_EES_PS_UP)  eleOut->ees_ps_up = sf;
  //   else if(sys == NtSys_EES_PS_DN)  eleOut->ees_ps_dn = sf;
  //   else if(sys == NtSys_EES_LOW_UP) eleOut->ees_low_up = sf;
  //   else if(sys == NtSys_EES_LOW_DN) eleOut->ees_low_dn = sf;
  //   else if(sys == NtSys_EER_UP)     eleOut->eer_up = sf;
  //   else if(sys == NtSys_EER_DN)     eleOut->eer_dn = sf;
  // }
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::saveMuonSF(SusyNtSys sys)
{
#warning saveMuonSF not implemented
  // // Loop over preselected leptons and fill the systematic shifts
  // for(uint iLep=0; iLep < m_preLeptons.size(); iLep++){
  //   const LeptonInfo* lep = & m_preLeptons[iLep];
  //   if(lep->isElectron()) continue;

  //   // Systematic shifted energy
  //   float E_sys = lep->lv()->E() / GeV;

  //   // Try to find this muon in the list of SusyNt muons
  //   Susy::Muon* muOut = 0;
  //   for(uint iMu=0; iMu<m_susyNt.muo()->size(); iMu++){
  //     Susy::Muon* mu = & m_susyNt.muo()->at(iMu);
  //     if(mu->idx == lep->idx()){
  //       muOut = mu;
  //       break;
  //     }
  //   }

  //   // If muon not found, then we need to add it
  //   if(muOut == 0){
  //     addMissingMuon(lep, sys);
  //     muOut = & m_susyNt.muo()->back();
  //   }

  //   // Calculate systematic scale factor
  //   float sf = E_sys / muOut->E();
  //   if(sys == NtSys_MS_UP)      muOut->ms_up = sf;
  //   else if(sys == NtSys_MS_DN) muOut->ms_dn = sf;
  //   else if(sys == NtSys_ID_UP) muOut->id_up = sf;
  //   else if(sys == NtSys_ID_DN) muOut->id_dn = sf;
  //  } // end loop over leptons
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::saveJetSF(SusyNtSys sys)
{
#warning saveJetSF not implemented
  // // Loop over selected jets and fill the systematic shifts
  // for(uint iJet=0; iJet<m_preJets.size(); iJet++){
  //   uint jetIdx = m_preJets[iJet];

  //   // Systematic shifted energy
  //   float E_sys = m_susyObj[m_eleIDDefault]->GetJetTLV(jetIdx).E() / GeV;

  //   // Try to find this jet in the list of SusyNt jets
  //   Susy::Jet* jetOut = 0;
  //   for(uint iJ=0; iJ<m_susyNt.jet()->size(); ++iJ){
  //     Susy::Jet* jet = & m_susyNt.jet()->at(iJ);
  //     if(jet->idx == jetIdx){
  //       jetOut = jet;
  //       break;
  //     }
  //   }

  //   // If jet not found, then we need to add it
  //   if(jetOut == 0){
  //     addMissingJet(jetIdx, sys);
  //     jetOut = & m_susyNt.jet()->back();
  //   }

  //   // Calculate systematic scale factor
  //   float sf = E_sys / jetOut->E();
  //   if(sys == NtSys_JES_UP)      jetOut->jes_up = sf;
  //   else if(sys == NtSys_JES_DN) jetOut->jes_dn = sf;
  //   else if(sys == NtSys_JER)    jetOut->jer = sf;
  // } // end loop over jets in pre-jets
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::saveTauSF(SusyNtSys sys)
{
#warning saveTauSF not implemented
  // // Loop over preselected taus and fill systematic shifts
  // for(uint iTau=0; iTau<m_preTaus.size(); iTau++){
  //   uint tauIdx = m_preTaus[iTau];

  //   // Get the systematic shifted E, used to calculate a shift factor
  //   float E_sys = m_susyObj[m_eleIDDefault]->GetTauTLV(tauIdx).E() / GeV;

  //   // Try to find this tau in the list of SusyNt taus
  //   Susy::Tau* tauOut = 0;
  //   for(uint iT=0; iT<m_susyNt.tau()->size(); iT++){
  //     Susy::Tau* tau = & m_susyNt.tau()->at(iT);
  //     if(tau->idx == tauIdx){
  //       tauOut = tau;
  //       break;
  //     }
  //   }
  //   // If tau not found, then it was not nominally pre-selected and must be added now
  //   if(tauOut == 0){
  //     addMissingTau(tauIdx, sys);
  //     tauOut = & m_susyNt.tau()->back();
  //   }

  //   // Calculate systematic scale factor
  //   float sf = E_sys / tauOut->E();
  //   if(sys == NtSys_TES_UP) tauOut->tes_up = sf;
  //   if(sys == NtSys_TES_DN) tauOut->tes_dn = sf;
  // }
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::addMissingElectron(const LeptonInfo* lep, SusyNtSys sys)
{
  // // This electron did not pass nominal cuts, and therefore
  // // needs to be added, but with the correct TLV

  // // Reset the Nominal TLV
  // // NOTE: this overwrites the TLV in SUSYObjDef with the nominal variables,
  // // regardless of our current systematic.
  // const D3PDReader::ElectronD3PDObjectElement* element = lep->getElectronElement();
  // m_susyObj[m_eleIDDefault]->SetElecTLV(lep->idx(), element->eta(), element->phi(), element->cl_eta(), element->cl_phi(), element->cl_E(),
  //                      element->tracketa(), element->trackphi(), element->nPixHits(), element->nSCTHits(), SystErr::NONE);
  // // Now push it back onto to susyNt
  // fillElectronVars(lep);
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::addMissingMuon(const LeptonInfo* lep, SusyNtSys sys)
{
  // // This muon did not pass nominal cuts, and therefore
  // // needs to be added, but with the correct TLV
  // // Reset the Nominal TLV
  // // NOTE: this overwrites the TLV in SUSYObjDef with the nominal variables,
  // // regardless of our current systematic.
  // const D3PDReader::MuonD3PDObjectElement* element = lep->getMuonElement();
  // m_susyObj[m_eleIDDefault]->SetMuonTLV(lep->idx(), element->pt(), element->eta(), element->phi(),
  //                      element->me_qoverp_exPV(), element->id_qoverp_exPV(), element->me_theta_exPV(),
  //                      element->id_theta_exPV(), element->charge(), element->isCombinedMuon(),
  //                      element->isSegmentTaggedMuon(), SystErr::NONE);
  // fillMuonVars(lep);
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::addMissingJet(int index, SusyNtSys sys)
{
  // // Get the systematic shifted E, used to calculate a shift factor
  // //TLorentzVector tlv_sys = m_susyObj[m_eleIDDefault]->GetJetTLV(index);
  // //float E_sys = m_susyObj[m_eleIDDefault]->GetJetTLV(index).E();

  // // Reset the Nominal TLV
  // // NOTE: this overwrites the TLV in SUSYObjDef with the nominal variables,
  // // regardless of our current systematic.
  // const D3PDReader::JetD3PDObjectElement* jet = &m_event.jet_AntiKt4LCTopo[index];
  // m_susyObj[m_eleIDDefault]->FillJet(index, jet->pt(), jet->eta(), jet->phi(), jet->E(),
  //                   jet->constscale_eta(), jet->constscale_phi(), jet->constscale_E(), jet->constscale_m(),
  //                   jet->ActiveAreaPx(), jet->ActiveAreaPy(), jet->ActiveAreaPz(), jet->ActiveAreaE(),
  //                   m_event.Eventshape.rhoKt4LC(),
  //                   m_event.eventinfo.averageIntPerXing(),
  //                   m_event.vxp.nTracks());
  // fillJetVar(index);
  // // Set SF This should only be done in saveJetSF
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::addMissingTau(int index, SusyNtSys sys)
{
  // // This tau did not pass nominal cuts, and therefore
  // // needs to be added, but with the correct TLV.
  // // Get the systematic shifted E, used to calculate a shift factor
  // //TLorentzVector tlv_sys = m_susyObj[m_eleIDDefault]->GetTauTLV(index);
  // //float E_sys = m_susyObj[m_eleIDDefault]->GetTauTLV(index).E();
  // // Grab the d3pd variables
  // const D3PDReader::TauD3PDObjectElement* element = & m_event.tau[index];
  // // Reset the Nominal TLV
  // // NOTE: this overwrites the TLV in SUSYObjDef with the nominal variables,
  // // regardless of our current systematic.
  // m_susyObj[m_eleIDDefault]->SetTauTLV(index, element->pt(), element->eta(), element->phi(), element->Et(), element->numTrack(),
  //                     element->leadTrack_eta(), SUSYTau::TauMedium, SystErr::NONE, true);
  // // Fill the tau vars for this guy
  // fillTauVar(index);
}
/*--------------------------------------------------------------------------------*/
//bool SusyNtMaker::isBuggyWwSherpaSample(const int &dsid)
//{
//  return (dsid==126892 || dsid==157817 || dsid==157818 || dsid==157819);
//}
/*--------------------------------------------------------------------------------*/
//bool SusyNtMaker::hasRadiativeBquark(const vint_t *pdg, const vint_t *status)
//{
//  if(!pdg || !status || pdg->size()!=status->size()) return false;
//  const vint_t &p = *pdg;
//  const vint_t &s = *status;
//  const int pdgB(5), statRad(3);
//  for(size_t i=0; i<p.size(); ++i) if(abs(p[i])==pdgB && s[i]==statRad) return true;
//  return false;
//}
//----------------------------------------------------------
SusyNtMaker& SusyNtMaker::initializeOuputTree()
{
    m_outTreeFile = new TFile("susyNt.root", "recreate");
    m_outTree = new TTree("susyNt", "susyNt");
    m_outTree->SetAutoSave(10000000); // DG-2014-08-15 magic numbers, ask Steve
    m_outTree->SetMaxTreeSize(3000000000u);
    m_susyNt.SetActive();
    m_susyNt.WriteTo(m_outTree);
    return *this;
}
//----------------------------------------------------------
SusyNtMaker& SusyNtMaker::initializeCutflowHistograms()
{
  h_rawCutFlow = makeCutFlow("rawCutFlow", "rawCutFlow;Cuts;Events");
  h_genCutFlow = makeCutFlow("genCutFlow", "genCutFlow;Cuts;Events");
  return *this;
}
//----------------------------------------------------------
bool SusyNtMaker::guessWhetherIsWhSample(const TString &samplename)
{
    return (samplename.Contains("simplifiedModel_wA_noslep_WH") ||
            samplename.Contains("Herwigpp_sM_wA_noslep_notauhad_WH"));
}
//----------------------------------------------------------
SusyNtMaker& SusyNtMaker::saveOutputTree()
{
    m_outTreeFile = m_outTree->GetCurrentFile();
    m_outTreeFile->Write(0, TObject::kOverwrite);
    cout<<"susyNt tree saved to "<<m_outTreeFile->GetName()<<endl;
    m_outTreeFile->Close();
    return *this;
}
//----------------------------------------------------------
std::string SusyNtMaker::timerSummary() /*const*/ // TStopwatch::<*>Time is not const
{
  double realTime = m_timer.RealTime();
  double cpuTime  = m_timer.CpuTime();
  int hours = int(realTime / 3600);
  realTime -= hours * 3600;
  int min   = int(realTime / 60);
  realTime -= min * 60;
  int sec   = int(realTime);
  int nEventInput = m_cutstageCounters.front();
  int nEventOutput = m_outTree ? m_outTree->GetEntries() : -1;
  float speed = nEventInput/m_timer.RealTime()/1000;
  TString line1; line1.Form("Real %d:%02d:%02d, CPU %.3f", hours, min, sec, cpuTime);
  TString line2; line2.Form("[kHz]: %2.3f",speed);
  ostringstream oss;
  oss<<"---------------------------------------------------\n"
     <<" Number of events processed: "<<nEventInput<<endl
     <<" Number of events saved:     "<<nEventOutput<<endl
     <<"\t Analysis time: "<<line1<<endl
     <<"\t Analysis speed "<<line2<<endl
     <<"---------------------------------------------------"<<endl
     <<endl;
  return oss.str();
}
//----------------------------------------------------------
std::string SusyNtMaker::counterSummary() const
{
  ostringstream oss;
  oss<<"Object counter"<<endl
     <<"  BaseEle   "<<n_base_ele   <<endl
     <<"  BaseMuo   "<<n_base_muo   <<endl
     <<"  BaseTau   "<<n_base_tau   <<endl
     <<"  BaseJet   "<<n_base_jet   <<endl
     <<"  SigEle    "<<n_sig_ele    <<endl
     <<"  SigMuo    "<<n_sig_muo    <<endl
     <<"  SigTau    "<<n_sig_tau    <<endl
     <<"  SigJet    "<<n_sig_jet    <<endl
     <<endl;

  oss<<"Event counter"<<endl;
  vector<string> labels = SusyNtMaker::cutflowLabels();
  struct shorter { bool operator()(const string &a, const string &b) { return a.size() < b.size(); } };
  size_t max_label_length = max_element(labels.begin(), labels.end(), shorter())->size();

  for(size_t i=0; i<m_cutstageCounters.size(); ++i)
      oss<<"  "<<setw(max_label_length+2)<<std::left<<labels[i]<<m_cutstageCounters[i]<<endl;
  oss<<endl;
  return oss.str();
}
//----------------------------------------------------------
struct FillCutFlow { ///< local function object to fill the cutflow histograms
    TH1 *raw, *gen, *perProcess; ///< ptr to histos with counters
    int iCut; ///< index of the sequential cut (must match bin labels, see SusyNtMaker::makeCutFlow())
    bool passAll; ///< whether we've survived all cuts so far
    bool includeThisCut_; ///< whether this cut should be used when computing passAll
    vector< size_t > *counters;
    FillCutFlow(TH1 *r, TH1* g, TH1* p, vector< size_t > *cs) :
        raw(r), gen(g), perProcess(p), iCut(0), passAll(true), includeThisCut_(true), counters(cs) {}
    void operator()(bool thisEventDoesPassThisCut, float weight) {
        if(thisEventDoesPassThisCut && passAll) {
            if(raw       ) raw       ->Fill(iCut);
            if(gen       ) gen       ->Fill(iCut, weight);
            if(perProcess) perProcess->Fill(iCut, weight);
            counters->at(iCut) += 1;
        } else {
            if(includeThisCut_) passAll = false;
        }
        iCut++;
    }
    FillCutFlow& includeThisCut(bool v) { includeThisCut_ = v; return *this; }
};
//----------------------------------------------------------
bool SusyNtMaker::passEventlevelSelection()
{
    TH1F* h_procCutFlow = getProcCutFlow(m_susyFinalState);
    float w = m_susyNt.evt()->w;

    FillCutFlow fillCutFlow(h_rawCutFlow, h_genCutFlow, h_procCutFlow, &m_cutstageCounters);

    assignEventCleaningFlags();
    bool keep_all_events(!m_filter);
    bool pass_susyprop(!m_hasSusyProp);
    bool pass_grl(m_cutFlags & ECut_GRL), pass_lar(m_cutFlags & ECut_LarErr), pass_tile(m_cutFlags & ECut_TileErr);
    bool pass_ttc(m_cutFlags & ECut_TTC), pass_goodpv(m_cutFlags & ECut_GoodVtx), pass_tiletrip(m_cutFlags & ECut_TileTrip);
    bool pass_wwfix(true); //(!m_isMC || (m_susyObj[m_eleIDDefault]->Sherpa_WW_veto())); // DG-2014-08-16 sherpa ww bugfix probably obsolete

    fillCutFlow(true, w); // initial bin
    fillCutFlow.includeThisCut(false); // susyProp just counts (for normalization), doesn't drop
    fillCutFlow(pass_susyprop, w);
    fillCutFlow.includeThisCut(true);
    fillCutFlow(pass_grl, w);
    fillCutFlow(pass_lar, w);
    fillCutFlow(pass_tile, w);
    fillCutFlow(pass_ttc, w);
    fillCutFlow(pass_goodpv, w);
    fillCutFlow(pass_wwfix, w);
    if(m_dbg>=5 &&  !(keep_all_events || fillCutFlow.passAll) ) 
      cout << "SusyNtMaker fail passEventlevelSelection " 
	   << keep_all_events << " " << fillCutFlow.passAll <<  endl;
    return (keep_all_events || fillCutFlow.passAll);
}
//----------------------------------------------------------
bool SusyNtMaker::passObjectlevelSelection()
{
  SusyNtSys sys=NtSys::NOM;
    selectObjects(sys);
    // buildMet(sys);
    // assignObjectCleaningFlags();

    n_base_ele += m_baseElectrons.size();
    n_base_muo += m_baseMuons.size();
    n_base_tau += m_baseTaus.size();
    n_base_jet += m_baseJets.size();
    n_sig_ele += m_sigElectrons.size();
    n_sig_muo += m_sigMuons.size();
    n_sig_tau += m_sigTaus.size();
    n_sig_jet += m_sigJets.size();
    TH1F* h_procCutFlow = getProcCutFlow(m_susyFinalState);
    FillCutFlow fillCutFlow(h_rawCutFlow, h_genCutFlow, h_procCutFlow, &m_cutstageCounters);
    fillCutFlow.iCut = 8; // we've filled up to 'Buggy WWSherpa' in passEventlevelSelection
    bool pass_hotspot(true), pass_basdjet(true), pass_badmuon(true), pass_cosmic(true); // dummy values, todo
    bool pass_ge1l(1<=(m_sigElectrons.size()+m_sigMuons.size()));
    bool pass_ge2l(2<=(m_sigElectrons.size()+m_sigMuons.size()));
    bool pass_eq3l(3==(m_sigElectrons.size()+m_sigMuons.size()));
    float w = m_susyNt.evt()->w;
    fillCutFlow(pass_hotspot, w);
    fillCutFlow(pass_basdjet, w);
    fillCutFlow(pass_badmuon, w);
    fillCutFlow(pass_cosmic, w);
    fillCutFlow(pass_ge1l, w);
    fillCutFlow(pass_ge2l, w);
    fillCutFlow(pass_eq3l, w);
    bool has_at_least_one_lepton = pass_ge1l;
    if(m_dbg>=5 && !has_at_least_one_lepton)
      cout << "SusyNtMaker: fail passObjectlevelSelection " << endl;
    return has_at_least_one_lepton;
}
//----------------------------------------------------------
#undef GeV
