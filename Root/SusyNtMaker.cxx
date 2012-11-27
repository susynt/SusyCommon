#include "egammaAnalysisUtils/CaloIsoCorrection.h"
//#include "SUSYTools/MV1.h"

#include "MultiLep/MuonTools.h"
#include "MultiLep/ElectronTools.h"
#include "MultiLep/SusyGridCrossSectionTools.h"
#include "MultiLep/TruthTools.h"

#include "SusyCommon/SusyNtMaker.h"

using namespace std;

#define GeV 1000.

/*--------------------------------------------------------------------------------*/
// SusyNtMaker Constructor
/*--------------------------------------------------------------------------------*/
SusyNtMaker::SusyNtMaker() : m_fillNt(true)
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
  n_evt_ttcVeto=0;
  n_evt_larErr=0;
  n_evt_tileErr=0;
  n_evt_larHole=0;
  n_evt_hotSpot=0;
  n_evt_badJet=0;
  n_evt_goodVtx=0;
  n_evt_badMu=0;
  n_evt_cosmic=0;
  n_evt_1Lep=0;
  n_evt_2Lep=0;
  n_evt_3Lep=0;
  n_evt_saved=0;
}
/*--------------------------------------------------------------------------------*/
// Destructor
/*--------------------------------------------------------------------------------*/
SusyNtMaker::~SusyNtMaker()
{
}
/*--------------------------------------------------------------------------------*/
// The Begin() function is called at the start of the query.
// When running with PROOF Begin() is only called on the client.
// The tree argument is deprecated (on PROOF 0 is passed).
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::Begin(TTree* /*tree*/)
{
  SusyD3PDAna::Begin(0);
  if(m_dbg) cout << "SusyNtMaker::Begin" << endl;

  if(m_fillNt){

    // Open the output tree
    m_outTreeFile = new TFile("susyNt.root", "recreate");
    m_outTree = new TTree("susyNt", "susyNt");

    // Set autosave size (determines how often tree writes to disk)
    m_outTree->SetAutoSave(10000000);
    // Max tree size determines when a new file and tree are written
    m_outTree->SetMaxTreeSize(3000000000u);
    // Set all branches active for writing, for now.
    // Later, add switch for systematics
    m_susyNt.SetActive();
    m_susyNt.WriteTo(m_outTree);

  }

  m_isSusySample = m_sample.Contains("DGemt") || m_sample.Contains("DGstau") || m_sample.Contains("RPV");

  // create histogram for cutflow
  //h_cutFlow = new TH1F("cutFlow","Histogram storing cuts applied upstream", 5, -0.5, 3.5);
  //h_cutFlow->GetXaxis()->SetBinLabel(1, "total");
  //h_cutFlow->GetXaxis()->SetBinLabel(2, "GRL");
  //h_cutFlow->GetXaxis()->SetBinLabel(3, "LAr Error");
  //h_cutFlow->GetXaxis()->SetBinLabel(4, "TTC Veto");
  //h_cutFlow->GetXaxis()->SetBinLabel(5, "Good Vertex");

  // Raw event weights
  h_rawCutFlow = makeCutFlow("rawCutFlow", "rawCutFlow;Cuts;Events");
  // Generator event weights
  h_genCutFlow = makeCutFlow("genCutFlow", "genCutFlow;Cuts;Events");

  // Start the timer
  m_timer.Start();
}
/*--------------------------------------------------------------------------------*/
TH1F* SusyNtMaker::makeCutFlow(const char* name, const char* title)
{
  //TH1F* h = new TH1F(name, title, 9, -0.5, 3.5);
  TH1F* h = new TH1F(name, title, 13, 0., 13.);
  h->GetXaxis()->SetBinLabel(1, "Initial");
  h->GetXaxis()->SetBinLabel(2, "GRL");
  h->GetXaxis()->SetBinLabel(3, "LAr Error");
  h->GetXaxis()->SetBinLabel(4, "Tile Error");
  h->GetXaxis()->SetBinLabel(5, "TTC Veto");
  h->GetXaxis()->SetBinLabel(6, "Good Vertex");
  h->GetXaxis()->SetBinLabel(7, "Hot Spot");
  h->GetXaxis()->SetBinLabel(8, "Bad Jet");
  h->GetXaxis()->SetBinLabel(9, "Bad Muon");
  h->GetXaxis()->SetBinLabel(10, "Cosmic");
  h->GetXaxis()->SetBinLabel(11, ">=1 lep");
  h->GetXaxis()->SetBinLabel(12, ">=2 lep");
  h->GetXaxis()->SetBinLabel(13, "==3 lep");
  return h;
}
/*--------------------------------------------------------------------------------*/
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

/*--------------------------------------------------------------------------------*/
// Main process loop function - This is just an example for testing
/*--------------------------------------------------------------------------------*/
Bool_t SusyNtMaker::Process(Long64_t entry)
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

  if(selectEvent() && m_fillNt){
    m_outTree->Fill(); 
    n_evt_saved++;
  }

  return kTRUE;
}

/*--------------------------------------------------------------------------------*/
// The Terminate() function is the last function to be called during
// a query. It always runs on the client, it can be used to present
// the results graphically or save the results to file.
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::Terminate()
{
  // Stop the timer
  m_timer.Stop();

  SusyD3PDAna::Terminate();
  if(m_dbg) cout << "SusyNtMaker::Terminate" << endl;

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
  cout << "  TileErr  " << n_evt_tileErr  << endl;
  cout << "  TTC Veto " << n_evt_ttcVeto << endl;
  cout << "  GoodVtx  " << n_evt_goodVtx << endl;
  //cout << "  LarHole  " << n_evt_larHole << endl;
  cout << "  HotSpot  " << n_evt_hotSpot << endl;
  cout << "  BadJet   " << n_evt_badJet  << endl;
  cout << "  BadMuon  " << n_evt_badMu   << endl;
  cout << "  Cosmic   " << n_evt_cosmic  << endl;
  cout << "  >=1 Lep  " << n_evt_1Lep    << endl;
  cout << "  >=2 Lep  " << n_evt_2Lep    << endl;
  cout << "  ==3 Lep  " << n_evt_3Lep    << endl;
  cout << endl;

  if(m_fillNt){

    // Save the output tree
    m_outTreeFile = m_outTree->GetCurrentFile();
    m_outTreeFile->Write(0, TObject::kOverwrite);
    cout << "susyNt tree saved to " << m_outTreeFile->GetName() << endl;
    m_outTreeFile->Close();

  }

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
  printf(" Number of events saved:     %d \n",n_evt_saved);
  printf("\t Analysis time: Real %d:%02d:%02d, CPU %.3f      \n", hours, min, sec, cpuTime);
  printf("\t Analysis speed [kHz]: %2.3f                     \n",speed);
  printf("---------------------------------------------------\n\n");
}

/*--------------------------------------------------------------------------------*/
// Select event
/*--------------------------------------------------------------------------------*/
bool SusyNtMaker::selectEvent()
{
  if(m_dbg) cout << "selectEvent" << endl;


  m_susyObj.Reset();
  clearObjects();
  m_susyNt.clear();

  // Susy final state
  m_susyFinalState = m_isSusySample? get_finalState(&d3pd.truth) : -1;
  TH1F* h_procCutFlow = m_isSusySample? getProcCutFlow(m_susyFinalState) : 0;
  float w = m_isMC? d3pd.truth.event_weight() : 1;

  // Cut index
  int cut = 0;

  // Use a preprocessor macro to assist in filling the cutflow hists
  #define FillCutFlow()             \
    do{                             \
      h_rawCutFlow->Fill(cut);      \
      h_genCutFlow->Fill(cut, w);   \
      if(m_isSusySample) h_procCutFlow->Fill(cut, w);  \
      cut++;                        \
    } while(0)

  n_evt_initial++;
  FillCutFlow();
 
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
  // Obj Independent checks

  checkEventCleaning();

  // grl
  //if(!passGRL()) return false;
  if((m_cutFlags & ECut_GRL) == 0) return false;
  FillCutFlow();
  n_evt_grl++;

  // larErr
  //if(!passLarErr()) return false;
  if((m_cutFlags & ECut_LarErr) == 0) return false;
  FillCutFlow();
  n_evt_larErr++;

  // tileErr
  if((m_cutFlags & ECut_TileErr) == 0) return false;
  FillCutFlow();
  n_evt_tileErr++;

  // Incomplete TTC event veto
  //if(!passTTCVeto()) return false;
  if((m_cutFlags & ECut_TTC) == 0) return false;
  FillCutFlow();
  n_evt_ttcVeto++;

  // primary vertex cut is actually filtered
  //if(!passGoodVtx()) return false; 
  if((m_cutFlags & ECut_GoodVtx) == 0) return false;
  FillCutFlow();
  n_evt_goodVtx++;

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//  
  // Get Nominal Objects
  
  selectObjects();
  buildMet();
  //evtCheck();
  checkObjectCleaning();


  // These next cuts are not used to filter the SusyNt because they depend on systematics.
  // Instead, they are simply used for the counters, for comparing the cutflow

  // Tile hot spot
  //if(m_evtFlag & PASS_HotSpot)
  if(m_cutFlags & ECut_HotSpot)
  {
    FillCutFlow();
    n_evt_hotSpot++;

    // Bad jet cut
    if(m_cutFlags & ECut_BadJet)
    {
      FillCutFlow();
      n_evt_badJet++;

      // Bad muon veto
      if(m_cutFlags & ECut_BadMuon)
      {
        FillCutFlow();
        n_evt_badMu++;

        // Cosmic muon veto
        if(m_cutFlags & ECut_Cosmic)
        {
          FillCutFlow();
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
          uint nSigLep = m_sigElectrons.size() + m_sigMuons.size();
          //cout << "nSigLep " << nSigLep << endl;
          if(nSigLep >= 1){
            FillCutFlow();
            n_evt_1Lep++;
            if(nSigLep >= 2){
              FillCutFlow();
              n_evt_2Lep++;
              if(nSigLep == 3){
                FillCutFlow();
                n_evt_3Lep++;
                //cout << "Event " << d3pd.evt.EventNumber() << endl;
              }
            }
          }
        }
      }
    }
  }

  // Match the triggers
  // Will this work for systematic leptons?
  // I think so.
  matchTriggers();

  // Setup reco truth matching
  if(m_isMC){
    m_recoTruthMatch = RecoTruthMatch(0.1, d3pd.truth.channel_number(), 
                                      d3pd.truth.n(), 
                                      d3pd.truth.barcode(), 
                                      d3pd.truth.status(), 
                                      d3pd.truth.pdgId(), 
                                      d3pd.truth.parents(),
                                      d3pd.truth.children(),
                                      d3pd.truth.pt(),
                                      d3pd.truth.eta(),
                                      d3pd.truth.phi(),
                                      d3pd.truth.m(), 0);
  }

  if(m_fillNt){

    // This will fill the pre selected
    // objects prior to overlap removal
    fillNtVars();

    // If it is mc and option for sys is set
    if(m_isMC && m_sys) doSystematic(); 
    
    // TODO: add a command line option for controlling this filtering, 
    // so that we don't keep committing conflicting changes...

    // For filling the output tree, require at least 2 pre-selected leptons (baseline before OR)
    // Now counting taus as well!
    if((m_susyNt.ele()->size() + m_susyNt.muo()->size() + m_susyNt.tau()->size()) < 2)  return false;
  
    // For Fake studies
    //if((m_susyNt.ele()->size() + m_susyNt.muo()->size()) < 1)  return false;

  }
  
  return true;
}

/*--------------------------------------------------------------------------------*/
// Fill SusyNt variables
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillNtVars()
{
  fillEventVars();
  fillLeptonVars();
  fillJetVars();
  fillTauVars();
  fillMetVars();
  fillPhotonVars();
}

/*--------------------------------------------------------------------------------*/
// Fill Event variables
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillEventVars()
{
  Susy::Event* evt = m_susyNt.evt();

  evt->run              = d3pd.evt.RunNumber();
  evt->event            = d3pd.evt.EventNumber();
  evt->lb               = d3pd.evt.lbn();
  evt->stream           = m_stream;

  evt->isMC             = m_isMC;
  evt->mcChannel        = m_isMC? d3pd.truth.channel_number() : 0;
  evt->w                = m_isMC? d3pd.truth.event_weight()   : 1;

  evt->larError         = d3pd.evt.larError();

  evt->nVtx             = getNumGoodVtx();
  evt->avgMu            = d3pd.evt.averageIntPerXing();

  evt->hfor             = m_isMC? getHFORDecision() : -1;

  // SUSY final state
  //evt->susyFinalState   = m_isSusySample? get_finalState(&d3pd.truth) : -1;
  evt->susyFinalState   = m_susyFinalState;

  evt->passMllForAlpgen = m_isMC? PassMllForAlpgen(&d3pd.truth) : true;

  evt->trigFlags        = m_evtTrigFlags;

  evt->wPileup          = m_isMC? getPileupWeight() : 1;
  evt->wPileupAB3       = m_isMC? getPileupWeightAB3() : 1;
  evt->wPileupAB        = m_isMC? getPileupWeightAB() : 1;
  evt->xsec             = m_isMC? getXsecWeight() : 1;
  evt->lumiSF           = m_isMC? getLumiWeight() : 1;             
  evt->sumw             = m_isMC? m_sumw : 1;

  evt->pdf_id1          = m_isMC? d3pd.gen.pdf_id1()->at(0)    : 0;
  evt->pdf_id2          = m_isMC? d3pd.gen.pdf_id2()->at(0)    : 0;
  evt->pdf_x1           = m_isMC? d3pd.gen.pdf_x1()->at(0)     : 0;
  evt->pdf_x2           = m_isMC? d3pd.gen.pdf_x2()->at(0)     : 0;
  evt->pdf_scale        = m_isMC? d3pd.gen.pdf_scale()->at(0)  : 0;

  evt->pdfSF            = m_isMC? getPDFWeight8TeV() : 1;

  m_susyNt.evt()->cutFlags[NtSys_NOM] = m_cutFlags;
}

/*--------------------------------------------------------------------------------*/
// Fill lepton variables
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillLeptonVars()
{
  // loop over preselected leptons and fill the output tree
  for(uint iLep=0; iLep < m_preLeptons.size(); iLep++){
    const LeptonInfo* lep = & m_preLeptons[iLep];
    if(lep->isElectron()) fillElectronVars(lep);
    else fillMuonVars(lep);
  }
}
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillElectronVars(const LeptonInfo* lepIn)
{
  if(m_dbg) cout << "fillElectronVars" << endl;
  m_susyNt.ele()->push_back( Susy::Electron() );
  Susy::Electron* eleOut = & m_susyNt.ele()->back();
  const ElectronElement* element = lepIn->getElectronElement();

  // LorentzVector
  const TLorentzVector* lv = lepIn->lv();
  float pt  = lv->Pt() / GeV;
  float eta = lv->Eta();
  float phi = lv->Phi();
  float m   = lv->M() / GeV;

  eleOut->SetPtEtaPhiM(pt, eta, phi, m);
  eleOut->pt            = pt;
  eleOut->eta           = eta;
  eleOut->phi           = phi;
  eleOut->m             = m;

  // TODO: clean this up, group things together, etc.

  eleOut->ptcone20      = element->ptcone20()/GeV;
  eleOut->ptcone30      = element->ptcone30()/GeV;
  eleOut->q             = element->charge();
  eleOut->mcType        = m_isMC? element->type() : 0;
  eleOut->mcOrigin      = m_isMC? element->origin() : 0;
  eleOut->clusEta       = element->cl_eta();
  eleOut->clusE         = element->cl_E();

  // Check for charge flip
  eleOut->isChargeFlip   = m_isMC? m_recoTruthMatch.isChargeFlip(*lv, element->charge()) : false;
  eleOut->truthMatchType = m_isMC? m_recoTruthMatch.fakeType(*lv, element->origin(), element->type()) : -1;

  // Need to recalculate these variables for p1032 only
  if(m_d3pdTag >= D3PD_p1181){
    eleOut->mediumPP      = element->mediumPP();
    eleOut->tightPP       = element->tightPP();
  }
  else{
    double DEmaxs1 = 0;
    double rHad = 0;
    double rHad1 = 0;
    double ele_et = element->cl_E()/cosh(element->etas2());
    if(element->emaxs1()+element->Emax2() !=0)
      DEmaxs1 = (element->emaxs1()-element->Emax2())/(element->emaxs1()+element->Emax2());
    if(ele_et !=0){
      rHad = element->Ethad()/ele_et;
      rHad1 = element->Ethad1()/ele_et;
    }
  
    eleOut->mediumPP = m_susyObj.ST_isMediumPlusPlus(element->etas2(), ele_et, element->f3(),
						     rHad, rHad1, element->reta(), element->weta2(),
						     element->f1(), element->wstot(), DEmaxs1, 
						     element->deltaeta1(), 
						     element->trackd0_physics(),
						     element->TRTHighTOutliersRatio(), 
						     element->nTRTHits(), element->nTRTOutliers(),
						     element->nSiHits(), 
						     element->nSCTOutliers()+
						     element->nPixelOutliers(), element->nPixHits(),
						     element->nPixelOutliers(), element->nBLHits(),
						     element->nBLayerOutliers(),
						     element->expectHitInBLayer());
    eleOut->tightPP = m_susyObj.ST_isTightPlusPlus(element->etas2(), ele_et, element->f3(),
						   rHad, rHad1, element->reta(), element->weta2(),
						   element->f1(), element->wstot(), DEmaxs1, 
						   element->deltaeta1(), 
						   element->trackd0_physics(),
						   element->TRTHighTOutliersRatio(), 
						   element->nTRTHits(), element->nTRTOutliers(),
						   element->nSiHits(), 
						   element->nSCTOutliers()+
						   element->nPixelOutliers(), element->nPixHits(),
						   element->nPixelOutliers(), element->nBLHits(),
						   element->nBLayerOutliers(),
						   element->expectHitInBLayer(),
						   element->cl_E()*fabs(element->trackqoverp()),
						   element->deltaphi2(), 
						   m_susyObj.GetConvBit(element->isEM()));
  }

  eleOut->d0            = element->trackd0pv();
  eleOut->errD0         = element->tracksigd0pv();
  eleOut->z0            = element->trackz0pv();
  eleOut->errZ0         = element->tracksigz0pv();

  eleOut->d0Unbiased    = element->trackIPEstimate_d0_unbiasedpvunbiased();
  eleOut->errD0Unbiased = element->trackIPEstimate_sigd0_unbiasedpvunbiased();
  eleOut->z0Unbiased    = element->trackIPEstimate_z0_unbiasedpvunbiased();
  eleOut->errZ0Unbiased = element->trackIPEstimate_sigz0_unbiasedpvunbiased();

  // Get d0 from track - old procedure before the d3pd branches were directly available
  //int trkIdx            = get_electron_track( &d3pd.ele, lepIn->idx(), &d3pd.trk );
  //if(trkIdx!=-99){
    //eleOut->d0          = d3pd.trk.d0_wrtPV()->at(trkIdx);
    //eleOut->errD0       = sqrt(d3pd.trk.cov_d0_wrtPV()->at(trkIdx));
  //}

  // New iso variables!! 
  // Corrected topo iso is available in the susy d3pd, apparently calculated using the cluster E.
  // However, the CaloIsoCorrection header says to use the energy after scaling/smearing...
  // So, which should we use?
  // TODO: come back to this, open a discussion somewhere.
  // For now, I will just use the cluster E, which I suspect people will use anyway (even if mistaken)

  // Corrected etcone has Pt and ED corrections
  //eleOut->etcone30Corr  = CaloIsoCorrection::GetPtEDCorrectedIsolation(element->Etcone40(), element->Etcone40_ED_corrected(), lv->E(), element->etas2(), 
  //                                                                     element->etap(), element->cl_eta(), 0.3, m_isMC, element->Etcone30())/GeV;
  eleOut->etcone30Corr  = CaloIsoCorrection::GetPtEDCorrectedIsolation(element->Etcone40(), element->Etcone40_ED_corrected(), element->cl_E(), 
                                                                       element->etas2(), element->etap(), element->cl_eta(), 0.3, m_isMC, 
                                                                       element->Etcone30())/GeV;

  // Corrected topoEtcone has Pt and ED corrections.  Use D3PD branch for now
  //float topo            = CaloIsoCorrection::GetPtEDCorrectedTopoIsolation(element->ED_median(), element->cl_E(), element->etas2(), element->etap(), 
  //                                                                         element->cl_eta(), 0.3, m_isMC, element->topoEtcone30())/GeV;
  eleOut->topoEtcone30Corr      = element->topoEtcone30_corrected()/GeV;
  
  // Trigger flags
  eleOut->trigFlags     = m_eleTrigFlags[ lepIn->idx() ];

  // Efficiency scale factor.  For now, use tightPP if electrons is tightPP, otherwise mediumPP
  int set               = eleOut->tightPP? 7 : 6;
  //int rel               = m_isAF2 ? 9 : 8;
  eleOut->effSF         = m_isMC? m_susyObj.GetSignalElecSF   ( element->cl_eta(), lepIn->lv()->Pt(), set ) : 1;
  eleOut->errEffSF      = m_isMC? m_susyObj.GetSignalElecSFUnc( element->cl_eta(), lepIn->lv()->Pt(), set ) : 1;

  // Do we need this??
  eleOut->idx           = lepIn->idx();
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillMuonVars(const LeptonInfo* lepIn)
{
  if(m_dbg) cout << "fillMuonVars" << endl;
  m_susyNt.muo()->push_back( Susy::Muon() );
  Susy::Muon* muOut = & m_susyNt.muo()->back();
  const MuonElement* element = lepIn->getMuonElement();

  // need truthMuon for type and origin - not anymore!
  //const TruthMuonElement* trueMuon = m_isMC? getMuonTruth( &d3pd.muo, lepIn->idx(), &d3pd.truthMu ) : 0;

  // LorentzVector
  const TLorentzVector* lv = lepIn->lv();
  float pt  = lv->Pt() / GeV;
  float eta = lv->Eta();
  float phi = lv->Phi();
  float m   = lv->M() / GeV;
  muOut->SetPtEtaPhiM(pt, eta, phi, m);
  muOut->pt  = pt;
  muOut->eta = eta;
  muOut->phi = phi;
  muOut->m   = m;

  muOut->q              = element->charge();
  muOut->ptcone20       = element->ptcone20()/GeV;
  muOut->ptcone30       = element->ptcone30()/GeV;
  muOut->etcone30       = element->etcone30()/GeV;

  muOut->ptcone30ElStyle= m_d3pdTag>=D3PD_p1181? element->ptcone30_trkelstyle()/GeV : 0;

  muOut->d0             = element->d0_exPV();
  muOut->errD0          = sqrt(element->cov_d0_exPV());
  muOut->z0             = element->z0_exPV();
  muOut->errZ0          = sqrt(element->cov_z0_exPV());

  muOut->d0Unbiased     = element->trackIPEstimate_d0_unbiasedpvunbiased();
  muOut->errD0Unbiased  = element->trackIPEstimate_sigd0_unbiasedpvunbiased();
  muOut->z0Unbiased     = element->trackIPEstimate_z0_unbiasedpvunbiased();
  muOut->errZ0Unbiased  = element->trackIPEstimate_sigz0_unbiasedpvunbiased();

  muOut->isCombined     = element->isCombinedMuon();

  muOut->truthMatchType = m_isMC? m_recoTruthMatch.fakeType(*lv, element->origin(), element->type()) : -1;
  
  // theta_exPV.  Not sure if necessary.
  muOut->thetaPV        = element->theta_exPV();

  //muOut->mcType         = trueMuon? trueMuon->type()   : 0;
  //muOut->mcOrigin       = trueMuon? trueMuon->origin() : 0;
  muOut->mcType         = m_isMC? element->type()   : 0;
  muOut->mcOrigin       = m_isMC? element->origin() : 0;

  muOut->trigFlags      = m_muoTrigFlags[ lepIn->idx() ];

  muOut->effSF          = m_susyObj.GetSignalMuonSF(lepIn->idx());
  muOut->errEffSF       = m_susyObj.GetSignalMuonSFUnc(lepIn->idx(), SystErr::MEFFUP);
  
  // Do we need this??
  muOut->idx            = lepIn->idx();
}

/*--------------------------------------------------------------------------------*/
// Fill jet variables
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillJetVars()
{
  if(m_dbg) cout << "fillJetVars" << endl;
  // Loop over selected jets and fill output tree
  for(uint iJet=0; iJet<m_preJets.size(); iJet++){
    int jetIndex = m_preJets[iJet];  

    fillJetVar(jetIndex);
  }
}
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillJetVar(int jetIdx)
{

  const JetElement* element = & d3pd.jet[jetIdx];
  m_susyNt.jet()->push_back( Susy::Jet() );
  Susy::Jet* jetOut = & m_susyNt.jet()->back();
  
  const TLorentzVector* lv = & m_susyObj.GetJetTLV(jetIdx);
  float pt  = lv->Pt() / GeV;
  float eta = lv->Eta();
  float phi = lv->Phi();
  float m   = lv->M() / GeV;
  jetOut->SetPtEtaPhiM(pt, eta, phi, m);
  jetOut->pt  = pt;
  jetOut->eta = eta;
  jetOut->phi = phi;
  jetOut->m   = m;
  
  jetOut->idx           = jetIdx;
  jetOut->jvf           = element->jvtxf();
  jetOut->truthLabel    = m_isMC? element->flavor_truth_label() : 0;

  // truth jet matching
  jetOut->matchTruth    = m_isMC? matchTruthJet(jetIdx) : false;
  //cout << "matchTruth: " << jetOut->matchTruth << endl;

  // btag weights
  jetOut->sv0           = element->flavor_weight_SV0();
  jetOut->combNN        = element->flavor_weight_JetFitterCOMBNN();
  jetOut->mv1           = element->flavor_weight_MV1();

}

/*--------------------------------------------------------------------------------*/
// Fill Photon variables
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillPhotonVars()
{
  if(m_dbg) cout << "fillPhotonVars" << endl;

  // Loop over photons
  for(uint iPh=0; iPh<m_sigPhotons.size(); iPh++){
    int phIndex = m_sigPhotons[iPh];  

    fillPhotonVar(phIndex);
  }
}
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillPhotonVar(int phIdx)
{
  if(m_dbg) cout << "fillPhotonVars" << endl;
  m_susyNt.pho()->push_back( Susy::Photon() );
  Susy::Photon* phoOut = & m_susyNt.pho()->back();
  const PhotonElement* element = & d3pd.pho[phIdx];


  // Set TLV
  const TLorentzVector* phTLV = & m_susyObj.GetPhotonTLV(phIdx);
  float pt  = phTLV->Pt() / GeV;
  float E   = phTLV->E()  / GeV;
  float eta = phTLV->Eta();
  float phi = phTLV->Phi();
  
  phoOut->SetPtEtaPhiE(pt, eta, phi, E);
  phoOut->pt  = pt;
  phoOut->eta = eta;
  phoOut->phi = phi;
  phoOut->m   = 0.;
  
  // Save conversion info
  phoOut->isConv = element->isConv();

  // Miscellaneous 
  phoOut->idx    = phIdx;
}

/*--------------------------------------------------------------------------------*/
// Fill Tau variables
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillTauVars()
{
  if(m_dbg) cout << "fillTauVars" << endl;

  // Loop over selected taus
  for(uint iTau=0; iTau<m_preTaus.size(); iTau++){
    int tauIdx = m_preTaus[iTau];  

    fillTauVar(tauIdx);
  }
}
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillTauVar(int tauIdx)
{
  if(m_dbg) cout << "fillTauVars" << endl;
  m_susyNt.tau()->push_back( Susy::Tau() );
  Susy::Tau* tauOut = & m_susyNt.tau()->back();
  const TauElement* element = & d3pd.tau[tauIdx];

  // Set TLV
  const TLorentzVector* tauLV = & m_tauLVs.at(tauIdx);
  float pt  = tauLV->Pt() / GeV;
  float eta = tauLV->Eta();
  float phi = tauLV->Phi();
  float m   = tauLV->M() / GeV;
  
  tauOut->SetPtEtaPhiM(pt, eta, phi, m);
  tauOut->pt    = pt;
  tauOut->eta   = eta;
  tauOut->phi   = phi;
  tauOut->m     = m;

  tauOut->q                     = element->charge();
  tauOut->author                = element->author();
  tauOut->nTrack                = element->numTrack();
  tauOut->eleBDT                = element->BDTEleScore();
  tauOut->jetBDT                = element->BDTJetScore();

  tauOut->jetBDTSigLoose        = element->JetBDTSigLoose();
  tauOut->jetBDTSigMedium       = element->JetBDTSigMedium();
  tauOut->jetBDTSigTight        = element->JetBDTSigTight();
  tauOut->eleBDTLoose           = element->EleBDTLoose();
  tauOut->eleBDTMedium          = element->EleBDTMedium();
  tauOut->eleBDTTight           = element->EleBDTTight();

  tauOut->muonVeto              = element->muonVeto();
  tauOut->trueTau               = m_isMC? element->trueTauAssocSmall_matched() : false;

  tauOut->trigFlags             = m_tauTrigFlags[tauIdx];
  
  tauOut->idx   = tauIdx;
}

/*--------------------------------------------------------------------------------*/
// Fill MET variables
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::fillMetVars(SusyNtSys sys)
{
  if(m_dbg) cout << "fillMetVars" << endl;

  // Just fill the lv for now
  double Et  = m_met.Et()/GeV;
  double phi = m_met.Phi(); 

  //double px = m_met.Px()/GeV;
  //double py = m_met.Py()/GeV;
  //double pz = m_met.Pz()/GeV;
  //double E  = m_met.E()/GeV;
  
  m_susyNt.met()->push_back( Susy::Met() );
  Susy::Met* metOut = & m_susyNt.met()->back();
  metOut->Et  = Et;
  metOut->phi = phi;
  metOut->sys = sys;  

  // MET comp terms 
  // Need to save these for the MET systematics as well.
  // Use the sys enum to determine which argument to pass to SUSYTools
  METUtil::Systematics metSys = METUtil::None;

  // I guess these are the only ones we need to specify, the ones specified in SUSYTools...
  // All the rest should be automatic (I think), e.g. JES
  if(sys == NtSys_SCALEST_UP) metSys = METUtil::ScaleSoftTermsUp;
  else if(sys == NtSys_SCALEST_DN) metSys = METUtil::ScaleSoftTermsDown;
  else if(sys == NtSys_RESOST) metSys = METUtil::ResoSoftTermsUp;

  // Save the MET terms
  metOut->refEle        = m_susyObj.computeMETComponent(METUtil::RefEle, metSys).Mod()/GeV;
  metOut->refMuo        = m_susyObj.computeMETComponent(METUtil::MuonTotal, metSys).Mod()/GeV;
  metOut->refJet        = m_susyObj.computeMETComponent(METUtil::RefJet, metSys).Mod()/GeV;
  metOut->refGamma      = m_susyObj.computeMETComponent(METUtil::RefGamma, metSys).Mod()/GeV;
  metOut->softJet       = m_susyObj.computeMETComponent(METUtil::SoftJets, metSys).Mod()/GeV;
  metOut->refCell       = m_susyObj.computeMETComponent(METUtil::CellOutEflow, metSys).Mod()/GeV;
}

/*--------------------------------------------------------------------------------*/
// Handle Systematic 
/*--------------------------------------------------------------------------------*/
void SusyNtMaker::doSystematic()
{
  // Loop over the systematics:
  // Start at 1, nominal saved
  for(int i = 1; i < NtSys_N; i++){

    SusyNtSys sys = (SusyNtSys) i;

    // Reset Objects 
    m_susyObj.Reset();
    clearObjects();

    selectObjects(sys);
    buildMet(sys);                   
    //evtCheck();

    checkEventCleaning();
    checkObjectCleaning();

    // Lepton Specific sys    
    if( isElecSys(sys) )
      saveElectronSF(sys);
    else if( isMuonSys(sys) )
      saveMuonSF(sys);
    else if( isJetSys(sys) )
      saveJetSF(sys);


    // Fill the Met for this sys
    fillMetVars(sys);

    // Add the event flag for this event
    //addEventFlag(sys, m_evtFlag);
    m_susyNt.evt()->cutFlags[sys] = m_cutFlags;

  }// end loop over systematic
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::saveElectronSF(SusyNtSys sys)
{
  // loop over preselected leptons and fill the output tree
  for(uint iLep=0; iLep < m_preLeptons.size(); iLep++){
    const LeptonInfo* lep = & m_preLeptons[iLep];
    if(!lep->isElectron()) continue;
    
    bool match = false;

    // Match Indices and save sf
    for(uint iE=0; iE<m_susyNt.ele()->size(); ++iE){
      Susy::Electron* ele = & m_susyNt.ele()->at(iE);
      if( ele->idx == lep->idx() ){
	match = true;
	
	float sf = lep->lv()->E() / GeV / ele->E();
	//if(sys == NtSys_EES_UP)      ele->ees_up = sf;
	//else if(sys == NtSys_EES_DN) ele->ees_dn = sf;
	if(sys == NtSys_EES_Z_UP)        ele->ees_z_up = sf;
	else if(sys == NtSys_EES_Z_DN)   ele->ees_z_dn = sf;
	else if(sys == NtSys_EES_MAT_UP) ele->ees_mat_up = sf;
	else if(sys == NtSys_EES_MAT_DN) ele->ees_mat_dn = sf;
	else if(sys == NtSys_EES_PS_UP)  ele->ees_ps_up = sf;
	else if(sys == NtSys_EES_PS_DN)  ele->ees_ps_dn = sf;
	else if(sys == NtSys_EES_LOW_UP) ele->ees_low_up = sf;
	else if(sys == NtSys_EES_LOW_DN) ele->ees_low_dn = sf;
	else if(sys == NtSys_EER_UP)     ele->eer_up = sf;
	else if(sys == NtSys_EER_DN)     ele->eer_dn = sf;
	
      }// end if electron matches
    }// end loop over electrons in susyNt
    
    // The dreaded case...
    if( !match )
      addMissingElectron(lep, sys);

  }// end loop over leptons
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::saveMuonSF(SusyNtSys sys)
{
  // loop over preselected leptons and fill the output tree
  for(uint iLep=0; iLep < m_preLeptons.size(); iLep++){
    const LeptonInfo* lep = & m_preLeptons[iLep];
    if(lep->isElectron()) continue;
    
    bool match = false;

    // Match Indices and save sf
    for(uint iM=0; iM<m_susyNt.muo()->size(); ++iM){
      Susy::Muon* mu = & m_susyNt.muo()->at(iM);
      if( mu->idx == lep->idx() ){
	match = true;
	
	float sf = lep->lv()->E() / GeV / mu->E();
	if(sys == NtSys_MS_UP)      mu->ms_up = sf;
	else if(sys == NtSys_MS_DN) mu->ms_dn = sf;
	else if(sys == NtSys_ID_UP) mu->id_up = sf;
	else if(sys == NtSys_ID_DN) mu->id_dn = sf;

      }// end if muon matches
    }// end loop over muons in SusyNt
    
    // The dreaded case...
    if(!match) addMissingMuon(lep, sys);

  }// end loop over leptons
      
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::saveJetSF(SusyNtSys sys)
{
  // Loop over selected jets and fill output tree
  for(uint iJet=0; iJet<m_preJets.size(); iJet++){
    uint jetIdx = m_preJets[iJet];
    //const JetElement* element = & d3pd.jet[jetIdx];
    const TLorentzVector jetTLV = m_susyObj.GetJetTLV(jetIdx);
    
    bool match = false;

    // Check Indices and save sf
    for(uint iJ=0; iJ<m_susyNt.jet()->size(); ++iJ){
      Susy::Jet* jet = & m_susyNt.jet()->at(iJ);
      if( jet->idx == jetIdx ){
	match = true;
	
	float sf = jetTLV.E() / GeV / jet->E();
	if(sys == NtSys_JES_UP)      jet->jes_up = sf;
	else if(sys == NtSys_JES_DN) jet->jes_dn = sf;
	else if(sys == NtSys_JER)    jet->jer = sf;

      }// end if the leptons match	
    }// end loop over what we have saved
    
    // The dreaded case...
    if(!match) addMissingJet(jetIdx, sys);
    
  }// end loop over jets in pre-jets

}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::addMissingElectron(const LeptonInfo* lep, SusyNtSys sys)
{
  //cout<<"Adding an electron that was missing!"<<endl;

  // This electron did not pass nominal cuts, and therefore
  // needs to be added, but with the correct TLV
  
  TLorentzVector tlv_sys =m_susyObj.GetElecTLV(lep->idx());

  D3PDReader::ElectronD3PDObject* e = & d3pd.ele;
  int index      = lep->idx();
  float cl_eta   = e->cl_eta()->at(index);
  float cl_phi   = e->cl_phi()->at(index);
  float cl_E     = e->cl_E()->at(index);
  float trk_eta  = e->tracketa()->at(index);
  float trk_phi  = e->trackphi()->at(index);
  float nPixHits = e->nPixHits()->at(index);
  float nSCTHits = e->nSCTHits()->at(index);
  bool isData    = !m_isMC;

  // Reset the Nominal TLV
  m_susyObj.SetElecTLV(index, cl_eta, cl_phi, cl_E, trk_eta, trk_phi, nPixHits, nSCTHits, isData, d3pd.evt.RunNumber(), SystErr::NONE);

  // Now push it back onto to susyNt
  fillElectronVars(lep);
  
  // Set the sf
  Susy::Electron* ele = & m_susyNt.ele()->back();
  float sf = tlv_sys.E() / GeV / ele->E();
  if(sys == NtSys_EES_Z_UP)        ele->ees_z_up = sf;
  else if(sys == NtSys_EES_Z_DN)   ele->ees_z_dn = sf;
  else if(sys == NtSys_EES_MAT_UP) ele->ees_mat_up = sf;
  else if(sys == NtSys_EES_MAT_DN) ele->ees_mat_dn = sf;
  else if(sys == NtSys_EES_PS_UP)  ele->ees_ps_up = sf;
  else if(sys == NtSys_EES_PS_DN)  ele->ees_ps_dn = sf;
  else if(sys == NtSys_EES_LOW_UP) ele->ees_low_up = sf;
  else if(sys == NtSys_EES_LOW_DN) ele->ees_low_dn = sf;
  else if(sys == NtSys_EER_UP)     ele->eer_up = sf;
  else if(sys == NtSys_EER_DN)     ele->eer_dn = sf;
  
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::addMissingMuon(const LeptonInfo* lep, SusyNtSys sys)
{
  //cout<<"Adding a muon that was missing!"<<endl;

  // This muon did not pass nominal cuts, and therefore
  // needs to be added, but with the correct TLV
  
  TLorentzVector tlv_sys = m_susyObj.GetMuonTLV(lep->idx());

  D3PDReader::MuonD3PDObject* m = & d3pd.muo;
  int index               = lep->idx();
  float pt                = m->pt()->at(index);
  float eta               = m->eta()->at(index);
  float phi               = m->phi()->at(index);
  float E                 = m->E()->at(index);
  float me_qoverp_exPV    = m->me_qoverp_exPV()->at(index);
  float id_qoverp_exPV    = m->id_qoverp_exPV()->at(index);
  float me_theta_exPV     = m->me_theta_exPV()->at(index);
  float id_theta_exPV     = m->id_theta_exPV()->at(index);
  float charge            = m->charge()->at(index);
  int isCombined          = m->isCombinedMuon()->at(index);
  int isSegTag            = m->isSegmentTaggedMuon()->at(index);
  bool isData             = !m_isMC;

  // Reset the Nominal TLV
  //m_susyObj.SetMuonTLV(index, pt, eta, phi, E, me_qoverp_exPV, id_qoverp_exPV, me_theta_exPV, id_theta_exPV, isCombined, isData, SystErr::NONE);
  m_susyObj.SetMuonTLV(index, pt, eta, phi, E, me_qoverp_exPV, id_qoverp_exPV, me_theta_exPV, 
                       id_theta_exPV, charge, isCombined, isSegTag, isData, SystErr::NONE);
  
  // Now push it back onto to susyNt
  fillMuonVars(lep);
  
  // Set the sf
  Susy::Muon* mu = & m_susyNt.muo()->back();
  float sf = tlv_sys.E() / GeV / mu->E();
  if(sys == NtSys_MS_UP) mu->ms_up = sf;
  else if(sys == NtSys_MS_DN) mu->ms_dn = sf;
  else if(sys == NtSys_ID_UP) mu->id_up = sf;
  else if(sys == NtSys_ID_DN) mu->id_dn = sf;
}

/*--------------------------------------------------------------------------------*/
void SusyNtMaker::addMissingJet(int index, SusyNtSys sys)
{
  // Get the tlv with the sys
  TLorentzVector tlv_sys = m_susyObj.GetJetTLV(index);
  
  D3PDReader::JetD3PDObject* jets = & d3pd.jet;

  float pt  = jets->pt()->at(index);
  float eta = jets->eta()->at(index);
  float phi = jets->phi()->at(index);
  float E   = jets->E()->at(index);
  
  // Reset TLV
  m_susyObj.SetJetTLV(index, pt, eta, phi, E);

  // Fill the Jet vars for this guy
  fillJetVar(index);

  // Set SF
  Susy::Jet* j = & m_susyNt.jet()->back();
  float sf = tlv_sys.E() / GeV / j->E();
  if(sys == NtSys_JER)         j->jer = sf;
  else if(sys == NtSys_JES_UP) j->jes_up = sf;
  else if(sys == NtSys_JES_DN) j->jes_dn = sf;
}

#undef GeV
