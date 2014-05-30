#ifndef SusyCommon_SusyD3PDAna_h
#define SusyCommon_SusyD3PDAna_h


#include <iostream>

#include "GoodRunsLists/TGoodRunsList.h"
#include "GoodRunsLists/TGoodRunsListReader.h"
#include "SUSYTools/SUSYObjDef.h"
#include "SUSYTools/FakeMetEstimator.h"
#include "SUSYTools/SUSYCrossSection.h"
#include "SUSYTools/HforToolD3PD.h"
#include "PileupReweighting/TPileupReweighting.h"
#include "LeptonTruthTools/RecoTauMatch.h"


/* #ifdef USEPDFTOOL */
/* #include "MultiLep/PDFTool.h" */
/* #endif */

#include "SusyCommon/LeptonInfo.h"
#include "SusyCommon/SusyD3PDInterface.h"

/*

    SusyD3PDAna - a class for performing object selections and event cleaning on susy d3pds

*/

class SusyD3PDAna : public SusyD3PDInterface
{

  public:

    // Constructor and destructor
    SusyD3PDAna();
    virtual ~SusyD3PDAna();
    
    // Begin is called before looping on entries
    virtual void    Begin(TTree *tree);
    // Main event loop function
    virtual Bool_t  Process(Long64_t entry);
    // Terminate is called after looping is finished
    virtual void    Terminate();

    //
    // Object selection
    // Selected leptons have kinematic and cleaning cuts (no overlap removal)
    // Baseline leptons = selected + overlap removed
    // 

    // Full object selection
    void selectObjects(SusyNtSys sys = NtSys_NOM){
      selectBaselineObjects(sys);
      selectSignalObjects();
      if(m_selectTruth) selectTruthObjects();
    }
    void selectBaselineObjects(SusyNtSys sys = NtSys_NOM);
    void selectSignalObjects();
    void performOverlapRemoval();
    void selectSignalPhotons();
    void selectTruthObjects();

    // MissingEt
    void buildMet(SusyNtSys sys = NtSys_NOM);

    // Clear object selection
    void clearObjects();


    //
    // Trigger - check matching for all baseline leptons
    //
    void resetTriggers(){
      m_evtTrigFlags = 0;
      m_eleTrigFlags.clear();
      m_muoTrigFlags.clear();
      m_tauTrigFlags.clear();
    }
    void matchTriggers(){
      fillEventTriggers();
      matchElectronTriggers();
      matchMuonTriggers();
      matchTauTriggers();
    }
    void fillEventTriggers();
    void matchElectronTriggers();
    bool matchElectronTrigger(const TLorentzVector* lv, std::vector<int>* trigBools);
    void matchMuonTriggers();
    bool matchMuonTrigger(const TLorentzVector* lv, std::vector<int>* trigBools);
    void matchTauTriggers();
    bool matchTauTrigger(const TLorentzVector* lv, std::vector<int>* trigBools);



    //
    // Event cleaning
    //

    // grl
    void setGRLFile(TString fileName) { m_grlFileName = fileName; }
    bool passGRL() { return m_isMC || m_grl.HasRunLumiBlock(d3pd.evt.RunNumber(), d3pd.evt.lbn()); }
    // incomplete TTC event veto
    bool passTTCVeto() { return (d3pd.evt.coreFlags() & 0x40000) == 0; }
    // Tile error
    bool passTileErr() { return m_isMC || (d3pd.evt.tileError()!=2); }
    // lar error
    bool passLarErr() { return m_isMC || (d3pd.evt.larError()!=2); }
    // lar hole veto
    bool passLarHoleVeto();
    // tile hot spot
    bool passTileHotSpot();
    // bad jet
    bool passBadJet();
    // good vertex
    bool passGoodVtx();
    // tile trip
    bool passTileTrip();
    // bad muon veto
    bool passBadMuon();
    // cosmic veto
    bool passCosmic();

    // Sherpa WW fix for radiative b quarks
    // See thread "Diboson MC Truth Discrepancy" atlas-phys-susy-d3pd.cern.ch, Mar 2013
    // UPDATE: now using central SUSYTools method: Sherpa_WW_veto
    /*
    bool isBuggyWWSherpaSample(const int& dsid){
      return dsid==126892 || dsid==157817 || dsid==157818 || dsid==157819;
    }
    bool hasRadiativeBQuark(const std::vector<int>* pdg, const std::vector<int>* status);
    */

    // Event level cleaning cuts
    void checkEventCleaning();
    // Object level cleaning cuts; these depend on sys
    void checkObjectCleaning();

    //
    // Event weighting
    //

    // Full event weight includes generator, xsec, pileup, and lumi weights.
    // Default weight uses A-E lumi.
    // You can supply a different integrated luminosity, 
    // but the the pileup weights will still correspond to A-E.
    float getEventWeight(float lumi = LUMI_A_E);
    // This function will give the MC weight corresponding to the A-B3 unblinded dataset
    // with correct pileup weights 
    float getEventWeightAtoB3();
    // MC weight corresponding to A-B dataset (5.83/fb)
    float getEventWeightAtoB();

    // event weight (xsec*kfac) 
    float getXsecWeight();
    // lumi weight (lumi/sumw) normalized to 4.7/fb
    float getLumiWeight();
    // luminosity to normalize to (in 1/pb)
    void setLumi(float lumi) { m_lumi = lumi; }
    // sum of mc weights for sample
    void setSumw(float sumw) { m_sumw = sumw; }
    // user cross section, overrides susy cross section
    void setXsec(float xsec) { m_xsec = xsec; }
    // user cross section uncert
    void setErrXsec(float err) { m_errXsec = err; }

    // pileup weight for full dataset: currently A-L
    float getPileupWeight();
    float getPileupWeightUp();
    float getPileupWeightDown();
    // pileup weight for A-B3 (1.037/fb)
    float getPileupWeightAB3();
    // pileup weight for A-B (5.83/fb)
    float getPileupWeightAB();
    // pileup weight for A-E HCP dataset
    float getPileupWeightAE();
    // PDF reweighting of 7TeV -> 8TeV
    float getPDFWeight8TeV();

    // Lepton efficiency SF
    float getLepSF(const std::vector<LeptonInfo>& leptons);

    // BTag efficiency SF
    float getBTagSF(const std::vector<int>& jets);


    //
    // Utility methods
    //

    // calculate random run/lb numbers for MC
    void calcRandomRunLB();

    // Mass helpers
    //float Mll();
    //bool isZ();
    //bool hasZ();

    // HF overlap removal decision
    int getHFORDecision();

    // Count number of good vertices
    uint getNumGoodVtx();

    // Match a reco jet to a truth jet
    bool matchTruthJet(int iJet);

    //
    // Running conditions
    //

    // Sample name - used to set isMC flag
    TString sample() { return m_sample; }
    void setSample(TString s) { m_sample = s; }

    // AF2 flag
    void setAF2(bool isAF2=true) { m_isAF2 = isAF2; }

    // Set MC Production flag
    void setMCProduction(MCProduction prod) { m_mcProd = prod; }

    // Set SUSY D3PD tag to know which branches are ok
    void setD3PDTag(D3PDTag tag) { m_d3pdTag = tag; }

    // Set sys run
    void setSys(bool sysOn){ m_sys = sysOn; };
    
    // Toggle photon selection
    void setSelectPhotons(bool doIt) { m_selectPhotons = doIt; }

    // Toggle tau selection and overlap removal
    void setSelectTaus(bool doIt) { m_selectTaus = doIt; }

    // Set-Get truth selection
    void setSelectTruthObjects(bool doIt) { m_selectTruth = doIt; }
    bool getSelectTruthObjects(         ) { return m_selectTruth; }

    // Set MET flavor - at the moment, only STVF and STVF_JVF are available.
    // Anything else will raise an error.
    void setMetFlavor(std::string metFlav);
    void setDoMetMuonCorrection(bool doMetMuCorr) { m_doMetMuCorr = doMetMuCorr; }
    void setDoMetFix(bool doMetFix) { m_doMetFix = doMetFix; }
    //void setUseMetMuons(bool useMetMu) { m_useMetMuons = useMetMu; }

    //
    // Event dumps
    //
    void dumpEvent();
    void dumpBaselineObjects();
    void dumpSignalObjects();


  protected:

    TString                     m_sample;       // sample name
    DataStream                  m_stream;       // data stream enum, taken from sample name
    bool                        m_isAF2;        // flag for ATLFastII samples
    MCProduction                m_mcProd;       // MC production campaign

    bool                        m_isSusySample; // is susy grid sample
    int                         m_susyFinalState;// susy subprocess

    D3PDTag                     m_d3pdTag;      // SUSY D3PD tag

    bool                        m_selectPhotons;// Toggle photon selection
    bool                        m_selectTaus;   // Toggle tau selection and overlap removal
    bool                        m_selectTruth;  // Toggle truth selection

    SUSYMet::met_definition     m_metFlavor;    // MET flavor enum (e.g. STVF, STVF_JVF)
    bool                        m_doMetMuCorr;  // Control MET muon Eloss correction in SUSYTools
    bool                        m_doMetFix;     // Control MET Egamma-jet overlap fix in SUSYTools
    //bool                      m_useMetMuons;  // Use appropriate muons for met

    //
    // Object collections (usually just vectors of indices)
    //

    // "container" objects pass minimal selection cuts
    std::vector<int>            m_contTaus;     // container taus

    // "selected" objects pass kinematic cuts, but no overlap removal applied
    std::vector<int>            m_preElectrons; // selected electrons
    std::vector<int>            m_preMuons;     // selected muons
    std::vector<LeptonInfo>     m_preLeptons;   // selected leptons
    std::vector<int>            m_preJets;      // selected jets
    std::vector<int>            m_preTaus;      // selected taus
    std::vector<int>            m_metMuons;     // selected muons with larger eta cut for met calc.
    
    // "baseline" objects pass selection + overlap removal
    std::vector<int>            m_baseElectrons;// baseline electrons
    std::vector<int>            m_baseMuons;    // baseline muons
    std::vector<LeptonInfo>     m_baseLeptons;  // baseline leptonInfos
    std::vector<int>            m_baseTaus;     // baseline taus
    std::vector<int>            m_baseJets;     // baseline jets

    // "signal" objects pass baseline + signal selection (like iso)
    std::vector<int>            m_sigElectrons; // signal electrons
    std::vector<int>            m_sigMuons;     // signal muons
    std::vector<LeptonInfo>     m_sigLeptons;   // signal leptonInfos
    std::vector<int>            m_sigTaus;      // signal taus
    std::vector<int>            m_sigPhotons;   // signal photons
    std::vector<int>            m_sigJets;      // signal jets

    // MET
    TLorentzVector              m_met;          // fully corrected MET

    // Truth Objects
    std::vector<int>            m_truParticles; // selected truth particles
    std::vector<int>            m_truJets;      // selected truth jets
    TLorentzVector              m_truMet;       // Truth MET

    long long                   m_evtTrigFlags; // Event trigger flags
    
    // Trigger object matching maps
    // Key: d3pd index, Val: trig bit word
    std::map<int, long long>    m_eleTrigFlags; // electron trigger matching flags
    std::map<int, long long>    m_muoTrigFlags; // muon trigger matching flags
    std::map<int, long long>    m_tauTrigFlags; // tau trigger matching flags
    
    //
    // Event quantities
    //

    float                       m_lumi;         // normalized luminosity (defaults to 4.7/fb)
    float                       m_sumw;         // sum of mc weights for normalization, must be set by user
    float                       m_xsec;         // optional user cross section, to override susy xsec usage
    float                       m_errXsec;      // user cross section uncertainty

    uint                        m_mcRun;        // Random run number for MC from pileup tool
    uint                        m_mcLB;         // Random lb number for MC from pileup tool

    bool                        m_sys;          // True if you want sys for MC, must be set by user. 

    uint                        m_cutFlags;     // Event cleaning cut flags

    //
    // Tools
    //

    SUSYObjDef                  m_susyObj;      // SUSY object definitions
    Root::TElectronEfficiencyCorrectionTool* m_eleMediumSFTool;

    TString                     m_grlFileName;  // grl file name
    Root::TGoodRunsList         m_grl;          // good runs list

    FakeMetEstimator            m_fakeMetEst;   // fake met estimator for lar hole veto

    Root::TPileupReweighting*   m_pileup;       // pileup reweighting
    Root::TPileupReweighting*   m_pileup_up;    // pileup reweighting
    Root::TPileupReweighting*   m_pileup_dn;    // pileup reweighting

    Root::TPileupReweighting*   m_pileupAB3;    // pileup reweighting for 2012 A-B3 only
    Root::TPileupReweighting*   m_pileupAB;     // pileup reweighting for 2012 A-B
    Root::TPileupReweighting*   m_pileupAE;     // pileup reweighting for 2012 A-H (HCP dataset)

    // The SUSY CrossSectionDB has its own map for retrieving xsec info, but
    // it has a lot of entries so lookup is slow.  Save our own xsec map

    SUSY::CrossSectionDB*                       m_susyXsec;     // SUSY cross section database
    std::map<int,SUSY::CrossSectionDB::Process> m_xsecMap;      // our own xsec map for faster lookup times

    HforToolD3PD                m_hforTool;     // heavy flavor overlap removal tool

    #ifdef USEPDFTOOL
    PDFTool*                    m_pdfTool;      // PDF reweighting tool (In MultiLep pkg)
    #endif

    //RecoTruthMatch            m_recoTruthMatch;       // Lepton truth matching tool
    RecoTauMatch                m_recoTruthMatch;       // Lepton truth matching tool

};

#endif
