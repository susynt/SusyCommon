#ifndef SusyCommon_SusyD3PDAna_h
#define SusyCommon_SusyD3PDAna_h


#include <iostream>

#include "GoodRunsLists/TGoodRunsList.h"
#include "GoodRunsLists/TGoodRunsListReader.h"
#include "SUSYTools/SUSYObjDef.h"
#include "SUSYTools/FakeMetEstimator.h"
#include "SUSYTools/SUSYCrossSection.h"
#include "PileupReweighting/TPileupReweighting.h"
#include "MultiLep/LeptonInfo.h"

#ifdef USEPDFTOOL
#include "MultiLep/PDFTool.h"
#endif

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

    // Sample name - used to set isMC flag
    TString sample() { return m_sample; }
    void setSample(TString s) { m_sample = s; }

    //
    // Object selection
    // Selected leptons have kinematic and cleaning cuts (no overlap removal)
    // Baseline leptons = selected + overlap removed
    // 

    // Full object selection
    void selectObjects(SusyNtSys sys = NtSys_NOM){
      selectBaselineObjects(sys);
      selectSignalObjects();
    }
    void selectBaselineObjects(SusyNtSys sys = NtSys_NOM);
    void selectSignalObjects();
    void performOverlapRemoval();
    void selectSignalPhotons();

    // MissingEt
    void buildMet(SusyNtSys sys = NtSys_NOM);

    // Clear object selection
    void clearObjects();

    // Count number of good vertices
    uint getNumGoodVtx();

    //
    // Event cleaning
    //

    // grl
    void setGRLFile(TString fileName) { m_grlFileName = fileName; }
    bool passGRL(){ return m_isMC || m_grl.HasRunLumiBlock(d3pd.evt.RunNumber(), d3pd.evt.lbn()); }
    // lar error
    bool passLarErr(){ return m_isMC || (d3pd.evt.larError()==0); }
    // lar hole veto
    bool passLarHoleVeto();
    // tile hot spot
    bool passTileHotSpot();
    // bad jet
    bool passBadJet();
    // good vertex
    bool passGoodVtx();
    // bad muon veto
    bool passBadMuon();
    // cosmic veto
    bool passCosmic();

    // Poorly named, checks some event cleaning cuts
    void evtCheck();

    //
    // Event weighting
    //

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
    // pileup weight, not included in event weight above
    float getPileupWeight();
    float getPileupWeight1fb();
    // PDF reweighting of 7TeV -> 8TeV
    float getPDFWeight8TeV();

    //
    // Running conditions
    //

    // Set sys run
    void setSys(bool sysOn){ m_sys = sysOn; };
    
    // Toggle photon selection
    //void setSavePhotons(bool phOn){ m_savePh = phOn; };
    void setSelectPhotons(bool doIt) { m_selectPhotons = doIt; }

    // Toggle tau selection and overlap removal
    void setSelectTaus(bool doIt) { m_selectTaus = doIt; }
    
    //
    // Trigger - check matching for all baseline leptons
    //
    void resetTriggers(){
      m_evtTrigFlags = 0;
      m_eleTrigFlags.clear();
      m_muoTrigFlags.clear();
    }
    void matchTriggers(){
      fillEventTriggers();
      matchElectronTriggers();
      matchMuonTriggers();
    }
    void fillEventTriggers();
    void matchElectronTriggers();
    bool matchElectronTrigger(const TLorentzVector* lv, std::vector<int>* trigBools);
    void matchMuonTriggers();
    bool matchMuonTrigger(const TLorentzVector* lv, std::vector<int>* trigBools);


    // Debugging method
    void dump();

  protected:

    TString                     m_sample;       // sample name
    DataStream                  m_stream;       // data stream enum, taken from sample name

    //std::string                 m_metCalib;     // Calibration string for MET (and jets)

    //
    // Object collections (usually just vectors of indices)
    //

    // "selected" objects pass kinematic cuts, but no overlap removal applied
    std::vector<int>            m_preElectrons; // selected electrons
    std::vector<int>            m_preMuons;     // selected muons
    std::vector<LeptonInfo>     m_preLeptons;   // selected leptons
    std::vector<int>            m_preJets;      // selected jets
    std::vector<int>            m_preTaus;      // selected taus
    
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

    // Tau TLorentzVectors, because they are not stored in SUSYTools
    std::vector<TLorentzVector> m_tauLVs;

    // MET
    TLorentzVector              m_met;          // fully corrected MET

    uint                        m_evtTrigFlags; // Event trigger flags
    
    // Trigger object matching maps
    // Key: d3pd index, Val: trig bit word
    std::map<int, uint>         m_eleTrigFlags; // electron trigger matching flags
    std::map<int, uint>         m_muoTrigFlags; // muon trigger matching flags
    
    float                       m_lumi;         // normalized luminosity (defaults to 4.7/fb)
    float                       m_sumw;         // sum of mc weights for normalization, must be set by user
    float                       m_xsec;         // optional user cross section, to override susy xsec usage
    bool                        m_sys;          // True if you want sys for MC, must be set by user. 
    //bool                        m_savePh;       // True if want to save photons
    bool                        m_selectPhotons;// Toggle photon selection
    bool                        m_selectTaus;   // Toggle tau selection and overlap removal

    uint                        m_evtFlag;      // Reset after each evt

    //
    // Tools
    //

    SUSYObjDef                  m_susyObj;      // SUSY object definitions

    TString                     m_grlFileName;  // grl file name
    Root::TGoodRunsList         m_grl;          // good runs list

    FakeMetEstimator            m_fakeMetEst;   // fake met estimator for lar hole veto

    Root::TPileupReweighting*   m_pileup;       // pileup reweighting
    Root::TPileupReweighting*   m_pileup1fb;    // pileup reweighting for 2012 A-B3 only

    // The SUSY CrossSectionDB has its own map for retrieving xsec info, but
    // it has a lot of entries so lookup is slow.  Save our own xsec map

    SUSY::CrossSectionDB*                       m_susyXsec;     // SUSY cross section database
    std::map<int,SUSY::CrossSectionDB::Process> m_xsecMap;      // our own xsec map for faster lookup times

    #ifdef USEPDFTOOL
    PDFTool*                    m_pdfTool;      // PDF reweighting tool (In MultiLep pkg)
    #endif

};

#endif
