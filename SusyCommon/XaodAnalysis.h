#ifndef SusyCommon_XaodAnalysis_h
#define SusyCommon_XaodAnalysis_h


#include <iostream>

#include "TSelector.h"
#include "TTree.h"

// Infrastructure include(s):
#ifdef ROOTCORE
#   include "xAODRootAccess/Init.h"
#   include "xAODRootAccess/TEvent.h"
#endif // ROOTCORE

#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "LeptonTruthTools/RecoTauMatch.h"

#include "SusyCommon/LeptonInfo.h"
#include "SusyNtuple/SusyDefs.h"

#include "xAODEventInfo/EventInfo.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODCore/ShallowCopy.h"


namespace susy {

///  a class for performing object selections and event cleaning on xaod
class XaodAnalysis : public TSelector
{

  public:
    XaodAnalysis();
    virtual ~XaodAnalysis();

    virtual Bool_t  Process(Long64_t entry);
    virtual void    Terminate();
    virtual void    Init(TTree *tree); ///< Init is called every time a new TTree is attached
    virtual void    SlaveBegin(TTree *tree);

    virtual Bool_t  Notify() { return kTRUE; } /// Called at the first entry of a new file in a chain
    virtual void    SlaveTerminate(){};
    /// Due to ROOT's stupid design, need to specify version >= 2 or the tree will not connect automatically
    virtual Int_t   Version() const { return 2; }
    virtual XaodAnalysis& setDebug(int debugLevel) { m_dbg = debugLevel; return *this; }
    XaodAnalysis& initSusyTools(); ///< initialize SUSYObjDef_xAOD
    bool processingMc12b() const { return m_mcProd == MCProd_MC12b; }
    /// access the event info
    /**
      \todo: all these xaod getters should be cached (if they are a bottleneck) DG-2014-08-16 to be checked
     */
    virtual const xAOD::EventInfo* xaodEventInfo();
    /// access the default collection of muons from the D3PDReader
    /**
       By default this function returns a pointer to
       mu_staco. However, if we always call this function (rather than
       accessing directly m_event.mu_staco), one can decide to
       override this member function, and easily switch to another
       muon collection.
       In addition we can call here all the functions needed to access
       the auxilliary information.
     */
    virtual xAOD::MuonContainer* xaodMuons();
    /// access the default collection of electrons from the D3PDReader
    /**
       By default returns m_event.el; for its motivation, see XaodAnalysis::xaodMuons().
       \todo In this case there might be some ambiguity to be sorted out when calling SUSYObjDef::GetMET().
     */
    virtual xAOD::ElectronContainer* xaodElectrons();
    /// access the default collection of taus from the D3PDReader
    /**
       By default returns m_event.tau; for its motivation, see XaodAnalysis::xaodMuons().
     */
    virtual xAOD::TauJetContainer* xaodTaus();
    /// access the default collection of jets from the D3PDReader
    /**
       By default returns m_event.jet_AntiKt4LCTopo; for its motivation, see XaodAnalysis::xaodMuons().
       \todo In this case there might be some ambiguity to be sorted out when calling SUSYObjDef::GetMET().
     */
    virtual xAOD::JetContainer* xaodJets();
    /// access the default collection of photons
    virtual const xAOD::PhotonContainer* xaodPhothons();
    /// access the truth event
    virtual const xAOD::TruthEventContainer* xaodTruthEvent();
    /// access the truth particles
    virtual const xAOD::TruthParticleContainer* xaodTruthParticles();
    /// delete the objects created when accessing the aux info with shallow copies
    void clearShallowCopies();

    /// fill collections of baseline+signal+truth objects
    /**
           Object selection
           Selected leptons have kinematic and cleaning cuts (no overlap removal)
           Baseline leptons = selected + overlap removed
    */
    void selectObjects(SusyNtSys sys);
    void selectBaselineObjects(SusyNtSys sys);
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
    bool matchElectronTrigger(const TLorentzVector &lv, std::vector<int>* trigBools);
    void matchMuonTriggers();
    bool matchMuonTrigger(const TLorentzVector &lv, std::vector<int>* trigBools);
    void matchTauTriggers();
    bool matchTauTrigger(const TLorentzVector &lv, std::vector<int>* trigBools);



    //
    // Event cleaning
    //

    // grl
    XaodAnalysis& setGRLFile(TString fileName);
    static std::string defauldGrlFile();
    bool initGrlTool();
    bool passGRL(const xAOD::EventInfo* eventinfo); ///< good run list
    bool passTTCVeto(); ///< incomplete TTC event veto
    bool passTileErr(const xAOD::EventInfo* eventinfo); ///< Tile error
    bool passLarErr(); ///< lar error

    bool passLarHoleVeto(); ///< lar hole veto
    bool passTileHotSpot(); ///< tile hot spot
    bool passBadJet(); ///< bad jet
    bool passGoodVtx(); ///< good vertex
    bool passTileTrip(); ///< tile trip
    bool passBadMuon(); ///< bad muon veto
    bool passCosmic(); ///< cosmic veto

    void assignEventCleaningFlags(); ///< Event level cleaning cuts
    void assignObjectCleaningFlags();// Object level cleaning cuts;  t/<hese depend on sys

    // Event weighting

    /// Full event weight includes generator, xsec, pileup, and lumi weights.
    /** Default weight uses A-E lumi.
     You can supply a different integrated luminosity,
     but the the pileup weights will still correspond to A-E.
    */
    float getEventWeight(float lumi = LUMI_A_E);
    float getXsecWeight(); ///< event weight (xsec*kfac)
    float getLumiWeight(); ///< lumi weight (lumi/sumw) normalized to 4.7/fb
    void setLumi(float lumi) { m_lumi = lumi; } ///< luminosity to normalize to (in 1/pb)
    void setSumw(float sumw) { m_sumw = sumw;  } ///< sum of mc weights for sample
    void setXsec(float xsec) { m_xsec = xsec;  } ///< user cross section, overrides susy cross section
    void setErrXsec(float err) { m_errXsec = err;  } ///< user cross section uncert
    float getPileupWeight(); ///< pileup weight for full dataset: currently A-L
    float getPileupWeightUp();
    float getPileupWeightDown();
    float getPDFWeight8TeV(); ///< PDF reweighting of 7TeV -> 8TeV
    float getLepSF(const std::vector<LeptonInfo>& leptons); ///< Lepton efficiency SF
    float getBTagSF(const std::vector<int>& jets); ///< BTag efficiency SF

    // Utility methods
    void calcRandomRunLB(); ///< calculate random run/lb numbers for MC
    int getHFORDecision(); ///< HF overlap removal decision (DG obsolete?)
    uint getNumGoodVtx(); ///< Count number of good vertices
    bool matchTruthJet(int iJet); ///< Match a reco jet to a truth jet

    // Running conditions
    TString sample() { return m_sample; } ///< Sample name - used to set isMC flag
    XaodAnalysis& setSample(TString s) { m_sample = s; return *this; }
    void setAF2(bool isAF2=true) { m_isAF2 = isAF2; } ///< AF2 flag
    void setMCProduction(MCProduction prod) { m_mcProd = prod; } ///< Set MC Production flag
    void setD3PDTag(D3PDTag tag) { m_d3pdTag = tag; } ///< Set SUSY D3PD tag to know which branches are ok
    void setSys(bool sysOn){ m_sys = sysOn; }; ///< Set sys run
    void setSelectPhotons(bool doIt) { m_selectPhotons = doIt; } ///< Toggle photon selection
    void setSelectTaus(bool doIt) { m_selectTaus = doIt; } ///< Toggle tau selection and overlap removal
    void setSelectTruthObjects(bool doIt) { m_selectTruth = doIt; } ///< Set-Get truth selection
    bool getSelectTruthObjects(         ) { return m_selectTruth; }
    void setMetFlavor(std::string metFlav); ///< only STVF and STVF_JVF are available (anything else will raise an error)
    void setDoMetMuonCorrection(bool doMetMuCorr) { m_doMetMuCorr = doMetMuCorr; }
    void setDoMetFix(bool doMetFix) { m_doMetFix = doMetFix; }
    /// whether the options specified by the user are consistent with the event info
    /**
       This function should be called when the first event is being
       read. Leave it up to the user to decide whether aborting in
       case of inconsistent options (need to check that all the input
       branches from the D3PD are filled in correctly).
       In the class inheriting from XaodAnalysis, one should have:
       \code{.cpp}
       Process() {
           GetEntry(entry);
           if(!m_flagsHaveBeenChecked) {
               m_flagsAreConsistent = runningOptionsAreValid();
               m_flagsHaveBeenChecked=true;
           }
       ...
       }
       \endcode
     */
    bool runningOptionsAreValid();
    //void setUseMetMuons(bool useMetMu) { m_useMetMuons = useMetMu; }

    //
    // Event dumps
    //
    void dumpEvent();
    void dumpBaselineObjects();
    void dumpSignalObjects();
    // helpers
    static DataStream streamFromSamplename(const TString &s, bool isdata); ///< guess data stream from sample name
    static bool isDataFromSamplename(const TString &s); ///< guess from sample name whether it's data sample
    static bool isSimuFromSamplename(const TString &s); ///< guess from sample name whether it's a simulated sample


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

    /* Root::TElectronEfficiencyCorrectionTool* m_eleMediumSFTool; */

    TString                     m_grlFileName;  // grl file name
    GoodRunsListSelectionTool*   m_grl;         // good runs list

    /* Root::TPileupReweighting*   m_pileup;       // pileup reweighting */
    /* Root::TPileupReweighting*   m_pileup_up;    // pileup reweighting */
    /* Root::TPileupReweighting*   m_pileup_dn;    // pileup reweighting */

    // The SUSY CrossSectionDB has its own map for retrieving xsec info, but
    // it has a lot of entries so lookup is slow.  Save our own xsec map
    /* SUSY::CrossSectionDB*                       m_susyXsec;     // SUSY cross section database */
    /* std::map<int,SUSY::CrossSectionDB::Process> m_xsecMap;      // our own xsec map for faster lookup times */

    RecoTauMatch                m_recoTruthMatch;       // Lepton truth matching tool

    TTree* m_tree;              // Current tree
    Long64_t m_entry;           // Current entry in the current tree (not chain index!)
    int m_dbg;                  // debug level
    bool m_isMC;                // is MC flag
    bool m_flagsAreConsistent;  ///< whether the cmd-line flags are consistent with the event
    bool m_flagsHaveBeenChecked;///< whether the cmd-line have been checked

    xAOD::TEvent m_event;
    xAOD::TStore m_store;
    ST::SUSYObjDef_xAOD m_susyObj;      // SUSY object definitions
    /// electrons from the xaod
    /**
       Note: one can only access collections as const. One can make
       modifiable shallow copies, but you have to remember to delete
       them at each event. (see SUSYToolsTester.cxx)
     */
    xAOD::MuonContainer* m_xaodMuons;
    xAOD::ShallowAuxContainer* m_xaodMuonsAux; ///< muon aux info
    xAOD::ElectronContainer* m_xaodElectrons;
    xAOD::ShallowAuxContainer* m_xaodElectronsAux; ///< electron aux info
    xAOD::TauJetContainer* m_xaodTaus;
    xAOD::ShallowAuxContainer* m_xaodTausAux; ///< tau aux info
    xAOD::JetContainer* m_xaodJets;
    xAOD::ShallowAuxContainer* m_xaodJetsAux; ///< jet aux info

    /// cleanup shallow copies and aux containers
    /**
       They are created when retrieving the collections with
       XaodAnalysis::xaodMuons etc; it's the user's responsibility to
       clean things up.
    */
    XaodAnalysis& deleteShallowCopies();
};

} // susy

#endif
