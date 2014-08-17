#ifndef SusyCommon_SusyNtMaker_h
#define SusyCommon_SusyNtMaker_h

#include "SusyNtuple/SusyNtObject.h" // DG-2014-08-15 note to self: reversed include order breaks things (VarHandle bug?)
#include "SusyCommon/XaodAnalysis.h"

#include "TStopwatch.h"

#include <iostream>
#include <string>

/*

    SusyNtMaker - a class for making SusyNt from Susy D3PDs

*/

namespace Root { class TElectronEfficiencyCorrectionTool; }

namespace susy {
class SusyNtMaker : public XaodAnalysis
{

 public:
  typedef std::vector< int > vint_t;
  public:

    // Constructor and destructor
    SusyNtMaker();
    virtual ~SusyNtMaker();

    virtual void    SlaveBegin(TTree *tree);
    // Main event loop function
    virtual Bool_t  Process(Long64_t entry);
    // Terminate is called after looping is finished
    virtual void    Terminate();

    /// whether this event should be written to the output
    /**
       This selection includes the event-level criteria and the
       object-level ones
    */
    virtual bool    selectEvent();
    /// whether this event passes the event-level criteria
    /**
       These are the criteria that only depend on flags, and not on
       objects. Also increment the counters/histos used for
       bookkeeping and normalization.
     */
    virtual bool passEventlevelSelection();
    /// whether this event passes the object-level criteria
    /**
       These are the criteria that depend on the reconstructed
       objects. Also increment the counters/histos used for
       bookkeeping and normalization.
     */
    virtual bool passObjectlevelSelection();

    // Initialize a cutflow histo
    TH1F* makeCutFlow(const char* name, const char* title);
    TH1F* getProcCutFlow(int signalProcess);

    //
    // SusyNt Fill methods
    //

    void fillNtVars();
    void fillEventVars();
    void fillLeptonVars();
    void fillElectronVars(const LeptonInfo* lepIn);
    void fillMuonVars(const LeptonInfo* lepIn);
    void fillJetVars();
    void fillJetVar(int jetIdx);
    void fillPhotonVars();
    void fillPhotonVar(int phIdx);
    void fillTauVars();
    void fillTauVar(int tauIdx);
    void fillMetVars(SusyNtSys sys = NtSys_NOM);
    void fillTruthParticleVars();
    void fillTruthJetVars();
    void fillTruthMetVars();

    // Systematic Methods
    void doSystematic();

    void saveElectronSF(SusyNtSys sys);
    void saveMuonSF(SusyNtSys sys);
    void saveJetSF(SusyNtSys sys);
    void saveTauSF(SusyNtSys sys);

    // This should be updated, we have some duplicated code which is dangerous
    void addMissingElectron(const LeptonInfo*, SusyNtSys sys);
    void addMissingMuon(const LeptonInfo*, SusyNtSys sys);
    void addMissingJet(int index, SusyNtSys sys);
    void addMissingTau(int index, SusyNtSys sys);

    // Systematic enum checks
    bool isElecSys(SusyNtSys s){
      return (s == NtSys_EES_Z_UP   || s == NtSys_EES_Z_DN ||
	      s == NtSys_EES_MAT_UP || s == NtSys_EES_MAT_DN ||
	      s == NtSys_EES_PS_UP  || s == NtSys_EES_PS_DN ||
	      s == NtSys_EES_LOW_UP || s == NtSys_EES_LOW_DN ||
	      s == NtSys_EER_UP     || s == NtSys_EER_DN);
    };
    bool isMuonSys(SusyNtSys s){
      return (s == NtSys_MS_UP || s == NtSys_MS_DN || s == NtSys_ID_UP || s == NtSys_ID_DN);
    };
    bool isJetSys(SusyNtSys s){
      return (s == NtSys_JES_UP || s == NtSys_JES_DN || s == NtSys_JER);
    };
    bool isTauSys(SusyNtSys s){
      return (s == NtSys_TES_UP || s == NtSys_TES_DN);
    }

    //void addEventFlag(SusyNtSys s, int eventFlag){
      //m_susyNt.evt()->evtFlag[s] = eventFlag;
    //};

    // Toggle SusyNt file writing
    void setFillNt(bool fill=true) { m_fillNt = fill; }

    // Toggle filtering
    void setFilter(bool filter=true) { m_filter = filter; }

    // Set light lepton filter
    void setNLepFilter(uint nLep) { m_nLepFilter = nLep; }
    // Set light lepton + tau filter
    void setNLepTauFilter(uint nLepTau) { m_nLepTauFilter = nLepTau; }

    // Toggle trigger filtering
    void setFilterTrigger(bool filter=true) { m_filterTrigger = filter; }

    // Toggle saving container taus instead of selected taus
    void setSaveContTaus(bool saveContTaus=true) { m_saveContTaus = saveContTaus; }
    static bool guessWhetherIsWhSample(const TString &samplename);
    std::string timerSummary();
    std::string counterSummary() const;
 protected:
    SusyNtMaker& initializeOuputTree();
    SusyNtMaker& saveOutputTree();
    SusyNtMaker& initializeCutflowHistograms();
 private:
    //static bool isBuggyWwSherpaSample(const int &dsid); //!< see thread "Diboson MC Truth Discrepancy" atlas-phys-susy-d3pd.cern.ch, Mar2013
    //static bool hasRadiativeBquark(const vint_t *pdg, const vint_t *status);

/*     // Function to get reco + ID efficiency scale factors. */
/*     // Copied from SUSYTools in order to provide the extra flexibility to use more than one set of electron ID. */
/*     // DG May2014: this function should really be somewhere else. Do we still need it? */
/*     void get_electron_eff_sf(float& sf, float& uncert, */
/*                              const float el_cl_eta, const float pt, */
/*                              bool recoSF, bool idSF, bool triggerSF, bool isAF2, */
/*                              Root::TElectronEfficiencyCorrectionTool* electronRecoSF, */
/*                              Root::TElectronEfficiencyCorrectionTool* electronIDSF, */
/*                              Root::TElectronEfficiencyCorrectionTool* electronTriggerSF, */
/*                              int RunNumber); */

 protected:

    TFile*              m_outTreeFile;  // output tree file
    TTree*              m_outTree;      // output tree

    Susy::SusyNtObject  m_susyNt;       // SusyNt interface

    // Control flags
    bool                m_fillNt;       // Flag to turn off Nt filling (for fast cutflow checks)
    bool                m_filter;       // Flag to turn off filtering for signal samples
    uint                m_nLepFilter;   // Number of light leptons to filter on.
    uint                m_nLepTauFilter;// Number of leptons (light+tau) to filter on.
    bool                m_filterTrigger;// Only save events that pass any of our triggers
    bool                m_saveContTaus; // Save container taus instead of selected taus

    // Some useful flags
    bool                m_isWhSample;   // is WH sample
    int                 m_hDecay;       // higgs decay type (see WhTruthExtractor::Hdecays)
    bool                m_hasSusyProp;  // whether this event is affected by the susy propagator bug (only for c1c1)

    // Some object counts
    uint                n_base_ele;
    uint                n_base_muo;
    uint                n_base_tau;
    uint                n_base_jet;
    uint                n_sig_ele;
    uint                n_sig_muo;
    uint                n_sig_tau;
    uint                n_sig_jet;

    // Some event counts
    uint                n_evt_initial;
    uint                n_evt_grl;
    uint                n_evt_ttcVeto;
    uint                n_evt_WwSherpa;
    uint                n_evt_tileTrip;
    uint                n_evt_larErr;
    uint                n_evt_tileErr;
    uint                n_evt_larHole;
    uint                n_evt_hotSpot;
    uint                n_evt_badJet;
    uint                n_evt_goodVtx;
    uint                n_evt_badMu;
    uint                n_evt_cosmic;
    uint                n_evt_susyProp;
    uint                n_evt_1Lep;
    uint                n_evt_2Lep;
    uint                n_evt_3Lep;
    uint                n_evt_saved;    // number of events save in the SusyNt

    // histogram to save cutflow
    //TH1F*               h_cutFlow;

    // We are currently changing the procedure for counting events in the cutflow!
    // We would like to have a histogram with total raw numbers, as above, but in
    // addition we would like a similar histo filled with the generator weights.
    // Finally, for samples with multiple signal processes, we want to be able to keep
    // track of the weighted number of events for each process!

    // So, we'll use a map of histos for signal processes
    TH1F*               h_rawCutFlow;           // cutflow filled always with weight=1
    TH1F*               h_genCutFlow;           // cutflow filled with generator weights
    std::map<int,TH1F*> m_procCutFlows;         // cutflows, one for each subprocess

    // Timer
    TStopwatch          m_timer;

};

} // susy

#endif
