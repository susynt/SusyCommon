// Dear emacs, this is -*- c++ -*-
#ifndef SusyCommon_SusyNtMaker_h
#define SusyCommon_SusyNtMaker_h

#include "SusyNtuple/SusyNtObject.h" // DG-2014-08-15 note to self: reversed include order breaks things (VarHandle bug?)
#include "SusyCommon/XaodAnalysis.h"
#include "SusyCommon/SystematicMapping.h"

#include "SusyNtuple/TriggerTools.h"


#include "TStopwatch.h"

#include <iostream>
#include <string>
#include <vector>

/*

    SusyNtMaker - a class for making SusyNt from Susy D3PDs

*/

namespace Root { class TElectronEfficiencyCorrectionTool; }


namespace Susy {
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

    /// fill histogram of event-level triggers that this event triggered
    virtual void fillTriggerHisto(); //dantrim trig

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
    
    /// labels of the cut stages used in the selection
    /**
       They are used to book the histograms and the counters.
     */
    static const std::vector< std::string > cutflowLabels();
    // Initialize a cutflow histo
    TH1F* makeCutFlow(const char* name, const char* title);
    TH1F* getProcCutFlow(int signalProcess);

    //
    // SusyNt Fill methods
    //

    void fillNtVars();
    void fillEventVars();
    void fillElectronVars();
    void fillMuonVars();
    void fillJetVars();
    void fillTauVars();
    void fillPhotonVars();
    void fillTruthParticleVars();
    void storeElectron(const xAOD::Electron &in, int ele_idx);
    void storeMuon(const xAOD::Muon &in, const xAOD::MuonContainer &muons);
    void storeJet(const xAOD::Jet &in);
    void storeTau(const xAOD::TauJet &in);
    void storePhoton(const xAOD::Photon &in);
    void storeTruthParticle(const xAOD::TruthParticle &in);
    void fillMetVars(SusyNtSys sys = NtSys::NOM);
    void fillTrackMetVars(SusyNtSys sys = NtSys::NOM);
    void fillTruthJetVars();
    void fillTruthMetVars();
  
    void doSystematic();

    void storeElectronKinSys(ST::SystInfo sysInfo, SusyNtSys sys);
    void storeMuonKinSys(ST::SystInfo sysInfo, SusyNtSys sys);
    void storeJetKinSys(ST::SystInfo sysInfo, SusyNtSys sys);
    void storeTauKinSys(ST::SystInfo sysInfo, SusyNtSys sys);

    // Toggle SusyNt file writing
    void setFillNt(bool fill=true) { m_fillNt = fill; }

    // Toggle filtering
    void setFilter(bool filter=true) { m_filter = filter; }

    // Set light lepton filter
    void setNLepFilter(uint nLep) { m_nLepFilter = nLep; }
    // Set light lepton + tau filter
    void setNLepTauFilter(uint nLepTau) { m_nLepTauFilter = nLepTau; }

    // Toggle saving container taus instead of selected taus
    void setSaveContTaus(bool saveContTaus=true) { m_saveContTaus = saveContTaus; }
    static bool guessWhetherIsWhSample(const TString &samplename);
    std::string timerSummary();
    std::string counterSummary() const;
//    /**
//     * \defgroup SusyNt-specific metadata. Must be specified before writing the output.
//     * @{
//     */
//    std::string         m_inputContainerName;  ///< name of the dq2 input container
//    std::string         m_outputContainerName; ///< name of the dq2 output container
//    std::string         m_productionTag;       ///< SusyNtuple production tag
//    std::string         m_productionCommand;   ///< command used to create the ntuple (with all options and flags)
//    /**@}*/
    // Check particle pdgID to see if BSM particle    
    bool isBSM(const xAOD::TruthParticle *part);

 protected:
    SusyNtMaker& initializeOuputTree();
    SusyNtMaker& saveOutputTree();
    SusyNtMaker& initializeCutflowHistograms();
    /// write to the output file our SusyNt-specific metadata
    SusyNtMaker& writeMetadata();
 private:
    //static bool isBuggyWwSherpaSample(const int &dsid); //!< see thread "Diboson MC Truth Discrepancy" atlas-phys-susy-d3pd.cern.ch, Mar2013
    //static bool hasRadiativeBquark(const vint_t *pdg, const vint_t *status);


 protected:

    TFile*              m_outTreeFile;  // output tree file
    TTree*              m_outTree;      // output tree

    Susy::SusyNtObject  m_susyNt;       // SusyNt interface

    // Control flags
    bool                m_fillNt;       // Flag to turn off Nt filling (for fast cutflow checks)
    bool                m_filter;       // Flag to turn off filtering for signal samples
    uint                m_nLepFilter;   // Number of light leptons to filter on.
    uint                m_nLepTauFilter;// Number of leptons (light+tau) to filter on.
    bool                m_saveContTaus; // Save container taus instead of selected taus

    // Some useful flags
    bool                m_isWhSample; // is WH signal sample
    int                 m_hDecay;       // higgs decay type (see WhTruthExtractor::Hdecays)
    bool                m_hasSusyProp;  // whether this event is affected by the susy propagator bug (only for c1c1)

    // Some object counts
    uint                n_pre_ele;
    uint                n_pre_muo;
    uint                n_pre_tau;
    uint                n_pre_jet;
    uint                n_base_ele;
    uint                n_base_muo;
    uint                n_base_tau;
    uint                n_base_jet;
    uint                n_sig_ele;
    uint                n_sig_muo;
    uint                n_sig_tau;
    uint                n_sig_jet;

    // We are currently changing the procedure for counting events in the cutflow!
    // We would like to have a histogram with total raw numbers, as above, but in
    // addition we would like a similar histo filled with the generator weights.
    // Finally, for samples with multiple signal processes, we want to be able to keep
    // track of the weighted number of events for each process!

    // So, we'll use a map of histos for signal processes
    TH1F*               h_rawCutFlow;           // cutflow filled always with weight=1
    TH1F*               h_genCutFlow;           // cutflow filled with generator weights
    std::map<int,TH1F*> m_procCutFlows;         // cutflows, one for each subprocess
    std::vector< size_t > m_cutstageCounters; ///< used to print the summary cutflow table

    TH1F*               h_passTrigLevel;        // histogram storing event-level fired triggers // dantrim trig

    // Timer
    TStopwatch          m_timer;
};

} // susy

#endif
