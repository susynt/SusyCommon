#ifndef SusyCommon_SusyNtMaker_h
#define SusyCommon_SusyNtMaker_h


#include <iostream>

#include "TStopwatch.h"

#include "SusyCommon/SusyD3PDAna.h"
#include "SusyNtuple/SusyNtObject.h"


/*

    SusyNtMaker - a class for making SusyNt from Susy D3PDs

*/

class SusyNtMaker : public SusyD3PDAna
{

  public:

    // Constructor and destructor
    SusyNtMaker();
    virtual ~SusyNtMaker();

    // Begin is called before looping on entries
    virtual void    Begin(TTree *tree);
    // Main event loop function
    virtual Bool_t  Process(Long64_t entry);
    // Terminate is called after looping is finished
    virtual void    Terminate();

    // Event selection - loose object/event cuts for filling tree
    virtual bool    selectEvent();

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
    void fillMetVars(SusyNtSys sys = NtSys_NOM);

    // Systematic Methods
    void doSystematic();

    void saveElectronSF(SusyNtSys sys);
    void saveMuonSF(SusyNtSys sys);
    void saveJetSF(SusyNtSys sys);

    void addMissingElectron(const LeptonInfo*, SusyNtSys sys);
    void addMissingMuon(const LeptonInfo*, SusyNtSys sys);
    void addMissingJet(int index, SusyNtSys sys);

    // Systematic enum checks
    bool isElecSys(SusyNtSys s){ 
      return (s == NtSys_EES_UP || s == NtSys_EES_DN || s == NtSys_EER_UP || s == NtSys_EER_DN); 
    }
    bool isMuonSys(SusyNtSys s){ 
      return (s == NtSys_MS_UP || s == NtSys_MS_DN || s == NtSys_ID_UP || s == NtSys_ID_DN); 
    }
    bool isJetSys(SusyNtSys s){ 
      return (s == NtSys_JES_UP || s == NtSys_JES_DN || s == NtSys_JER); 
    }
    
    void addEventFlag(SusyNtSys s, int eventFlag){ 
      m_susyNt.evt()->evtFlag[s] = eventFlag;
    };
    
 protected:
    
    TFile*              m_outTreeFile;  // output tree file
    TTree*              m_outTree;      // output tree

    Susy::SusyNtObject  m_susyNt;       // SusyNt interface

    // Some object counts
    uint                n_base_ele;
    uint                n_base_muo;
    uint                n_base_jet;

    // Some event counts
    uint                n_evt_initial;
    uint                n_evt_grl;
    uint                n_evt_larErr;
    uint                n_evt_larHole;
    uint                n_evt_badJet;
    uint                n_evt_goodVtx;
    uint                n_evt_badMu;
    uint                n_evt_cosmic;
    uint                n_evt_saved;

    // histogram to save cutflow 
    TH1F*               h_cutFlow;

    // Timer
    TStopwatch          m_timer;

};

#endif
