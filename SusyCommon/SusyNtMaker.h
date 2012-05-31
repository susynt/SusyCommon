#ifndef SusyCommon_SusyNtMaker_h
#define SusyCommon_SusyNtMaker_h


#include <iostream>

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
    void fillMetVars(SYSTEMATIC sys = NOM);

    // Systematic Methods
    void doSystematic();

    void saveElectronSF(SYSTEMATIC sys);
    void saveMuonSF(SYSTEMATIC sys);
    void saveJetSF(SYSTEMATIC sys);

    void addMissingElectron(const LeptonInfo*, SYSTEMATIC sys);
    void addMissingMuon(const LeptonInfo*, SYSTEMATIC sys);
    void addMissingJet(int index, SYSTEMATIC sys);

    bool isElecSys(SYSTEMATIC s){ return (s == EES_UP || s == EES_DN || s == EER_UP || s == EER_DN); };
    bool isMuonSys(SYSTEMATIC s){ return (s == MS_UP || s == MS_DN || s == ID_UP || s == ID_DN); };
    bool isJetSys(SYSTEMATIC s){ return (s == JES_UP || s == JES_DN || s == JER); };
    
    void addEventFlag(SYSTEMATIC s, int eventFlag){ 
      m_susyNt.evt()->evtFlag[s] =eventFlag;
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

    // histogram to save cutflow
    TH1F*               cutflow;

};

#endif
