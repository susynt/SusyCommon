#ifndef SusyCommon_SusyMetValidation_h
#define SusyCommon_SusyMetValidation_h


#include <iostream>

#include "TStopwatch.h"

#include "SusyCommon/SusyD3PDAna.h"


/*

    SusyMetValidation - a class for validating the MET from SUSYTools
    
*/

class SusyMetValidation : public SusyD3PDAna
{
  
  public:

    // Constructor and destructor
    SusyMetValidation();
    virtual ~SusyMetValidation();

    // Begin is called before looping on entries
    virtual void    Begin(TTree *tree);
    // Main event loop function
    virtual Bool_t  Process(Long64_t entry);
    // Terminate is called after looping is finished
    virtual void    Terminate();

    // Event selection 
    bool selectEvent();

    // Check the MET
    void checkMet();


 protected:
    
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
    uint                n_evt_larErr;
    uint                n_evt_larHole;
    uint                n_evt_hotSpot;
    uint                n_evt_badJet;
    uint                n_evt_goodVtx;
    uint                n_evt_badMu;
    uint                n_evt_cosmic;
    uint                n_evt_1Lep;
    uint                n_evt_2Lep;
    uint                n_evt_3Lep;
    uint                n_evt_saved;    // number of events save in the SusyNt

    // Timer
    TStopwatch          m_timer;
};

#endif
