#ifndef SusyCommon_SusyMetValidation_h
#define SusyCommon_SusyMetValidation_h


#include <iostream>

#include "TStopwatch.h"

#include "TriggerMatch/TriggerMatchMultiLep.h"
#include "SusyCommon/SusyD3PDSkimmer.h"


/*

    SusyMetValidation - a class for validating the MET from SUSYTools
    
*/

//class SusyMetValidation : public SusyD3PDAna
class SusyMetValidation : public SusyD3PDSkimmer
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

    // Trigger cut
    bool passTrigger();

    // Histograms
    void bookHistos();
    void saveHistos();

    // I'm still playing around with stuff, so things are changing a bit

    // MET flavor enum
    enum MetFlav {
      Met_unknown = -1,
      Met_susy_stvf = 0,
      Met_d3pd_stvf,
      Met_d3pd_reff,
      Met_N
    };
    std::string         metNames[Met_N];

    // Set the output histo file name
    void setHistFileName(std::string s) { m_histFileName = s; }


 protected:

    TLorentzVector      m_d3pdMet;              // MET directly from the D3PD

    // Trigger matching tool
    TriggerMatchMultiLep* m_triggerMatch;
    
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
    uint                n_evt_hfor;
    uint                n_evt_larErr;
    uint                n_evt_larHole;
    uint                n_evt_hotSpot;
    uint                n_evt_badJet;
    uint                n_evt_goodVtx;
    uint                n_evt_badMu;
    uint                n_evt_cosmic;
    uint                n_evt_3Lep;
    uint                n_evt_trig;
    uint                n_evt_sfos;
    uint                n_evt_zMass;
    uint                n_evt_met;
    uint                n_evt_bJet;
    uint                n_evt_mt;
    uint                n_evt_lepPt;
    uint                n_evt_2mu;

    float               n_evt_pileup;
    float               n_evt_lepSF;
    float               n_evt_bTagSF;
    float               n_evt_lumi;

    // Timer
    TStopwatch          m_timer;

    // Histograms
    // Do I need to split up histograms in any way?
    // I usually split by lepton channel and systematic
    // At the moment I only have the mmm channel, so maybe I can skip that for now

    std::string         m_histFileName;         // output histo file
    TFile*              m_histFile;             // output histo file

    // Met histograms
    TH1F*               h_met[Met_N];           // fully corrected met
    TH1F*               h_metEle[Met_N];        // electron term
    TH1F*               h_metMuo[Met_N];        // muon term
    TH1F*               h_metJet[Met_N];        // jet term
    TH1F*               h_metCell[Met_N];       // cell out term

    // Lepton histos
    TH1F*               h_nMetLep[Met_N];       // number of leptons in met
    TH1F*               h_nMetEl[Met_N];        // number of electrons in met
    TH1F*               h_nMetMu[Met_N];        // number of muons in each met

    TH1F*               h_metLepPt[Met_N];      // met lepton pt
    TH1F*               h_metMuPt[Met_N];       // met muon pt
    TH1F*               h_metElPt[Met_N];       // met electron pt

    TH1F*               h_metLepEta[Met_N];     // met lepton eta
    TH1F*               h_metMuEta[Met_N];      // met muon eta
    TH1F*               h_metElEta[Met_N];      // met electron eta

    // Met weights
    TH1F*               h_muWet[Met_N];         // Et met weight
    TH1F*               h_muWpx[Met_N];         // Px met weight
    TH1F*               h_muWpy[Met_N];         // Py met weight

    // Other histograms
    //TH1F*               h_nMuRaw;               // number of raw d3pd muons
    TH1F*               h_nMu;                  // Number of signal muons (same as lepton channel)

    TH1F*               h_rawMuPt;              // pt of raw muons
    TH1F*               h_rawMuEta;             // eta of raw muons

};

#endif
