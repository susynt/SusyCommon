#ifndef SusyCommon_SusyD3PDSkimmer_h
#define SusyCommon_SusyD3PDSkimmer_h

#include "SusyCommon/SusyD3PDAna.h"


/*

    SusyD3PDSkimmer - a class for skimming/slimming SUSY D3PDs

*/

class SusyD3PDSkimmer : public SusyD3PDAna
{

  public:

    // Constructor and destructor
    SusyD3PDSkimmer();
    virtual ~SusyD3PDSkimmer();

    // Begin is called before looping on entries
    virtual void    Begin(TTree *tree);
    // Main event loop function
    virtual Bool_t  Process(Long64_t entry);
    // Terminate is called after looping is finished
    virtual void    Terminate();

    // The capabilities of this class will grow with time
    // Initially, it will be simple.  Skim the d3pd based on basline object selection

    /*
        This class will need to
          - initialize output tree
          - set active branches
          - select objects
          - select event
          - mark event for saving
          - save event to output tree
    */

    // Initialize the output tree
    virtual void makeOutputTree(std::string fileName="NTUP_SUSY.root"){
      m_outputFile = new TFile(fileName.c_str(), "recreate");
      m_outputTree = new TTree("susy", "susy");
    }

    // Set the meta data chain
    virtual void setMetaChain(TChain* metaChain) { m_metaChain = metaChain; }

    bool selectEvent();

  protected:

    // Outputs
    TFile*      m_outputFile;           // output file
    TTree*      m_outputTree;           // skimmed output d3pd
    TChain*     m_metaChain;            // d3pd meta data to be written out

    // Skim cuts
    int         m_nBaseLepMin;          // minimum number of baseline leptons
    int         m_nSigLepMin;           // minimum number of signal leptons

    // Event counters
    uint        n_evt_initial;          // initial events read in
    uint        n_evt_nBaseLep;
    uint        n_evt_nSigLep;

};

#endif
