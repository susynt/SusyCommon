#ifndef SusyCommon_SusyD3PDInterface_h
#define SusyCommon_SusyD3PDInterface_h


#include <iostream>
#include <iomanip>

#include "TSelector.h"
#include "TTree.h"

#include "D3PDReader/EventInfoD3PDObject.h"
#include "D3PDReader/EventShapeD3PDObject.h"
#include "D3PDReader/ElectronD3PDObject.h"
#include "D3PDReader/MuonD3PDObject.h"
#include "D3PDReader/JetD3PDObject.h"
#include "D3PDReader/PhotonD3PDObject.h"
#include "D3PDReader/TauD3PDObject.h"
#include "D3PDReader/METD3PDObject.h"
#include "D3PDReader/MissingETTruthD3PDObject.h"
#include "D3PDReader/TrackD3PDObject.h"
#include "D3PDReader/triggerBitsD3PDObject.h"
#include "D3PDReader/EFElectronD3PDObject.h"
#include "D3PDReader/TrigMuonEFInfoD3PDObject.h"
#include "D3PDReader/TrigEFTauD3PDObject.h"
#include "D3PDReader/PrimaryVertexD3PDObject.h"
#include "D3PDReader/GenEventD3PDObject.h"
#include "D3PDReader/TruthParticleD3PDObject.h"
#include "D3PDReader/TruthMuonD3PDObject.h"
//#include "D3PDReader/TruthJetD3PDObject.h" // now in Event.h -> jet_AntiKt4Truth


#include "SusyNtuple/SusyDefs.h"


/*

    SusyD3PDContainer - A basic class for holding the D3PDObjects

*/


typedef D3PDReader::ElectronD3PDObjectElement ElectronElement;
typedef D3PDReader::MuonD3PDObjectElement MuonElement;
typedef D3PDReader::TauD3PDObjectElement TauElement;
typedef D3PDReader::JetD3PDObjectElement JetElement;
typedef D3PDReader::PhotonD3PDObjectElement PhotonElement;
typedef D3PDReader::TruthMuonD3PDObjectElement TruthMuonElement;
typedef D3PDReader::JetD3PDObjectElement TruthJetElement;

class SusyD3PDContainer
{
  public:

    // Constructor 
    SusyD3PDContainer(const Long64_t& entry);

    // Connect the objects to an input tree
    void ReadFrom( TTree* tree );

    // Connect the objects to an output tree
    void WriteTo( TTree* tree );

    // Read in all variables that we need to write out as well
    void ReadAllActive();

    D3PDReader::EventInfoD3PDObject     evt;
    D3PDReader::EventShapeD3PDObject    evtShape;
    D3PDReader::ElectronD3PDObject      ele;
    D3PDReader::MuonD3PDObject          muo;
    D3PDReader::JetD3PDObject           jet;
    D3PDReader::PhotonD3PDObject        pho;
    D3PDReader::TauD3PDObject           tau;
    D3PDReader::METD3PDObject           met;
    D3PDReader::MissingETTruthD3PDObject metTruth;
    D3PDReader::TrackD3PDObject         trk;
    D3PDReader::PrimaryVertexD3PDObject vtx;
    D3PDReader::triggerBitsD3PDObject   trig;
    D3PDReader::EFElectronD3PDObject    trigEfEl;
    D3PDReader::TrigMuonEFInfoD3PDObject trigEfMu;
    D3PDReader::TrigEFTauD3PDObject     trigEfTau;
    D3PDReader::GenEventD3PDObject      gen;
    D3PDReader::TruthParticleD3PDObject truth;
    D3PDReader::TruthMuonD3PDObject     truthMu;
    D3PDReader::JetD3PDObject           truthJet;
};



/*

    SusyD3PDInterface
    A class for reading SUSY D3PDs using the interfaces defined in the SUSY MultiLep common code package

*/

class SusyD3PDInterface : public TSelector
{

  public:

    // Constructor and destructor
    SusyD3PDInterface();
    virtual ~SusyD3PDInterface();

    //
    // TSelector methods
    //

    // Init is called every time a new TTree is attached
    virtual void    Init(TTree *tree);
    // Begin is called before looping on entries
    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree){};
    // Called at the first entry of a new file in a chain
    virtual Bool_t  Notify() { return kTRUE; }
    // Terminate is called after looping is finished
    virtual void    Terminate();
    virtual void    SlaveTerminate(){};
    // Due to ROOT's stupid design, need to specify version >= 2 or the tree will not connect automatically
    virtual Int_t   Version() const {
      return 2;
    }

    // Main event loop function
    virtual Bool_t  Process(Long64_t entry);

    // Get entry simply communicates the entry number from TSelector 
    // to this class and hence to all of the VarHandles
    virtual Int_t   GetEntry(Long64_t e, Int_t getall = 0) {
      m_entry=e;
      return kTRUE;
    }


    // Container for D3PD objects - see definition above
    SusyD3PDContainer d3pd;


    //
    // Other methods
    //

    // Debug level
    void setDebug(int dbg) { m_dbg = dbg; }
    int dbg() { return m_dbg; }

    // Is MC flag
    void setIsMC(bool isMC=true) { m_isMC = isMC; }
    bool isMC() { return m_isMC; }

    // Access tree
    TTree* getTree() { return m_tree; }

    ClassDef(SusyD3PDInterface, 1);

  protected:

    TTree* m_tree;              // Current tree

    Long64_t m_entry;           // Current entry in the current tree (not chain index!)

    int m_dbg;                  // debug level
    bool m_isMC;                // is MC flag

};

#endif
