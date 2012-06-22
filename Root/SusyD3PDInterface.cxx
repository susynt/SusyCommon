#include "SusyCommon/SusyD3PDInterface.h"

#include "MultiLep/MuonTools.h"

using namespace std;

#define GeV 1000.

/*--------------------------------------------------------------------------------*/
// SusyD3PDContainer constructor
/*--------------------------------------------------------------------------------*/
SusyD3PDContainer::SusyD3PDContainer(const Long64_t& entry) :
        evt(entry),
        ele(entry),
        muo(entry),
        jet(entry),
        met(entry),
        trk(entry),
        vtx(entry),
        trig(entry),
        gen(entry),
        truth(entry),
        truthMu(entry),
	pho(entry)
{
}
/*--------------------------------------------------------------------------------*/
// Connect tree to the D3PD objects in container
/*--------------------------------------------------------------------------------*/
void SusyD3PDContainer::ReadFrom(TTree* tree)
{
  evt.ReadFrom(tree);
  ele.ReadFrom(tree);
  muo.ReadFrom(tree);
  jet.ReadFrom(tree);
  met.ReadFrom(tree);
  trk.ReadFrom(tree);
  vtx.ReadFrom(tree);
  trig.ReadFrom(tree);
  gen.ReadFrom(tree);
  truth.ReadFrom(tree);
  truthMu.ReadFrom(tree);
  pho.ReadFrom(tree);
}

/*--------------------------------------------------------------------------------*/
// SusyD3PDInterface Constructor
/*--------------------------------------------------------------------------------*/
SusyD3PDInterface::SusyD3PDInterface() :
        d3pd(m_entry),
        m_entry(0),
        m_dbg(0),
        m_isMC(true)
{
}
/*--------------------------------------------------------------------------------*/
// Destructor
/*--------------------------------------------------------------------------------*/
SusyD3PDInterface::~SusyD3PDInterface()
{
}

/*--------------------------------------------------------------------------------*/
// Attach tree (or is it a chain???)
/*--------------------------------------------------------------------------------*/
void SusyD3PDInterface::Init(TTree* tree)
{
  if(m_dbg) cout << "SusyD3PDInterface::Init" << endl;
  m_tree = tree;
  d3pd.ReadFrom(tree);
}

/*--------------------------------------------------------------------------------*/
// The Begin() function is called at the start of the query.
// When running with PROOF Begin() is only called on the client.
// The tree argument is deprecated (on PROOF 0 is passed).
/*--------------------------------------------------------------------------------*/
void SusyD3PDInterface::Begin(TTree* /*tree*/)
{
  if(m_dbg) cout << "SusyD3PDInterface::Begin" << endl;
}

/*--------------------------------------------------------------------------------*/
// Main process loop function - This is just an example for testing
/*--------------------------------------------------------------------------------*/
Bool_t SusyD3PDInterface::Process(Long64_t entry)
{
  // Communicate the entry number to the interface objects
  GetEntry(entry);

  if(m_dbg) cout << "____________________________________________________________" << endl;

  static Long64_t chainEntry = -1;
  chainEntry++;
  if(m_dbg || chainEntry%10000==0)
  {
    cout << "**** Processing entry " << setw(6) << chainEntry 
         << " run " << setw(6) << d3pd.evt.RunNumber()
         << " event " << setw(7) << d3pd.evt.EventNumber() << " ****" << endl;
  }

  if(m_isMC){
    if( d3pd.gen.weight()->at(0).size() < 1 ){
      cout << "Found strange MC event" << endl;
    }
  }

  if(m_dbg){
    // Loop over electrons
    cout << "Electrons" << endl;
    for(Int_t iEle=0; iEle<d3pd.ele.n(); iEle++){
      const ElectronElement* electron = & d3pd.ele[iEle];
      if(electron->pt()<10*GeV) continue;
      cout << "  " << iEle << " pt " << electron->pt()/GeV << " eta " << electron->eta() << endl;
    }
  
    // Loop over muons
    cout << "Muons" << endl;
    for(Int_t iMu=0; iMu<d3pd.muo.n(); iMu++){
      const MuonElement* muon = & d3pd.muo[iMu];
      if(muon->pt()<10*GeV) continue;
      cout << "  " << iMu << " pt " << muon->pt()/GeV << " eta " << muon->eta() << endl;

      // try to match muon to truthMuon
      //const TruthMuonElement* trueMuon = getMuonTruth( &d3pd.muo, iMu, &d3pd.truthMu );
    }

    // Loop over jets
    cout << "Jets" << endl;
    for(Int_t iJet=0; iJet<d3pd.jet.n(); iJet++){
      const JetElement* jet = & d3pd.jet[iJet];
      if(jet->pt()<20*GeV) continue;
      cout << "  " << iJet << " pt " << jet->pt()/GeV << " eta " << jet->eta() << endl;
    }
  }

  return kTRUE;
}

/*--------------------------------------------------------------------------------*/
// The Terminate() function is the last function to be called during
// a query. It always runs on the client, it can be used to present
// the results graphically or save the results to file.
/*--------------------------------------------------------------------------------*/
void SusyD3PDInterface::Terminate()
{
  if(m_dbg) cout << "SusyD3PDInterface::Terminate" << endl;
}

#undef GeV
