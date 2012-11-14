#include "SusyCommon/SusyD3PDSkimmer.h"

using namespace std;

/*--------------------------------------------------------------------------------*/
// SusyD3PDSkimmer Constructor
/*--------------------------------------------------------------------------------*/
SusyD3PDSkimmer::SusyD3PDSkimmer()
{
  m_nBaseLepMin = 3;
  m_nSigLepMin = 3;

  // Initialize counters
  n_evt_initial = 0;
  n_evt_nBaseLep = 0;
  n_evt_nSigLep = 0;
}
/*--------------------------------------------------------------------------------*/
// Destructor
/*--------------------------------------------------------------------------------*/
SusyD3PDSkimmer::~SusyD3PDSkimmer()
{}

/*--------------------------------------------------------------------------------*/
// The Begin() function is called at the start of the query.
// When running with PROOF Begin() is only called on the client.
// The tree argument is deprecated (on PROOF 0 is passed).
/*--------------------------------------------------------------------------------*/
void SusyD3PDSkimmer::Begin(TTree* /*tree*/)
{
  SusyD3PDAna::Begin(0);
  if(m_dbg) cout << "SusyD3PDSkimmer::Begin" << endl;

  // Book the output tree 
  m_outputFile = new TFile("NTUP_SUSY.root", "recreate");
  m_outputTree = new TTree("susy", "susy");

  // Setup the output branches
  // First, use all interfaced branches
  d3pd.evt.SetActive();
  d3pd.ele.SetActive();
  d3pd.muo.SetActive();
  d3pd.jet.SetActive();
  d3pd.pho.SetActive();
  d3pd.tau.SetActive();
  d3pd.met.SetActive();
  d3pd.trk.SetActive();
  d3pd.vtx.SetActive();
  d3pd.trig.SetActive();
  d3pd.gen.SetActive();
  if(m_isMC) d3pd.truth.SetActive();
  if(m_isMC) d3pd.truthMu.SetActive();
  if(m_isMC) d3pd.truthJet.SetActive();
  d3pd.WriteTo(m_outputTree);

  // Try saving the meta data to the output chain
  m_metaChain->Merge(m_outputFile, 0, "keep");
}

/*--------------------------------------------------------------------------------*/
// Main process loop function - This is just an example for testing
/*--------------------------------------------------------------------------------*/
Bool_t SusyD3PDSkimmer::Process(Long64_t entry)
{
  // Communicate the entry number to the interface objects
  GetEntry(entry);

  static Long64_t chainEntry = -1;
  chainEntry++;
  if(m_dbg || chainEntry%5000==0)
  {
    cout << "**** Processing entry " << setw(6) << chainEntry
         << " run " << setw(6) << d3pd.evt.RunNumber()
         << " event " << setw(7) << d3pd.evt.EventNumber() << " ****" << endl;
  }

  // Object selection
  m_susyObj.Reset();
  clearObjects();
  selectObjects();
  buildMet();

  if(selectEvent()){
    // save event to output
    d3pd.ReadAllActive();
    m_outputTree->Fill();
  }

  return kTRUE;
}

/*--------------------------------------------------------------------------------*/
// The Terminate() function is the last function to be called during
// a query. It always runs on the client, it can be used to present
// the results graphically or save the results to file.
/*--------------------------------------------------------------------------------*/
void SusyD3PDSkimmer::Terminate()
{
  SusyD3PDAna::Terminate();
  if(m_dbg) cout << "SusyD3PDSkimmer::Terminate" << endl;

  // Event counter
  cout << "Event counter" << endl;
  cout << "  Initial  " << n_evt_initial << endl;
  cout << "  nBaseLep " << n_evt_nBaseLep << endl;
  cout << "  nSigLep  " << n_evt_nSigLep << endl;
  cout << endl;

  // Write the output tree
  m_outputFile->Write(0, TObject::kOverwrite);
  cout << "SUSY D3PD saved to " << m_outputFile->GetName() << endl;
  m_outputFile->Close();
}

/*--------------------------------------------------------------------------------*/
// Event selection
/*--------------------------------------------------------------------------------*/
bool SusyD3PDSkimmer::selectEvent()
{
  n_evt_initial++;

  // Cut on number of baseline leptons
  int nBaseLep = m_baseLeptons.size();
  if(m_nBaseLepMin > 0 && nBaseLep < m_nBaseLepMin) return false;
  n_evt_nBaseLep++;

  // Cut on number of signal leptons
  int nSigLep = m_sigLeptons.size();
  if(m_nSigLepMin > 0 && nSigLep < m_nSigLepMin) return false;
  n_evt_nSigLep++;

  return true;
}
