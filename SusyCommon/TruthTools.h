//------------------------------------------
//   TruthTools.h
//   Definitions of truth particle tools
//
//   author: Steve Farrell <sfarrell@cern.ch>
//-------------------------------------------

#ifndef TruthTools_h
#define TruthTools_h

#include <vector>

#include <iostream>
#include <iomanip>
#include <string>
#include <istream>
#include <sstream>

using namespace std;

#include "TObject.h"

// DROP #include "D3PDReader/EventInfoD3PDObject.h"
// DROP #include "D3PDReader/TruthParticleD3PDObject.h"

// Function written by Cristophe and Mike for filtering the low mass alpgen DY
// in order to combine it with the Sherpa Z sample, which only has Mll>40 GeV.
//   returns true if could find Z or gamma -> l+l- with Mll<40 GeV in the MC truth
//   returns false if could find Z or gamma -> l+l- with Mll>40 GeV in the MC truth
//   returns false if could not find Z->l+l- or gamma->l+l-
Bool_t PassMllForAlpgen(Int_t mc_channel_numbers, Int_t mc_n,
                        vector<float>* mc_pt, vector<float>* mc_eta, vector<float>* mc_phi, vector<float>* mc_m,
                        vector<int>* mc_pdgId, vector<int>* mc_status, vector<int>* mc_barcode,
                        vector<vector<int> >* mc_parents, vector<vector<int> >* mc_children,
                        Bool_t DEBUG_PASSMLL_ALPGEN=false);

// DROP // Alternative function for MultiLep classes
// DROP Bool_t PassMllForAlpgen(D3PDReader::TruthParticleD3PDObject* truthParticles, Bool_t DEBUG_PASSMLL_ALPGEN=false);

//! mass of the Z boson from the Alpgen truth; -1.0 if there is no Z in the MC truth record
float MllForAlpgen(Int_t mc_n,
				   vector<float>* mc_pt,
				   vector<float>* mc_eta,
				   vector<float>* mc_phi,
				   vector<float>* mc_m,
				   vector<int>* mc_pdgId,
				   vector<int>* mc_status,
				   vector<int>* mc_barcode,
				   vector<vector<int> >* mc_parents,
				   vector<vector<int> >* mc_children,
				   Bool_t DEBUG_PASSMLL_ALPGEN=false);

// DROP //! Just a wrapper of the above for D3PDReader classes
// DROP float MllForAlpgen(D3PDReader::TruthParticleD3PDObject* truthParticles,
// DROP 				   Bool_t DEBUG_PASSMLL_ALPGEN=false);
//! mass of the Z boson from the Sherpa truth; -1.0 if there is no Z in the MC truth record
float MllForSherpa(vector<float>* mc_pt,
				   vector<float>* mc_eta,
				   vector<float>* mc_phi,
				   vector<float>* mc_m,
				   vector<int>* mc_pdgId,
				   vector<int>* mc_status,
				   bool verbose=false);

//! Find the indices of the two leptons from the Z decay
/*!
  The typical Sherpa MC record looks like this:

\verbatim
----------------------------------------
i       status       par      id    name
----------------------------------------
0             11              -1      /d
1             11               1       d
...
20            11               1       d
21            11              -2      /u
22             3  [...]       21       g
23             3  [...]       -1      /d
24             3  22, 23,      11      e-
25             3  22, 23,     -11      e+
26             3  22, 23,      -1      /d
27            11 [...]          4       c
28            11 [...]         -4      /c
...
\endverbatim

We identify the Z->ll leptons as the first two leptons that have status 3.
 */
bool findIndicesLeptonsZll(const vector<int>* mc_pdgId,
						   const vector<int>* mc_status,
						   std::pair<size_t, size_t> &result,
						   bool verbose=false);
//! Check if sample is Alpgen DY
bool IsAlpgenLowMass(UInt_t channel_number);
//! Check if sample is Alpgen Z->ll
bool IsAlpgenPythiaZll(UInt_t channel_number);
//! Check if sample is Sherpa Z->ll
bool IsSherpaZll(UInt_t datasetId);

// The code below between <begin> and <end> should be removed or put within a namespace
// DG Mar2013

// <begin>
std::string vecToString(const vector<int> &vec);
//----------------------------------
enum PdgIds{
  kPd=+1, kAd=-1,
  kPu=+2, kAu=-2,
  kPs=+3, kAs=-3,
  kPc=+4, kAc=-4,
  kPb=+5, kAb=-5,
  kPt=+6, kAt=-6,
  kPele=+11, kAele=-11,
  kPve=+12, kAve=-12,
  kPmu=+13, kAmu=-13,
  kPvmu=+14, kAvmu=-14,
  kPtau=+15, kAtau=-15,
  kPvtau=+16, kAvtau=-16,
  kPg=+21,
  kPgam=+22,
  kPz=+23,
  kPw=+24, kAw=-24,
  kPh=+25
};

const char* pdgidToString(const int &id);
void printEvent(const vector<int>* pdg,
				const vector<int>* status,
				const vector<vector<int> >* parents);
// <end>
//! Just a wrapper of the above for MultiLep classes
// DROP float MllForSherpa(D3PDReader::TruthParticleD3PDObject* truthParticles, Bool_t verbose=false);
// DROP 




#endif
