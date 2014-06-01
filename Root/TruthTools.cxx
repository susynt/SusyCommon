//------------------------------------------
//   TruthTools.cxx
//   Implementations of functions to deal with
//   mc truth particles for filtering, etc.
//
//   author: Steve Farrell <sfarrell@cern.ch
//-------------------------------------------

#include <iostream>
#include <cassert>

#include "TLorentzVector.h"

#include "SusyCommon/TruthTools.h"

using namespace std;

//----------------------------------
Bool_t PassMllForAlpgen(Int_t mc_channel_number, Int_t mc_n, 
                        vector<float>* mc_pt, vector<float>* mc_eta, vector<float>* mc_phi, vector<float>* mc_m,
                        vector<int>* mc_pdgId, vector<int>* mc_status, vector<int>* mc_barcode,
                        vector<vector<int> >* mc_parents, vector<vector<int> >* mc_children,
                        Bool_t DEBUG_PASSMLL_ALPGEN)
{
    // Only alpgen low mass samples can fail this check
    if(!IsAlpgenLowMass(mc_channel_number)) return true;
    float mll = MllForAlpgen(mc_n,
                             mc_pt, mc_eta, mc_phi, mc_m, 
                             mc_pdgId, mc_status, mc_barcode,
                             mc_parents, mc_children,
                             DEBUG_PASSMLL_ALPGEN);
    return (mll < 40.);
}
//----------------------------------
Bool_t PassMllForAlpgen(D3PDReader::EventInfoD3PDObject *eventInfo,
                        D3PDReader::TruthParticleD3PDObject* truthParticles, Bool_t DEBUG_PASSMLL_ALPGEN)
{
    return PassMllForAlpgen(eventInfo->mc_channel_number(), truthParticles->n(), 
                            truthParticles->pt(), truthParticles->eta(), truthParticles->phi(), truthParticles->m(),
                            truthParticles->pdgId(), truthParticles->status(), truthParticles->barcode(),
                            truthParticles->parents(), truthParticles->children(),
                            DEBUG_PASSMLL_ALPGEN);
}
//----------------------------------
float MllForAlpgen(Int_t mc_n, 
                   vector<float>* mc_pt, vector<float>* mc_eta, vector<float>* mc_phi, vector<float>* mc_m, 
                   vector<int>* mc_pdgId, vector<int>* mc_status, vector<int>* mc_barcode,
                   vector<vector<int> >* mc_parents, vector<vector<int> >* mc_children,
                   Bool_t DEBUG_PASSMLL_ALPGEN)
{
    Int_t ZpdgId=23;
    float MeV2GeV(0.001);
    for(Int_t imcpart=0; imcpart<mc_n; ++imcpart){
        if(mc_pdgId->at(imcpart) == ZpdgId){
            if(DEBUG_PASSMLL_ALPGEN){
                cout << endl;
                cout << "found a Z boson at the truth level "<< endl;
                cout << "Z status     = " << mc_status ->at(imcpart) << endl;
                cout << "Z barcode    = " << mc_barcode->at(imcpart) << endl;
                cout << "Z # children = " << mc_children->at(imcpart).size() << endl;
                cout << "Z mass       = " << mc_m->at(imcpart) << endl;
            } // if debug
            // Just take mass of the Z boson
            return mc_m->at(imcpart) * MeV2GeV;
        } // if Z
    } // end mcparticle loop
    return -1.0;
}
// Alternative function for MultiLep classes
//----------------------------------
float MllForAlpgen(D3PDReader::TruthParticleD3PDObject* truthParticles,
                   Bool_t DEBUG_PASSMLL_ALPGEN)
{
    D3PDReader::TruthParticleD3PDObject *tp = truthParticles;
    return MllForAlpgen(tp->n(), 
                        tp->pt(), tp->eta(), tp->phi(), tp->m(),
                        tp->pdgId(), tp->status(), tp->barcode(),
                        tp->parents(), tp->children(),
                        DEBUG_PASSMLL_ALPGEN);
}
//----------------------------------
bool IsAlpgenLowMass(UInt_t datasetId)
{
    return (146830 <= datasetId && datasetId <= 146835) || // AlpgenJimmy Zee Mll10to60
           (146840 <= datasetId && datasetId <= 146845) || // AlpgenJimmy Zmumu Mll10to60
           (146850 <= datasetId && datasetId <= 146855);   // AlpgenJimmy Ztautau Mll10to60
}
//----------------------------------
bool IsAlpgenPythiaZll(UInt_t datasetId)
{
    // Somewhat cleaner/faster implementation
    return (110805 <= datasetId && datasetId <= 110828) || // AlpgenPythia_P2011C_Z bb/cc
           (117650 <= datasetId && datasetId <= 117655) || // AlpgenPythia_P2011C_Zee*
           (117660 <= datasetId && datasetId <= 117665) || // AlpgenPythia_P2011C_Zmumu*
           (117670 <= datasetId && datasetId <= 117675);   // AlpgenPythia_P2011C_Ztautau*
    // Davide, 2013-03-12 : ugly implementation, should use regexp,
    // discuss with Anyes and Steve.
    // The list below is from 
    // > dq2-ls "mc12_8TeV.*.AlpgenPythia_P2011C_Z*_p1328/" | cut | sort | uniq
    /*
    return (datasetId==110805 || // AlpgenPythia_P2011C_ZeeccNp0
            datasetId==110806 || // AlpgenPythia_P2011C_ZeeccNp1
            datasetId==110807 || // AlpgenPythia_P2011C_ZeeccNp2
            datasetId==110808 || // AlpgenPythia_P2011C_ZeeccNp3
            datasetId==110809 || // AlpgenPythia_P2011C_ZmumuccNp0
            datasetId==110810 || // AlpgenPythia_P2011C_ZmumuccNp1
            datasetId==110811 || // AlpgenPythia_P2011C_ZmumuccNp2
            datasetId==110812 || // AlpgenPythia_P2011C_ZmumuccNp3
            datasetId==110813 || // AlpgenPythia_P2011C_ZtautauccNp0
            datasetId==110814 || // AlpgenPythia_P2011C_ZtautauccNp1
            datasetId==110815 || // AlpgenPythia_P2011C_ZtautauccNp2
            datasetId==110816 || // AlpgenPythia_P2011C_ZtautauccNp3
            datasetId==110817 || // AlpgenPythia_P2011C_ZeebbNp0
            datasetId==110818 || // AlpgenPythia_P2011C_ZeebbNp1
            datasetId==110819 || // AlpgenPythia_P2011C_ZeebbNp2
            datasetId==110820 || // AlpgenPythia_P2011C_ZeebbNp3
            datasetId==110821 || // AlpgenPythia_P2011C_ZmumubbNp0
            datasetId==110822 || // AlpgenPythia_P2011C_ZmumubbNp1
            datasetId==110823 || // AlpgenPythia_P2011C_ZmumubbNp2
            datasetId==110824 || // AlpgenPythia_P2011C_ZmumubbNp3
            datasetId==110825 || // AlpgenPythia_P2011C_ZtautaubbNp0
            datasetId==110826 || // AlpgenPythia_P2011C_ZtautaubbNp1
            datasetId==110827 || // AlpgenPythia_P2011C_ZtautaubbNp2
            datasetId==110828 || // AlpgenPythia_P2011C_ZtautaubbNp3
            datasetId==117650 || // AlpgenPythia_P2011C_ZeeNp0
            datasetId==117651 || // AlpgenPythia_P2011C_ZeeNp1
            datasetId==117652 || // AlpgenPythia_P2011C_ZeeNp2
            datasetId==117653 || // AlpgenPythia_P2011C_ZeeNp3
            datasetId==117654 || // AlpgenPythia_P2011C_ZeeNp4
            datasetId==117655 || // AlpgenPythia_P2011C_ZeeNp5
            datasetId==117660 || // AlpgenPythia_P2011C_ZmumuNp0
            datasetId==117661 || // AlpgenPythia_P2011C_ZmumuNp1
            datasetId==117662 || // AlpgenPythia_P2011C_ZmumuNp2
            datasetId==117663 || // AlpgenPythia_P2011C_ZmumuNp3
            datasetId==117664 || // AlpgenPythia_P2011C_ZmumuNp4
            datasetId==117665 || // AlpgenPythia_P2011C_ZmumuNp5
            datasetId==117670 || // AlpgenPythia_P2011C_ZtautauNp0
            datasetId==117671 || // AlpgenPythia_P2011C_ZtautauNp1
            datasetId==117672 || // AlpgenPythia_P2011C_ZtautauNp2
            datasetId==117673 || // AlpgenPythia_P2011C_ZtautauNp3
            datasetId==117674 || // AlpgenPythia_P2011C_ZtautauNp4
            datasetId==117675 ); // AlpgenPythia_P2011C_ZtautauNp5
    */
}
//-------------------------------------------
bool IsSherpaZll(UInt_t datasetId)
{
    // List of samples from Anyes (email 2013/03/09).
    // Check whether we need to add more, or implement with a regexp,
    // filtering what comes out of:
    // > dq2-ls "mc12_8TeV.*.Sherpa_CT10_Z*.*.NTUP_SUSY*_p1328/"
    return (datasetId==147770 || // Sherpa_CT10_Zee
            datasetId==147771 || // Sherpa_CT10_Zmumu
            datasetId==147772 || // Sherpa_CT10_Ztautau
            datasetId==128975 || // Sherpa_CT10_ZeeHeavyJets
            datasetId==128976 || // Sherpa_CT10_ZmumuHeavyJets
            datasetId==128977 || // Sherpa_CT10_ZtautauHeavyJets
            datasetId==146820 || // Sherpa_CT10_ZeeLightJets
            datasetId==146821 || // Sherpa_CT10_ZmumuLightJets
            datasetId==146822 || // Sherpa_CT10_ZtautauLightJets
            datasetId==173041 || // Sherpa_CT10_DYeeM08to15
            datasetId==173042 || // Sherpa_CT10_DYeeM15to40
            datasetId==173043 || // Sherpa_CT10_DYmumuM08to15
            datasetId==173044 || // Sherpa_CT10_DYmumuM15to40
            datasetId==173045 || // Sherpa_CT10_DYtautauM08to15
            datasetId==173046 ); // Sherpa_CT10_DYtautauM15to40
}
//-------------------------------------------
std::string vecToString(const vector<int> &vec)
{
  std::stringstream ss;
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<int>(ss,", "));
  return ss.str();
}
//-------------------------------------------
bool findIndicesLeptonsZll(const vector<int>* pdg, const vector<int>* status,
                           std::pair<size_t, size_t> &result,
                           bool verbose)
{
    bool invalidInput(!pdg || !status || pdg->size()!=status->size());
    if(invalidInput){
        if(verbose)
            cout << "findIndicesLeptonsZll: invalid inputs. "
                 << "pdg = "    << pdg <<    ", size " << (pdg ? pdg->size() : -1)
                 << "status = " << status << ", size " << (status ? status->size() : -1)
                 << endl;
        return false;
    } // end if(invalidInput)
    const int kInteresingStatus(3);
    const int kPele(+11), kAele(-11), kPmu(+13), kAmu(-13), kPtau(+15), kAtau(-15);
    for(size_t iP=0; iP<pdg->size(); ++iP){
        const int &p = pdg->at(iP);
        bool goodLepton = (status->at(iP)==kInteresingStatus
                           &&
                           (p==kPele || p==kAele || p==kPmu || p==kAmu || p==kPtau || p==kAtau));
        if(goodLepton){
            assert(iP); // this would invalidate the if's below; however, good leptons cannot have
                        // index 0 because there must be incoming particles too
            if      (!result.first)  { result.first = iP;         }
            else if (!result.second) { result.second = iP; break; }
        } // end if(goodLepton)
    } // end for(iP)
    return true;
}
//-------------------------------------------
float MllForSherpa(vector<float>* pt, vector<float>* eta, vector<float>* phi, vector<float>* m, 
                   vector<int>* pdg, vector<int>* status,
                   bool verbose)
{
    float mll(-1.0), MeV2GeV(0.001);
    bool invalidInput(!pt || !eta || !phi || !m || !pdg || !status
                      || pt->size() != eta->size()
                      || eta->size()!= phi->size()
                      || phi->size()!= m->size()
                      || m->size()  != pdg->size()
                      || pdg->size()!= status->size());
    if(invalidInput){
        if(verbose) cout<<"MllForSherpa : invalid inputs"<<endl; 
        return mll;
    } // end if(invalidInput)
    std::pair<size_t, size_t> indices2l(0, 0);
    if(!findIndicesLeptonsZll(pdg, status, indices2l, verbose)) {
        if(verbose) cout<<"MllForSherpa : cannot find lepton indices"<<endl;
        return mll;
    } // end if(!findIndicesLeptonsZll)
    size_t iL0(indices2l.first), iL1(indices2l.second);
    TLorentzVector l0,l1;
    l0.SetPtEtaPhiM(pt->at(iL0), eta->at(iL0), phi->at(iL0), m->at(iL0));
    l1.SetPtEtaPhiM(pt->at(iL1), eta->at(iL1), phi->at(iL1), m->at(iL1));
    mll = (l0+l1).M() * MeV2GeV;
    return mll;
}
//----------------------------------
const char* pdgidToString(const int &id)
{
  switch(id) {
  case kAd     :  return "/d"   ;
  case kPd     :  return "d"    ;
  case kAu     :  return "/u"   ;
  case kPu     :  return "u"    ;
  case kAs     :  return "/s"   ;
  case kPs     :  return "s"    ;
  case kAc     :  return "/c"   ;
  case kPc     :  return "c"    ;
  case kAb     :  return "/b"   ;
  case kPb     :  return "b"    ;
  case kAt     :  return "/t"   ;
  case kPt     :  return "t"    ;
  case kAele   :  return "e+"   ;
  case kPele   :  return "e-"   ;
  case kAve    :  return "ve"   ;
  case kPve    :  return "ve"   ;
  case kAmu    :  return "mu+"  ;
  case kPmu    :  return "mu-"  ;
  case kAvmu   :  return "vmu"  ;
  case kPvmu   :  return "vmu"  ;
  case kAtau   :  return "tau+" ;
  case kPtau   :  return "tau-" ;
  case kAvtau  :  return "vtau" ;
  case kPvtau  :  return "vtau" ;
  case kPg     :  return "g"    ;
  case kPgam   :  return "gamma";
  case kPz     :  return "Z"    ;
  case kPw     :  return "W-"   ;
  case kAw     :  return "W+"   ;
  case kPh     :  return "h"    ;
  default                 :  return "unkn" ;
  } // end switch(id)
}
//----------------------------------
void printEvent(const vector<int>* pdg,
				const vector<int>* status,
				const vector<vector<int> >* parents)
{
  using std::left;
  using std::right;
  size_t maxNpartToPrint=30;
  maxNpartToPrint = (pdg->size() < maxNpartToPrint
					 ?
					 pdg->size() : maxNpartToPrint);
  int colW=8;
  cout
    <<"--------------------------------"<<endl
    << left  << setw(colW)<<"i"
    << left  << setw(colW)<<"status"
    << right << setw(colW)<<"par"
    << right << setw(colW)<<"id"
    << right << setw(colW)<<"name"
    << endl
    <<"--------------------------------"<<endl;

  for(size_t iP=0; iP < maxNpartToPrint; ++iP){
    int id = pdg->at(iP);
    cout
      << left  << setw(colW)<<iP
	  << right << setw(colW)<<status->at(iP)
      << right << setw(colW)<<vecToString(parents->at(iP))
      << right << setw(colW)<<id
      << right << setw(colW)<<pdgidToString(id)
      << endl;
  } // end for(iP)
}
//----------------------------------
float MllForSherpa(D3PDReader::TruthParticleD3PDObject* truthParticles, Bool_t verbose)
{
  D3PDReader::TruthParticleD3PDObject *tp = truthParticles;
  //if(verbose) printEvent(tp->pdgId(), tp->status(), tp->parent_index());
  return MllForSherpa(tp->pt(), tp->eta(), tp->phi(), tp->m(),
					  tp->pdgId(), tp->status(),
					  verbose);
}
//----------------------------------
