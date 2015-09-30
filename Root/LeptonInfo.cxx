//------------------------------------------------------------------------------------------------
//   LeptonInfo.h
//   definition of LeptonInfo class for lepton interface convenience :)
//
//
//   author: Steve Farrell <sfarrell@cern.ch>
//------------------------------------------------------------------------------------------------

#include <iostream>
#include <iomanip>

#include "SusyCommon/LeptonInfo.h"

#define GeV 1000.

using namespace std;

//------------------------------------------------------------------------------------------------
// LeptonInfo constructor
//------------------------------------------------------------------------------------------------
LeptonInfo::LeptonInfo(bool isElectron, unsigned int index, TLorentzVector* lv, TObject* d3pdObj)
        : m_isEle(isElectron), m_idx(index), m_lv(lv), m_obj(d3pdObj)
{
    // Enforce non-NULL lorentz vector points, for now
    // this is because I don't have pointer protection in the '>' operators
    if(m_lv==NULL){
        cout << "LeptonInfo: Warning, LorentzVector pointer is null!" << endl;
    }
}

//------------------------------------------------------------------------------------------------
// LeptonInfo destructor
//------------------------------------------------------------------------------------------------
LeptonInfo::~LeptonInfo()
{}

//------------------------------------------------------------------------------------------------
// Get charge of lepton
//------------------------------------------------------------------------------------------------
float LeptonInfo::charge() const
{
    // if(m_isEle) return getElectronElement()->charge();
    // else return getMuonElement()->charge();
    return 0.0;
}

//------------------------------------------------------------------------------------------------
// Access D3PD object element
//------------------------------------------------------------------------------------------------
// DROP D3PDReader::ElectronD3PDObjectElement* LeptonInfo::getElectronElement() const
// DROP {
// DROP     //cout << "getElectronElement" << endl;
// DROP     //cout << m_obj << endl;
// DROP     //m_obj->Print();
// DROP     if(!m_isEle) cout << "LeptonInfo: Warning, requesting electron variables for a muon!" << endl;
// DROP     //return (D3PDReader::ElectronD3PDObjectElement*) m_obj;
// DROP     //return dynamic_cast<D3PDReader::ElectronD3PDObjectElement*>(m_obj);
// DROP     D3PDReader::ElectronD3PDObject* obj = dynamic_cast<D3PDReader::ElectronD3PDObject*>(m_obj);
// DROP     if(obj) return & (*obj)[m_idx];
// DROP     else return NULL;
// DROP }

//------------------------------------------------------------------------------------------------
// Access D3PD object element
//------------------------------------------------------------------------------------------------
// DROP D3PDReader::MuonD3PDObjectElement* LeptonInfo::getMuonElement() const
// DROP {
// DROP     if(m_isEle) cout << "LeptonInfo: Warning, requesting muon variables for an electron!" << endl;
// DROP     //return (D3PDReader::MuonD3PDObjectElement*) m_obj;
// DROP     //return dynamic_cast<D3PDReader::MuonD3PDObjectElement*>(m_obj);
// DROP     D3PDReader::MuonD3PDObject* obj = dynamic_cast<D3PDReader::MuonD3PDObject*>(m_obj);
// DROP     if(obj) return & (*obj)[m_idx];
// DROP     else return NULL;
// DROP }

//------------------------------------------------------------------------------------------------
// Print lepton variables
//------------------------------------------------------------------------------------------------
void LeptonInfo::print() const
{
    cout.precision(2);
    cout << std::fixed << (m_isEle? "El" : "Mu") << " : idx " << m_idx << " q " << (int)charge()
         << " pt " << setw(6) << m_lv->Pt()/GeV << " eta " << setw(5) << m_lv->Eta() << endl;
    cout.precision(6);
    cout.unsetf(std::ios_base::fixed);
}

//------------------------------------------------------------------------------------------------
// Build LeptonInfo vector (and sort it)
//------------------------------------------------------------------------------------------------
// vector<LeptonInfo> buildLeptonInfos(D3PDReader::ElectronD3PDObject* electrons, vector<int> & elecIndices,
//                                     D3PDReader::MuonD3PDObject* muons, vector<int> & muonIndices, SUSYObjDef & susyObj)
// {
//     vector<LeptonInfo> lepInfos;

//     // add the electrons
//     for(unsigned int iEle = 0; iEle < elecIndices.size(); iEle++){
//         int idx = elecIndices[iEle];
//         //D3PDReader::ElectronD3PDObjectElement* element = & (*electrons)[idx];
//         //cout << element << endl;
//         //element->Print();
//         lepInfos.push_back( LeptonInfo(true, idx, &susyObj.GetElecTLV(idx), electrons) );
//     }
//     // add the muons
//     for(unsigned int iMu = 0; iMu < muonIndices.size(); iMu++){
//         int idx = muonIndices[iMu];
//         //D3PDReader::MuonD3PDObjectElement* element = & (*muons)[idx];
//         //element->Print();
//         lepInfos.push_back( LeptonInfo(false, idx, &susyObj.GetMuonTLV(idx), muons) );
//     }

//     // sort the leptons by pt
//     std::sort(lepInfos.begin(), lepInfos.end(), std::greater<LeptonInfo>());
//     return lepInfos;
// }

