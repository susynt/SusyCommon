//--------------------------------------------------------------------------------
//   LeptonInfo.h
//   definition of LeptonInfo class for lepton interface convenience :)
//
//
//   author: Steve Farrell <sfarrell@cern.ch>
//--------------------------------------------------------------------------------

#ifndef LeptonInfo_h
#define LeptonInfo_h

#include <vector>

#include "TLorentzVector.h"

#include "D3PDReader/ElectronD3PDObject.h"
#include "D3PDReader/MuonD3PDObject.h"



/**

    LeptonInfo - a transient class used to consolidate useful information 
                 about a single lepton.

    Constructor argument:
        isElectron              flavor bool
        index                   index of element in d3pd vectors
        lv                      pointer to corrected lorentz vector
        d3pdObj                 either ElectronD3PDObject* or MuonD3PDObject*

    Accessing variables:
        Use lv() to access the corrected 4-momentum

*/

class LeptonInfo
{

  public:

    // build object by specifying the inputs
    LeptonInfo(bool isElectron, unsigned int index, TLorentzVector* lv, TObject* d3pdObj);
    //LeptonInfo(bool isElectron, unsigned int index, MultiLepD3PDAnalysis* ana);
    ~LeptonInfo();

    // Lepton flavor
    //int pdg() const;
    bool isElectron() const { return  m_isEle; }
    bool isMuon()     const { return !m_isEle; }

    // SUSYObjDef TLorentzVector representing the fully corrected 4-momentum
    TLorentzVector* lv() const { return m_lv; }

    // Index of this lepton in corresponding d3pd arrays
    unsigned int idx() const { return m_idx; }

    // For convenience, a method to access charge that works for both electrons and muons
    virtual float charge() const;

    // Access the d3pd variables for this lepton
    // Obviously only one of these will return a non-NULL pointer for a given lepton
    // If you try to grab the wrong one, it will print a warning.
    D3PDReader::ElectronD3PDObjectElement* getElectronElement() const;
    D3PDReader::MuonD3PDObjectElement*     getMuonElement() const;

    // Print info
    virtual void print() const;

    // Comparison operators for sorting leptons by pt
    inline bool operator > (const LeptonInfo & other) const
    {
        return m_lv->Pt() > other.m_lv->Pt();
    }
    inline bool operator < (const LeptonInfo & other) const
    {
        return m_lv->Pt() < other.m_lv->Pt();
    }

    // LeptonInfos are equal if they represent the same lepton
    inline bool operator == (const LeptonInfo & other) const
    {
        if(m_isEle == other.m_isEle){
            return m_idx==other.m_idx;
        }
        return false;
    }
    inline bool operator != (const LeptonInfo & other) const
    {
        return !( (*this) == other );
    }

    // Helper function to find this LeptonInfo in a collection
    bool findThisIn(const std::vector<LeptonInfo>& leps) const
    {
        for(unsigned int i=0; i<leps.size(); i++){
            if( (*this) == leps[i] ) return true;
        }
        return false;
    }

  private:

    bool m_isEle;               // is electron (otherwise a muon)
    unsigned int m_idx;         // index of this lepton in d3pd arrays
    TLorentzVector* m_lv;       // pointer to corresponding LV in SUSYObjDef
    TObject* m_obj;             // pointer to *D3PDObject for reading this lepton

};


// Helper function for building a sorted collection of LeptonInfo
std::vector<LeptonInfo> buildLeptonInfos(D3PDReader::ElectronD3PDObject* electrons, std::vector<int> & elecIndices, 
                                         D3PDReader::MuonD3PDObject* muons, std::vector<int> & muonIndices, SUSYObjDef & susyObj);



#endif
