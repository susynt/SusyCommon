#include "SusyCommon/XaodAnalysis.h"


//xAOD
#include "xAODBase/IParticleHelpers.h" // setOriginalObjectLink
#include "xAODEgamma/EgammaxAODHelpers.h"
#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

//SusyNtuple
//#include "SusyNtuple/RecoTruthClassificiation.h" // is this needed?

//Tools
//#include "ElectronPhotonSelectorTools/AsgElectronChargeIDSelectorTool.h"

//std/stl
#include <limits>
#include <algorithm> // copy_if, transform
#include <iterator> // back_inserter
#include <numeric> // accumulate
#include <iostream>
using namespace std;

using Susy::XaodAnalysis;

#undef CHECK
#define CHECK( ARG )                                                \
    do {                                                            \
        const bool result = ARG;                                    \
        if( ! result ) {                                            \
            ::Error( "XaodAnalysis", "Failed to execute: \"%s\"",   \
                     #ARG );                                        \
            exit(-1);                                               \
        }                                                           \
    } while( false )

// useful macro to initialize asg::AnaToolHandles
#define SET_DUAL_TOOL( TOOLHANDLE, TOOLTYPE, TOOLNAME )             \
    ASG_SET_ANA_TOOL_TYPE(TOOLHANDLE, TOOLTYPE);                    \
    TOOLHANDLE.setName(TOOLNAME);                                   \

//////////////////////////////////////////////////////////////////////////////
XaodAnalysis::XaodAnalysis() :
    m_dbg(0),
    m_isMC(false),
    m_input_chain(0),
    m_input_container_name(""),
    m_output_container_name(""),
    m_mc_type(MCType::MCInvalid)
{
}
/////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::set_debug(int dbg_level)
{
    m_dbg = dbg_level;
    cout << "XaodAnalysis::set_debug    Setting debug level to " << dbg_level << endl;
}
/////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::set_chain(TChain* chain)
{
    if(!chain) {
        cout << "XaodAnalysis::set_chain   Provided TChain is null, cannot continue" << endl;
        return false;
        
    }
    m_input_chain = chain;
    if(dbg())
        cout << "XaodAnalysis::set_chain    Loading chain with "
                    << m_input_chain->GetEntries() << " entries" << endl;
    return true;
}
/////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::data_or_mc_from_name(const TString &s)
{
    cout << "data_or_mc_from_name  " << endl;
    // check MC
    bool is_mc15 = false;
    bool is_mc16 = false;
    is_mc15 = s.Contains("mc15", TString::kIgnoreCase);
    is_mc16 = s.Contains("mc16", TString::kIgnoreCase);

    if(is_mc15 && is_mc16) {
        cout << "XaodAnalysis::data_or_mc_from_name    Sample seen as mc15 AND mc16!" << endl;
        return false;
    }
    if(is_mc15 || is_mc16) {
        m_isMC = true;
    }

    // check data
    bool is_data = false;
    bool is_data15 = false;
    bool is_data16 = false;
    is_data15 = s.Contains("data15", TString::kIgnoreCase);
    is_data16 = s.Contains("data16", TString::kIgnoreCase);
    if(is_data15 && is_data16) {
        cout << "XaodAnalysis::data_or_mc_from_name    Sample seen as data15 AND data16!" << endl;
        return false;
    }
    if(is_data15 || is_data16) {
        is_data = true;
    }

    if(m_isMC && is_data) {
        cout << "XaodAnalysis::data_or_mc_from_name    Sample seen as MC AND data!" << endl;
        return false;
    }
    if(is_data) m_isMC = false;

    if(!m_isMC && !is_data) {
        cout << "XaodAnalysis::data_or_mc_from_name    Cannot determine from input name if sample is MC or data!" << endl;
        return false;
    }

    return true;
}
/////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::set_input_container(std::string name)
{
    m_input_container_name = name;

    if(dbg())
        cout << "XaodAnalysis::set_input_container    Input container name: " << name << endl;
    bool type_found = data_or_mc_from_name(name);
    if(!type_found)
        cout << "XaodAnalysis::set_input_container    Input container '" << name << "' invalid" << endl;

    return type_found;
}
/////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::set_output_container(std::string name)
{
    m_output_container_name = name;
}
/////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::set_mctype(MCType type)
{

    if(type==MCType::MCInvalid) {
        cout << "XaodAnalysis::set_mctype    Provided MCType is Invalid, cannot continue" << endl;
        return false;
    }
    m_mc_type = type;
    cout << "XaodAnalysis::set_mctype    Treating input MC sample as " << MCType2str(type) << endl;
    return true;
}
/////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::SlaveBegin(TTree* tree)
{


    return;
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::Init(TTree* tree)
{


}
//////////////////////////////////////////////////////////////////////////////
Bool_t XaodAnalysis::Process(Long64_t entry)
{


    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::Terminate()
{


}
