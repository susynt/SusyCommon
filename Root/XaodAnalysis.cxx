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
    m_is_af2(false),
    m_input_chain(0),
    m_input_container_name(""),
    m_output_container_name(""),
    m_mc_type(MCType::MCInvalid),
    m_is_derivation(false),
    m_stream(Stream_Unknown),
    m_is_data15(false),
    m_is_data16(false),
    // xAOD EDM
    m_event(xAOD::TEvent::kClassAccess), // kAthenaAccess
    m_store(),
    m_nEventsProcessed(0),
    m_sumOfWeights(0),
    m_sumOfWeightsSquared(0),
    m_run_oneST(false),
    m_eleIDDefault(eleTightLLH)
{
}
/////////////////////////////////////////////////////////////////////////////
XaodAnalysis::~XaodAnalysis()
{

}
//////////////////////////////////////////////////////////////////////////////
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
    //cout << "data_or_mc_from_name  " << endl;
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
        m_is_data15 = is_data15;
        m_is_data16 = is_data16;
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
    cout << "XaodAnalysis::SlaveBegin" << endl;



    return;
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::Init(TTree* tree)
{
    cout << "XaodAnalysis::Init" << endl;
    if(dbg()>=10) cout << "XaodAnalysis::Init" << endl;

    xAOD::Init("Susy::XaodAnalysis").ignore();
    
    get_sumw(tree);

    m_is_derivation = XaodAnalysis::is_derivation_from_metadata(tree); 
    m_stream = XaodAnalysis::stream_from_input_container(input_container(), mc());

    // get the data/ area
    char *tmparea = getenv("ROOTCOREBIN");
    char *test_area = getenv("TestArea"); // Athena, for the big dogs

    if(tmparea != NULL) {
        m_data_dir = tmparea;
        m_data_dir = m_data_dir + "/data/";
    }
    else if(test_area != NULL) {
        tmparea = test_area;
        m_data_dir = tmparea;
        m_data_dir = m_data_dir + "/";
    }
    else {
        cout << "XaodAnalysis::Init    ERROR RootCore area may not be set up, exiting" << endl;
        exit(1);
    }

    ///////////////////////////////////////////////////////////
    cout << "---------------------------------------------------------------------" << endl;
    cout << "XaodAnalysis::Init    Treating sample as "
        << (mc() ?  MCType2str(mc_type()) : "data") << endl;
    cout << "---------------------------------------------------------------------" << endl;

    // initialize all local tools
    initialize_local_tools();

    // initialize our instances of SUSYTools
    initialize_SUSYTools();

    // systematics
    if(mc() && sys()) get_systematic_list();
    else {
        ST::SystInfo infodef;
        infodef.affectsKinematics = false;
        infodef.affectsWeights = false;
        infodef.affectsType = ST::Unknown;
        systInfoList.push_back(infodef);
    }


}
//////////////////////////////////////////////////////////////////////////////
Bool_t XaodAnalysis::Process(Long64_t entry)
{
    cout << "XaodAnalysis::Process" << endl;


    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::Terminate()
{
    cout << "XaodAnalysis::Terminate" << endl;


}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::get_sumw(TTree* tree)
{
    if(!mc()) return;

    // get the sumw of all of the files
    TObjArray* chain_files = chain()->GetListOfFiles();
    TIter next(chain_files);
    TChainElement* ch_file = 0;
    int file_counter = 0;
    while (( ch_file = (TChainElement*)next() )) {
        cout << "XaodAnalsysis::get_sumw    CutBookkeepr info for: " << ch_file->GetTitle() << endl;
        TFile* f = TFile::Open(ch_file->GetTitle());
        m_event.readFrom(f);
        m_event.getEntry(0);
        if(!collect_cutbooks(m_event, file_counter)) exit(1);
        f->Close();
        f->Delete();
        file_counter++;
    } // while
    cout << "------------------------------------------------------------------" << endl;
    cout << "XaodAnalysis::get_sumw    CutBookkeeper totals: " << endl;
    cout << "XaodAnalysis::get_sumw      > # events processed   : " << m_nEventsProcessed << endl;
    cout << "XaodAnalysis::get_sumw      > sum of weights       : " << m_sumOfWeights << endl;
    cout << "XaodAnalysis::get_sumw      > sum of weights^2     : " << m_sumOfWeightsSquared << endl;
    cout << "------------------------------------------------------------------" << endl;

    // revert back
    m_event.readFrom(tree);
    m_event.getEntry(0);
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::collect_cutbooks(xAOD::TEvent& event, int file_counter)
{
    bool ok = true;
    const xAOD::CutBookkeeperContainer* completeCBC = 0;
    ok = event.retrieveMetaInput(completeCBC, "CutBookkeepers");
    if(!ok) {
        cout << "XaodAnalysis::collect_cutbooks    ERROR Failed to get CBK metdata for file "
                << file_counter << " in input chain, exiting" << endl;
        exit(1);
    }

    const xAOD::CutBookkeeper* allEventsCBK = 0;
    int maxcycle = -1;
    for(auto cbk : *completeCBC) {
        if(cbk->name() == "AllExecutedEvents" && cbk->inputStream()=="StreamAOD" && cbk->cycle() > maxcycle) {
            maxcycle = cbk->cycle();
            allEventsCBK = cbk;
        } // fi
    } // cbk
    uint64_t nevents_ = -1;
    double sumw_ = -1.0;
    double sumw2_ = -1.0;
    if(allEventsCBK) {
        nevents_ = allEventsCBK->nAcceptedEvents();
        sumw_ = allEventsCBK->sumOfEventWeights();
        sumw2_ = allEventsCBK->sumOfEventWeightsSquared();
    }
    else {
        cout << "XaodAnalysis::collect_cutbooks    ERROR \"AllExecutedEvents\" branch "
                << "not found in CBK metadata for file " << file_counter << " in input chain, exiting" << endl;
        ok = false;
    }

    m_nEventsProcessed += nevents_;
    m_sumOfWeights += sumw_;
    m_sumOfWeightsSquared += sumw2_;

    if(dbg()) {
        cout << "XaodAnalysis::collect_cutbooks    > file[" << file_counter << "] # accepted events: " << m_nEventsProcessed
            << "  sumw: " << m_sumOfWeights << "  sumw2: " << m_sumOfWeightsSquared << endl;
    }

    return ok;

}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::is_derivation_from_metadata(TTree* tree)
{
    bool is_derivation = false;
    TTree* metadata = nullptr;

    if(TDirectory* treedir = get_directory_from_chain(tree)) {
        if(treedir->Get("MetaData")) {
            metadata = dynamic_cast<TTree*>(treedir->Get("MetaData"));
        }
    } 
    if(metadata) {
        metadata->LoadTree(0);
        is_derivation = !metadata->GetBranch("StreamAOD");
    }
    else {
        cout << "XaodAnalysis::is_derivation_from_metadata    WARNING Cannot get MetaData tree" << endl;
    }

    if(is_derivation) cout << "XaodAnalysis::is_derivation_from_metadata    Treating input as a derivation" << endl;

    return is_derivation;
}
//////////////////////////////////////////////////////////////////////////////
TDirectory* XaodAnalysis::get_directory_from_chain(TTree* tree)
{
    TDirectory* dir = nullptr;
    if(tree) {
        dir = tree->GetDirectory(); // file?
        if(dir) {
            if(dbg()) cout << "XaodAnalysis::get_directory_from_chain    Got the directory from the input tree: "
                        << dir->GetName() << endl;
        }
        else {
            if(dbg()) cout << "XaodAnalysis::get_directory_from_chain    Trying to get the directory from a TChain" << endl;
            if(TChain* c = dynamic_cast<TChain*>(tree)) {
                if(TChainElement* ce = static_cast<TChainElement*>(c->GetListOfFiles()->First())) {
                    TFile* first_file = TFile::Open(ce->GetTitle());
                    dir = static_cast<TDirectory*>(first_file);
                }
            }
        }
    }
    if(dbg()) {
        cout << "XaodAnalysis::get_directory_from_chain   Got directory: " << (dir ? dir->GetName() : "NULL") << endl;
    }
    return dir;
}
//////////////////////////////////////////////////////////////////////////////
DataStream XaodAnalysis::stream_from_input_container(const TString &s, bool isMC)
{
    DataStream stream = Stream_Unknown;
    if(isMC) stream = Stream_MC;
    else if(s.Contains("main", TString::kIgnoreCase) && !isMC) stream = Stream_PhysicsMain;
    else {
        cout << "XaodAnalysis::stream_from_input_container    WARNING Cannot determinie stream "
            << "from input container '" << s << "', returning " << streamName(stream) << endl;
    }
    return stream;
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::initialize_local_tools()
{

    if(dbg())
    cout << "XaodAnalysis::initialize_local_tools    Local tool initialization begin" << endl;

    if(!initialize_GRL_tool()) exit(1);

    // electron tools
    initialize_electron_tools();

    // photon tools
    initialize_photon_tools();

    // muon tools
    initialize_muon_tools();

    // tau tools
    initialize_tau_tools();

    // isolation
    initialize_isolation_tools();

    if(dbg())
    cout << "XaodAnalysis::initialize_local_tools    Local tool initialization complete" << endl;
    

}
//////////////////////////////////////////////////////////////////////////////
string XaodAnalysis::default_grl_file()
{
    string grl_file = "";
    if(data15()) {
        grl_file = "$ROOTCOREBIN/data/SusyCommon/data15_13TeV.periodAllYear_DetStatus-v79-repro20-02_DQDefects-00-02-02_PHYS_StandardGRL_All_Good_25ns.xml";
    }
    else if(data16()) {
        grl_file = "$ROOTCOREBIN/data/SusyCommon/data16_13TeV.periodAllYear_DetStatus-v88-pro20-21_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml"; 
    }
    else {
        cout << "XaodAnalysis::default_grl_file    ERROR Inconsistent data flags, cannot get GRL" << endl;
        exit(1);
    }
    return grl_file;
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::initialize_GRL_tool()
{
    if(mc()) return true;

    if(dbg()>=5) cout << "XaodAnalysis::initialize_GRL_tool" << endl;
    bool ok = true;

    m_grl_tool = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    vector<string> grl_files;
    grl_files.push_back(XaodAnalysis::default_grl_file());
    cout << "XaodAnalysis::initialize_GRL_tool    Using GRL file: " << grl_files.at(0) << endl;
    m_grl_tool->setProperty("GoodRunsListVec", grl_files);
    m_grl_tool->setProperty("PassThrough", false);
    if(!m_grl_tool->initialize().isSuccess()) ok = false;
    return ok;
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::initialize_electron_tools()
{
    if(dbg()>=5) cout << "XaodAnalysis::initialize_electron_tools" << endl;

    // Initialize the electron likelihood ID tools for each of the
    // working points

    std::string wp_veryloose = "VeryLooseLHElectron";
    std::string wp_loose = "LooseLHElectron";
    std::string wp_loose_blayer = "LooseBLLHElectron";
    std::string wp_medium = "MediumLHElectron";
    std::string wp_tight = "TightLHElectron";

    //handle
    string tool_name = "";

    // veryloose
    if(!m_elecSelLikelihoodVeryLoose.isUserConfigured()) {
        tool_name = "SUSY_EleLH_" + wp_veryloose;
        SET_DUAL_TOOL(m_elecSelLikelihoodVeryLoose, AsgElectronLikelihoodTool, tool_name);
        CHECK( m_elecSelLikelihoodVeryLoose.setProperty("WorkingPoint", wp_veryloose) );
        CHECK( m_elecSelLikelihoodVeryLoose.retrieve() );
    } // configured

    // loose
    if(!m_elecSelLikelihoodLoose.isUserConfigured()) {
        tool_name = "SUSY_EleLH_" + wp_loose;
        SET_DUAL_TOOL(m_elecSelLikelihoodLoose, AsgElectronLikelihoodTool, tool_name);
        CHECK( m_elecSelLikelihoodLoose.setProperty("WorkingPoint", wp_loose) );
        CHECK( m_elecSelLikelihoodLoose.retrieve() );
    } // configured

    // loose + b-layer
    if(!m_elecSelLikelihoodLooseBLayer.isUserConfigured()) {
        tool_name = "SUSY_EleLH_" + wp_loose_blayer;
        SET_DUAL_TOOL(m_elecSelLikelihoodLooseBLayer, AsgElectronLikelihoodTool, tool_name);
        CHECK( m_elecSelLikelihoodLooseBLayer.setProperty("WorkingPoint", wp_loose_blayer) );
        CHECK( m_elecSelLikelihoodLooseBLayer.retrieve() );
    } // configured 

    // medium
    if(!m_elecSelLikelihoodMedium.isUserConfigured()) {
        tool_name = "SUSY_EleLH_" + wp_medium;
        SET_DUAL_TOOL(m_elecSelLikelihoodMedium, AsgElectronLikelihoodTool, tool_name);
        CHECK( m_elecSelLikelihoodMedium.setProperty("WorkingPoint", wp_medium) );
        CHECK( m_elecSelLikelihoodMedium.retrieve() );
    } // configured 

    // tight
    if(!m_elecSelLikelihoodTight.isUserConfigured()) {
        tool_name = "SUSY_EleLH_" + wp_tight;
        SET_DUAL_TOOL(m_elecSelLikelihoodTight, AsgElectronLikelihoodTool, tool_name);

        CHECK( m_elecSelLikelihoodTight.setProperty("WorkingPoint", wp_tight) ); 
        CHECK( m_elecSelLikelihoodTight.retrieve() );
    } // configured

}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::initialize_photon_tools()
{
    if(dbg()>=5) cout << "XaodAnalysis::initialize_photon_tools" << endl;

    // Initialize photon selection tools
    string tool_name = "";
    
    // loose
    if(!m_photonSelLoose.isUserConfigured()) {
        tool_name = "SUSY_PhoSel_Loose";
        SET_DUAL_TOOL(m_photonSelLoose, AsgPhotonIsEMSelector, tool_name);
        CHECK( m_photonSelLoose.setProperty("WorkingPoint", "LoosePhoton") );
        CHECK( m_photonSelLoose.retrieve() );
    } // configured
    
    // tight
    if(!m_photonSelTight.isUserConfigured()) {
        tool_name = "SUSY_PhoSel_Tight";
        SET_DUAL_TOOL(m_photonSelTight, AsgPhotonIsEMSelector, tool_name);
        CHECK( m_photonSelTight.setProperty("WorkingPoint", "TightPhoton") );
        CHECK( m_photonSelTight.retrieve() );
    } // configured
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::initialize_muon_tools()
{
    if(dbg()>=5) cout << "XaodAnalysis::initialize_muon_tools" << endl;

    // Initialize muon selection tools
    string tool_name = "";
    
    // very loose
    if(!m_muonSelectionToolVeryLoose.isUserConfigured()) {
        tool_name = "SUSY_MuonSelTool_VeryLoose";
        SET_DUAL_TOOL(m_muonSelectionToolVeryLoose, CP::MuonSelectionTool, tool_name);
        CHECK( m_muonSelectionToolVeryLoose.setProperty("MaxEta", 2.7) );
        CHECK( m_muonSelectionToolVeryLoose.setProperty("MuQuality", int(xAOD::Muon::VeryLoose)) );
        CHECK( m_muonSelectionToolVeryLoose.setProperty("TrtCutOff", false) ); // SUSYTools default
        CHECK( m_muonSelectionToolVeryLoose.retrieve() );
    } // configured
    
    // loose
    if(!m_muonSelectionToolLoose.isUserConfigured()) {
        tool_name = "SUSY_MuonSelTool_Loose";
        SET_DUAL_TOOL(m_muonSelectionToolLoose, CP::MuonSelectionTool, tool_name);
        CHECK( m_muonSelectionToolLoose.setProperty("MaxEta", 2.7) );
        CHECK( m_muonSelectionToolLoose.setProperty("MuQuality", int(xAOD::Muon::Loose)) );
        CHECK( m_muonSelectionToolLoose.setProperty("TrtCutOff", false) ); // SUSYTools default
        CHECK( m_muonSelectionToolLoose.retrieve() );
    } // configured
    // medium
    if(!m_muonSelectionToolMedium.isUserConfigured()) {
        tool_name = "SUSY_MuonSelTool_Medium";
        SET_DUAL_TOOL(m_muonSelectionToolMedium, CP::MuonSelectionTool, tool_name);
        CHECK( m_muonSelectionToolMedium.setProperty("MaxEta", 2.7) );
        CHECK( m_muonSelectionToolMedium.setProperty("MuQuality", int(xAOD::Muon::Medium)) );
        CHECK( m_muonSelectionToolMedium.setProperty("TrtCutOff", false) ); // SUSYTools default
        CHECK( m_muonSelectionToolMedium.retrieve() );
    } // configured
    
    // tight
    if(!m_muonSelectionToolTight.isUserConfigured()) {
        tool_name = "SUSY_MuonSelTool_Tight";
        SET_DUAL_TOOL(m_muonSelectionToolTight, CP::MuonSelectionTool, tool_name);
        CHECK( m_muonSelectionToolTight.setProperty("MaxEta", 2.7) );
        CHECK( m_muonSelectionToolTight.setProperty("MuQuality", int(xAOD::Muon::Tight)) );
        CHECK( m_muonSelectionToolTight.setProperty("TrtCutOff", false) ); // SUSYTools default
        CHECK( m_muonSelectionToolTight.retrieve() );
    } // configured
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::initialize_tau_tools()
{
    if(dbg()>=5) cout << "XaodAnalysis::initialize_tau_tools" << endl;

    // Let's select some taus
    string tool_name = "";
    string input_file = "";
    
    // dantrim (Feb 20 2017) : Recalculate tau-lep OR -- current samples (MC p-tag p2879, 20.7.8.2) are not AODFixed
    // and doing this gets rid of errors, etc... and recaulculates internal variables
    // see: https://its.cern.ch/jira/browse/ATLASG-1087
    bool tauRecalcOLR = true;
    
    // loose
    if(!m_tauSelToolLoose.isUserConfigured()) {
        tool_name = "SUSY_TauSel_Loose";
        input_file = "SUSYTools/tau_selection_loose.conf";
        SET_DUAL_TOOL(m_tauSelToolLoose, TauAnalysisTools::TauSelectionTool, tool_name);
        CHECK( m_tauSelToolLoose.setProperty("ConfigPath", input_file) );
    
        if(tauRecalcOLR) {
            CHECK( m_tauSelToolLoose.setProperty("IgnoreAODFixCheck", true) );
            CHECK( m_tauSelToolLoose.setProperty("RecalcEleOLR", true) );
        }
    
        CHECK( m_tauSelToolLoose.retrieve() );
    } // configured
    // medium
    if(!m_tauSelToolMedium.isUserConfigured()) {
        tool_name = "SUSY_TauSel_Medium";
        input_file = "SUSYTools/tau_selection_medium.conf";
        SET_DUAL_TOOL(m_tauSelToolMedium, TauAnalysisTools::TauSelectionTool, tool_name);
        CHECK( m_tauSelToolMedium.setProperty("ConfigPath", input_file) );
    
        if(tauRecalcOLR) {
            CHECK( m_tauSelToolMedium.setProperty("IgnoreAODFixCheck", true) );
            CHECK( m_tauSelToolMedium.setProperty("RecalcEleOLR", true) );
        }
    
        CHECK( m_tauSelToolMedium.retrieve() );
    } // configured
    
    // tight
    if(!m_tauSelToolTight.isUserConfigured()) {
        tool_name = "SUSY_TauSel_Tight";
        input_file = "SUSYTools/tau_selection_tight.conf";
        SET_DUAL_TOOL(m_tauSelToolTight, TauAnalysisTools::TauSelectionTool, tool_name);
        CHECK( m_tauSelToolTight.setProperty("ConfigPath", input_file) );
    
        if(tauRecalcOLR) {
            CHECK( m_tauSelToolTight.setProperty("IgnoreAODFixCheck", true) );
            CHECK( m_tauSelToolTight.setProperty("RecalcEleOLR", true) );
        }
    
        CHECK( m_tauSelToolTight.retrieve() );
    } // configured
    
    // Truth Matching for all of the infos
    m_tauTruthMatchingTool = new TauAnalysisTools::TauTruthMatchingTool("SusyCommonTauTruthMatchingTool");
    CHECK( m_tauTruthMatchingTool->initialize() );
    m_tauTruthMatchingTool->msg().setLevel( MSG::INFO );

}

//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::initialize_isolation_tools()
{
    if(dbg()>=5) cout << "XaodAnalysis::initialize_isolation_tools" << endl;

    // see: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/IsolationSelectionTool
    string tool_name = "";
    
    // Gradient Loose WP for leptons, FixedCutTight WP for photons
    if(!m_isoToolGradientLooseTight.isUserConfigured()) {
        tool_name = "SUSY_IsoTool_GLTight";
        SET_DUAL_TOOL(m_isoToolGradientLooseTight, CP::IsolationSelectionTool, tool_name);
        CHECK( m_isoToolGradientLooseTight.setProperty("ElectronWP", "GradientLoose") );
        CHECK( m_isoToolGradientLooseTight.setProperty("MuonWP",     "GradientLoose") );
        CHECK( m_isoToolGradientLooseTight.setProperty("PhotonWP",   "FixedCutTight") );
        CHECK( m_isoToolGradientLooseTight.retrieve() );
    } // configured
    
    // Gradient WP for leptons, FixedCutTightCaloOnly WP for photons
    if(!m_isoToolGradientTightCalo.isUserConfigured()) {
        tool_name = "SUSY_IsoTool_GCalo";
        SET_DUAL_TOOL(m_isoToolGradientTightCalo, CP::IsolationSelectionTool, tool_name);
        CHECK( m_isoToolGradientTightCalo.setProperty("ElectronWP", "Gradient") );
        CHECK( m_isoToolGradientTightCalo.setProperty("MuonWP",     "Gradient") );
        CHECK( m_isoToolGradientTightCalo.setProperty("PhotonWP",   "FixedCutTightCaloOnly") );
        CHECK( m_isoToolGradientTightCalo.retrieve() );                                                              } // configured
    
    // LooseTrackOnly WP, Loose WP for photons
    if(!m_isoToolLooseTrackOnlyLoose.isUserConfigured()) {
        tool_name = "SUSY_IsoTool_LooseTrackLoose";
        SET_DUAL_TOOL(m_isoToolLooseTrackOnlyLoose, CP::IsolationSelectionTool, tool_name);
        CHECK( m_isoToolLooseTrackOnlyLoose.setProperty("ElectronWP", "LooseTrackOnly") );
        CHECK( m_isoToolLooseTrackOnlyLoose.setProperty("MuonWP",     "LooseTrackOnly") );
        CHECK( m_isoToolLooseTrackOnlyLoose.setProperty("PhotonWP",   "FixedCutLoose") );
        CHECK( m_isoToolLooseTrackOnlyLoose.retrieve() );
    } // configured
    // Loose WP for leptons, FixedCutTight WP for photons
    if(!m_isoToolLoose.isUserConfigured()) {
        tool_name = "SUSY_IsoTool_Loose";
        SET_DUAL_TOOL(m_isoToolLoose, CP::IsolationSelectionTool, tool_name);
        CHECK( m_isoToolLoose.setProperty("ElectronWP", "Loose") );
        CHECK( m_isoToolLoose.setProperty("MuonWP",     "Loose") );
        CHECK( m_isoToolLoose.setProperty("PhotonWP",   "FixedCutTight") );
        CHECK( m_isoToolLoose.retrieve() );
    }  // configured
    
    // FixedCutTightTrackOnly WP for leptons, FixedCutTight WP for photons
    if(!m_isoToolTight.isUserConfigured()) {
        tool_name = "SUSY_IsoTool_Tight";
        SET_DUAL_TOOL(m_isoToolTight, CP::IsolationSelectionTool, tool_name);
        CHECK( m_isoToolTight.setProperty("ElectronWP", "FixedCutTightTrackOnly") );
        CHECK( m_isoToolTight.setProperty("MuonWP",     "FixedCutTightTrackOnly") );
        CHECK( m_isoToolTight.setProperty("PhotonWP",   "FixedCutTight") );
        CHECK( m_isoToolTight.retrieve() );
    } // configured
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::initialize_SUSYTools()
{
    for(int susyObjId : Susy::leptonIds()) {
        string idname = SusyObjId2str(static_cast<SusyObjId>(susyObjId));
        bool isEle = isEleObj(static_cast<SusyObjId>(susyObjId));
        string name = "SUSYObjDef_xAOD_" + idname;
        m_susyObj[susyObjId] = new ST::SUSYObjDef_xAOD(name);

        cout << "-------------------------------------------------------------------" << endl;
        cout << "XaodAnalysis::initialize_SUSYTools    " << name << endl;

        string config_type = isEle ? "ele" : "muo";
        string config_name = "SusyCommon/SUSYTools_SusyNt_" + config_type + idname + ".conf";
        cout << "XaodAnalysis::initialize_SUSYTools    > Using configuration: " << config_name << endl;
        cout << "-------------------------------------------------------------------" << endl;
        CHECK( m_susyObj[susyObjId]->setProperty("ConfigFile", config_name) );

        m_susyObj[susyObjId]->msg().setLevel(dbg() ? MSG::DEBUG : MSG::WARNING);

        ST::ISUSYObjDef_xAODTool::DataSource datasource = mc() ? ST::ISUSYObjDef_xAODTool::Data :
            ( af2() ? ST::ISUSYObjDef_xAODTool::AtlfastII : ST::ISUSYObjDef_xAODTool::FullSim);

        CHECK( m_susyObj[susyObjId]->setProperty("DataSource", datasource) ); 

        // pileup reweighting
        vector<string> prw_files;
        if(mc_type() != MCType::MC15b && mc_type() != MCType::MC15c) {
            prw_files.push_back("dev/SUSYTools/merged_prw.root");
        }
        else if(mc_type() == MCType::MC15b && mc_type() != MCType::MC15c) {
            prw_files.push_back("dev/SUSYTools/merged_prw_mc15b.root");
        }
        else if(mc_type() != MCType::MC15b && mc_type() == MCType::MC15c) {
            prw_files.push_back("dev/SUSYTools/merged_prw_mc15c_latest.root");

            // add our signal
            prw_files.push_back(m_data_dir + "SusyCommon/signal_prw.root");

            // add'l MC
            prw_files.push_back(m_data_dir + "SusyCommon/additional_mc_prw.root");
        }
        else {
            cout << "XaodAnalysis::initialize_SUSYTools    "
                << "Inconsistent MC options when selecting PRW config, exiting" << endl;
            exit(1);
        }

        m_susyObj[susyObjId]->setProperty("PRWConfigFiles", prw_files);

        vector<string> lumi_calc_files;
        lumi_calc_files.push_back(m_data_dir+"SusyCommon/ilumicalc_histograms_None_297730-311481_OflLumi-13TeV-008.root");
        lumi_calc_files.push_back(m_data_dir+"SusyCommon/ilumicalc_histograms_None_276262-284484.root");
        m_susyObj[susyObjId]->setProperty("PRWLumiCalcFiles", lumi_calc_files);
        cout << endl;
        cout << "XaodAnalysis::initialize_SUSYTools    Configuring PRW tool " << endl;
        cout << "XaodAnalysis::initialize_SUSYTools     + prw config files: " << endl;
        for(auto x : prw_files)
            cout << "XaodAnalysis::initialize_SUSYTools     > " << x << endl;
        cout << "XaodAnalysis::initialize_SUSYTools     + lumi calc files: " << endl;
        for(auto x : lumi_calc_files)
            cout << "XaodAnalysis::initialize_SUSYTools     > " << x << endl;
        cout << endl;

        if(m_susyObj[susyObjId]->initialize() != StatusCode::SUCCESS) {
            cout << "XaodAnalysis::initialize_SUSYTools    Cannot initialize SUSYTools, aborting" << endl;
            abort();
        }

        cout << "-----------------------------------------------------------------" << endl;
        cout << "XaodAnalysis::initialize_SUSYTools    Initialzed SUSYTools with properties: " << endl;
        for(auto& x:m_susyObj[susyObjId]->getPropertyMgr()->getProperties()){
            if(x.second->typeName()=="string"){
                string foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                cout << " Property << " << x.first << ": " << foo << endl;
            }
            else if(x.second->typeName()=="int"){
                int foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                cout << " Property << " << x.first << ": " << foo << endl;
            }
            else if(x.second->typeName()=="float"){
                float foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                cout << " Property << " << x.first << ": " << foo << endl;
            }
            else if(x.second->typeName()=="double"){
                double foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                cout << " Property << " << x.first << ": " << foo << endl;
            }
            else if(x.second->typeName()=="bool"){
                bool foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                string value = foo ? "True" : "False";
                cout << " Property << " << x.first << ": " << value << endl;
            }
        }

        if(m_run_oneST) break;

    } // susyObjId

}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::get_systematic_list()
{
    if(dbg()>=5) cout << "XaodAnalysis::get_systematic_list" << endl;
    systInfoList = m_susyObj[m_eleIDDefault]->getSystInfoList();
}
