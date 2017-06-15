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
    m_production_command(""),
    m_isMC(false),
    m_is_af2(false),
    m_write_ntuple(true),
    m_nt_tag(""),
    m_output_filename(""),
    m_sys(false),
    m_derivation("Unknown"),
    m_nlep_filter(1),
    m_filter_trig(false),
    m_filter(false),
    m_saveContTaus(false),
    m_store_truth(false),
    m_input_chain(0),
    m_input_container_name(""),
    m_output_container_name(""),
    m_mc_type(MCType::MCInvalid),
    m_is_derivation(false),
    m_stream(Stream_Unknown),
    m_is_data15(false),
    m_is_data16(false),
    m_nEventsProcessed(0),
    m_sumOfWeights(0),
    m_sumOfWeightsSquared(0),
    // tools
    m_grl_filename(""),
    m_grl_tool(0),
    m_elecSelLikelihoodVeryLoose(""),
    m_elecSelLikelihoodLoose(""),
    m_elecSelLikelihoodLooseBLayer(""),
    m_elecSelLikelihoodMedium(""),
    m_elecSelLikelihoodTight(""),
    m_electronChargeIDTool(""),
    m_photonSelLoose(""),
    m_photonSelTight(""),
    m_muonSelectionToolVeryLoose(""),
    m_muonSelectionToolLoose(""),
    m_muonSelectionToolMedium(""),
    m_muonSelectionToolTight(""),
    m_tauSelToolLoose(""),
    m_tauSelToolMedium(""),
    m_tauSelToolTight(""),
    m_tauTruthMatchingTool(0),
    m_isoToolGradientLooseTight(""),
    m_isoToolGradientTightCalo(""),
    m_isoToolLooseTrackOnlyLoose(""),
    m_isoToolLoose(""),
    m_isoToolTight(""),
    m_run_oneST(false),
    m_eleIDDefault(eleTightLLH),
    // xAOD EDM
    m_event(xAOD::TEvent::kClassAccess), // kAthenaAccess
    m_store(),
    // xAOD containers
    m_xaodEventInfo(nullptr),
    m_evtTrigBits(m_nTriggerBits)
{
    m_triggerNames.clear();
    clear_output_objects();
    clear_containers();
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
void XaodAnalysis::set_nlep_filter(int nlep)
{
    m_nlep_filter = nlep;
    if(nlep > 0) m_filter = true;
}
/////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::set_trig_filter(bool doit)
{
    m_filter_trig = doit;
    if(doit) m_filter = true;
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
    return;
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::Init(TTree* tree)
{
    if(dbg()) cout << "XaodAnalysis::Init" << endl;

    xAOD::Init("Susy::XaodAnalysis").ignore();
    
    get_sumw(tree);

    m_is_derivation = is_derivation_from_metadata(tree); 
    m_derivation = get_derivation_type(m_event);
    m_stream = stream_from_input_container(input_container(), mc());

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

    // initialize our instances of SUSYTools
    initialize_SUSYTools();

    // initialize all local tools
    initialize_local_tools();


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
    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::Terminate()
{

    // clean up
    delete m_tauTruthMatchingTool;

    for(int i : Susy::leptonIds()) {
        delete m_susyObj[i];
        if(m_run_oneST) break;
    }
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
        cout << "XaodAnalsysis::get_sumw    CutBookkeeper info for: " << ch_file->GetTitle() << endl;
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
TString XaodAnalysis::get_derivation_type(xAOD::TEvent& event)
{
    TString derivation = "Unknown";
    bool ok = true;
    const xAOD::CutBookkeeperContainer* completeCBC = 0;
    ok = event.retrieveMetaInput(completeCBC, "CutBookkeepers");
    if(!ok) {
        cout << "XaodAnalysis::get_derivation_type    Failed to retrieve CBK metadata, exiting" << endl;
        exit(1);
    }
    for(auto cbk : *completeCBC) {
        if(cbk->name() == "AllExecutedEvents" && TString(cbk->inputStream()).Contains("StreamDAOD")) {
            derivation = TString(cbk->inputStream()).ReplaceAll("Stream","");
            cout << "XaodAnalysis::get_derivation_type    Identified DAOD flavor = " << derivation << endl;
        }
    }
    return derivation;
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
    initialize_chargeflip_tagger();

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
void XaodAnalysis::initialize_chargeflip_tagger()
{
    if(!m_electronChargeIDTool.isUserConfigured()) {
        std::string tool_name = "ElectronChargeIDTool_medium";
        SET_DUAL_TOOL(m_electronChargeIDTool, AsgElectronChargeIDSelectorTool, tool_name);
    
        //default cut value for https://twiki.cern.ch/twiki/bin/view/AtlasProtected/ElectronChargeFlipTaggerTool
        float BDTcut = -0.28087; // medium 97%
        CHECK( m_electronChargeIDTool.setProperty("TrainingFile", "ElectronPhotonSelectorTools/ChargeID/ECIDS_20161125for2017Moriond.root"));
        CHECK( m_electronChargeIDTool.setProperty("CutOnBDT", BDTcut));
        CHECK( m_electronChargeIDTool.retrieve() );
    }
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
                cout << "XaodAnalysis::initialize_SUSYTools     - ["<<idname<<"] Property << " << x.first << ": " << foo << endl;
            }
            else if(x.second->typeName()=="int"){
                int foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                cout << "XaodAnalysis::initialize_SUSYTools     - ["<<idname<<"] Property << " << x.first << ": " << foo << endl;
            }
            else if(x.second->typeName()=="float"){
                float foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                cout << "XaodAnalysis::initialize_SUSYTools     - ["<<idname<<"] Property << " << x.first << ": " << foo << endl;
            }
            else if(x.second->typeName()=="double"){
                double foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                cout << "XaodAnalysis::initialize_SUSYTools     - ["<<idname<<"] Property << " << x.first << ": " << foo << endl;
            }
            else if(x.second->typeName()=="bool"){
                bool foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                string value = foo ? "True" : "False";
                cout << "XaodAnalysis::initialize_SUSYTools     - ["<<idname<<"] Property << " << x.first << ": " << value << endl;
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
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::fill_event_cleaning_flags()
{
    const xAOD::EventInfo* eventinfo = xaodEventInfo();
    if(passGRL(eventinfo)) m_cutFlags |= ECut_GRL;
    if(passTTCVeto(eventinfo)) m_cutFlags |= ECut_TTC;
    if(passLarErr(eventinfo)) m_cutFlags |= ECut_LarErr;
    if(passTileErr(eventinfo)) m_cutFlags |= ECut_TileErr;
    if(passSCTErr(eventinfo)) m_cutFlags |= ECut_SCTErr;
    if(passGoodVtx()) m_cutFlags |= ECut_GoodVtx;
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::passGRL(const xAOD::EventInfo* ei)
{
    return (mc() || m_grl_tool->passRunLB(ei->runNumber(), ei->lumiBlock()));
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::passTTCVeto(const xAOD::EventInfo* ei)
{
    bool eventPassesTTC = ei->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ? false : true;
    return eventPassesTTC;
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::passLarErr(const xAOD::EventInfo* ei)
{
    bool eventPassesLarErr = ei->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ? false : true;
    return eventPassesLarErr;
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::passTileErr(const xAOD::EventInfo* ei)
{
     bool eventPassesTileTrip = ei->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ? false : true;
    return eventPassesTileTrip;
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::passSCTErr(const xAOD::EventInfo* ei)
{
    bool passSCTerr = ei->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error ? false : true;
    return passSCTerr;
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::passGoodVtx()
{
    if(m_susyObj[m_eleIDDefault]->GetPrimVtx()==nullptr) return false;
    return true;
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::fill_object_cleaning_flags()
{
    // see: https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/SusyObjectDefinitionsr2013TeV
    // flags are based on nominal objects
    SusyNtSys sys=NtSys::NOM;
    ST::SystInfo sysInfo =  systInfoList[0]; 
    if(passBadJet(sysInfo, sys)) m_cutFlags |= ECut_BadJet;
    if(passBadMuon(sysInfo, sys)) m_cutFlags |= ECut_BadMuon;
    if(passCosmic(sysInfo, sys)) m_cutFlags |= ECut_Cosmic;
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::passBadJet(ST::SystInfo sysInfo, SusyNtSys sys)
{
    xAOD::JetContainer* jets = xaodJets(sysInfo, sys);
    bool pass_cleaning = true;
    for(auto &i : m_preJets) {
        bool is_bad = (bool)jets->at(i)->auxdata<char>("bad")==1;
        if(is_bad) pass_cleaning = false;
    }
    return pass_cleaning;
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::passBadMuon(ST::SystInfo sysInfo, SusyNtSys sys)
{
    xAOD::MuonContainer* muons = xaodMuons(sysInfo, sys);
    bool pass = true;
    for(auto &i : m_preMuons) {
        bool pass_baseline = (bool)muons->at(i)->auxdata<char>("baseline")==1;
        bool is_bad = (bool)muons->at(i)->auxdata<char>("bad")==1;
        if(pass_baseline && is_bad) pass = false;
    }
    return pass;
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::passCosmic(ST::SystInfo sysInfo, SusyNtSys sys)
{
    xAOD::MuonContainer* muons = xaodMuons(sysInfo, sys);
    bool pass = true;
    for(auto &i : m_preMuons) {
        bool pass_OR = (bool)muons->at(i)->auxdata<char>("passOR")==1;
        bool pass_baseline = (bool)muons->at(i)->auxdata<char>("baseline")==1;
        bool is_cosmic = (bool)muons->at(i)->auxdata<char>("cosmic")==1;
        if(pass_OR && pass_baseline && is_cosmic) pass = false;
    }
    return pass;
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::clear_output_objects(bool do_nominal)
{
    m_preElectrons.clear();
    m_preMuons.clear();
    m_preLeptons.clear();
    m_preJets.clear();
    m_preTaus.clear();
    m_contTaus.clear();
    m_prePhotons.clear(); 

    m_baseElectrons.clear();
    m_baseMuons.clear();
    m_baseLeptons.clear();
    m_baseJets.clear();
    m_baseTaus.clear();
    m_basePhotons.clear();

    m_sigElectrons.clear();
    m_sigMuons.clear();
    m_sigLeptons.clear();
    m_sigJets.clear();
    m_sigTaus.clear();
    m_sigPhotons.clear();

    m_cutFlags = 0;

    if(do_nominal) {
        m_preElectrons_nom.clear();
        m_preMuons_nom.clear();
        m_preLeptons_nom.clear();
        m_preJets_nom.clear();
        m_preTaus_nom.clear();
        m_contTaus_nom.clear();
        m_prePhotons_nom.clear();
    }
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::delete_shallow_copies(bool do_nominal)
{
    if(dbg()>=5) cout << "XaodAnalysis::delete_shallow_copies   (do nominal: " << do_nominal << ")" << endl;

    if(m_metContainer) delete m_metContainer;
    if(m_metAuxContainer) delete m_metAuxContainer;
    if(m_trackMetContainer) delete m_trackMetContainer;
    if(m_trackMetAuxContainer) delete m_trackMetAuxContainer;

    if(do_nominal)
        m_store.clear();

    clear_containers(do_nominal);
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::clear_containers(bool do_nominal)
{

    if(dbg()>=5) cout << "XaodAnalysis::clear_containers    (do nominal: " << do_nominal << ")" << endl;

    // electrons
    m_xaodElectrons = 0;
    m_xaodElectronsAux = 0;

    // muons
    m_xaodMuons = 0;
    m_xaodMuonsAux = 0;

    // jets
    m_xaodJets = 0;
    m_xaodJetsAux = 0;

    // taus
    m_xaodTaus = 0;
    m_xaodTausAux = 0;

    // photons
    m_xaodPhotons = 0;
    m_xaodPhotonsAux = 0;

    // met
    m_metContainer = 0;
    m_metAuxContainer = 0;

    // track met
    m_trackMetContainer = 0;
    m_trackMetAuxContainer = 0;

    if(do_nominal) {

        // electrons
        m_xaodElectrons_nom = 0;
        m_xaodElectronsAux_nom = 0;

        // muons
        m_xaodMuons_nom = 0;
        m_xaodMuonsAux_nom = 0;

        // jets
        m_xaodJets_nom = 0;
        m_xaodJetsAux_nom = 0;

        // taus
        m_xaodTaus_nom = 0;
        m_xaodTausAux_nom = 0;

        // photons
        m_xaodPhotons_nom = 0;
        m_xaodPhotonsAux_nom = 0;

        m_xaodTruthEvent = 0;
        m_xaodTruthParticles = 0;
        m_xaodTruthParticlesAux = 0;

        m_xaodEventInfo = 0;
        m_xaodVertices = 0;

    }
    
}
//////////////////////////////////////////////////////////////////////////////
vector<string> XaodAnalysis::xaodTriggers()
{
    if(m_triggerNames.size()==0) {
        m_triggerNames = getTrigNames("run2");
        return m_triggerNames;
    }
    else {
        return m_triggerNames;
    }
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::retrieve_xaod_collections()
{

    if(dbg()) cout << "XaodAnalysis::retrieve_xaod_collections" << endl;

    // EventInfo object
    xaodEventInfo();

    // reconstructed primary vertices
    xaodVertices();

    // electrons
    xaodElectrons(systInfoList[0]);

    // muons
    xaodMuons(systInfoList[0]);

    // jets
    xaodJets(systInfoList[0]);

    // taus
    xaodTaus(systInfoList[0]);

    // photons
    xaodPhotons(systInfoList[0]);

    // met
    retrieveXaodMet(systInfoList[0]);

    // track met
    retrieveXaodTrackMet(systInfoList[0]);

    // truth particles
    xaodTruthParticles();



}
//////////////////////////////////////////////////////////////////////////////
const xAOD::EventInfo* XaodAnalysis::retrieveEventInfo(xAOD::TEvent &e, bool dbg)
{
    const xAOD::EventInfo* ei = nullptr;
    e.retrieve(ei, "EventInfo");
    if(dbg) {
        if(ei) cout << "XaodAnalysis::retrieveEventInfo    EventInfo retrieved " << ei << endl;
        else   cout << "XaodAnalysis::retrieveEventInfo    WARNING EventInfo unable to be retrieved" << endl;
    }
    return ei;
}
//////////////////////////////////////////////////////////////////////////////
const xAOD::EventInfo* XaodAnalysis::xaodEventInfo()
{
    if(!m_xaodEventInfo) {
        m_xaodEventInfo = retrieveEventInfo(m_event, dbg());
    }
    return m_xaodEventInfo;
}
//////////////////////////////////////////////////////////////////////////////
const xAOD::VertexContainer* XaodAnalysis::retrieveVertices(xAOD::TEvent &e, bool dbg)
{
    const xAOD::VertexContainer* vtx = nullptr;
    e.retrieve(vtx, "PrimaryVertices");
    if(dbg) {
        if(vtx) cout << "XaodAnalysis::retrieveVertices    Vertex container retrieved " << vtx << "  (size: " << vtx->size() << ")" << endl;
        else    cout << "XaodAnalysis::retrieveVertices    WARNING Vertex container unable to be retrieved" << endl;
    }
    return vtx;
}
//////////////////////////////////////////////////////////////////////////////
const xAOD::VertexContainer* XaodAnalysis::xaodVertices()
{
    if(!m_xaodVertices) {
        m_xaodVertices = retrieveVertices(m_event, dbg());
    }
    return m_xaodVertices;
}
//////////////////////////////////////////////////////////////////////////////
xAOD::ElectronContainer* XaodAnalysis::xaodElectrons(ST::SystInfo sysInfo, SusyNtSys sys)
{
    bool syst_affects = ST::testAffectsObject(xAOD::Type::Electron, sysInfo.affectsType);
    if(sys!=NtSys::NOM && syst_affects) {
        if(!m_xaodElectrons) {
            CHECK( m_susyObj[m_eleIDDefault]->GetElectrons(m_xaodElectrons, m_xaodElectronsAux, true) );
            if(dbg()>=5) cout << "XaodAnalysis::xaodElectrons    Electrons (syst affected) retrieved (size: " << m_xaodElectrons->size() << ")" << endl;
        }
        return m_xaodElectrons;
    }
    else {
        if(!m_xaodElectrons_nom) {
            CHECK( m_susyObj[m_eleIDDefault]->GetElectrons(m_xaodElectrons_nom, m_xaodElectronsAux_nom, true) );
            if(dbg()) cout << "XaodAnalysis::xaodElectrons    Electrons (nominal) retrieved (size: " << m_xaodElectrons_nom->size() << ")" << endl;
        }
        return m_xaodElectrons_nom;
    }
    return nullptr;
}
//////////////////////////////////////////////////////////////////////////////
xAOD::MuonContainer* XaodAnalysis::xaodMuons(ST::SystInfo sysInfo, SusyNtSys sys)
{
    bool syst_affects = ST::testAffectsObject(xAOD::Type::Muon, sysInfo.affectsType);
    if(sys!=NtSys::NOM && syst_affects) {
        if(!m_xaodMuons) {
            CHECK( m_susyObj[m_eleIDDefault]->GetMuons(m_xaodMuons, m_xaodMuonsAux, true) );
            if(dbg()>=5) cout << "XaodAnalysis::xaodMuons    Muons (syst affected) retrieved (size: " << m_xaodMuons->size() << ")" << endl;
        }
        return m_xaodMuons;
    }
    else {
        if(!m_xaodMuons_nom) {
            CHECK( m_susyObj[m_eleIDDefault]->GetMuons(m_xaodMuons_nom, m_xaodMuonsAux_nom, true) );
            if(dbg()>=5) cout << "XaodAnalysis::xaodMuons    Muons (nominal) retrieved (size: " << m_xaodMuons_nom->size() << ")" << endl;
        }
        return m_xaodMuons_nom;
    }
    return nullptr;
}
//////////////////////////////////////////////////////////////////////////////
xAOD::JetContainer* XaodAnalysis::xaodJets(ST::SystInfo sysInfo, SusyNtSys sys)
{
    bool syst_affects = ST::testAffectsObject(xAOD::Type::Jet, sysInfo.affectsType);
    if(sys!=NtSys::NOM && syst_affects) {
        if(!m_xaodJets) {
            CHECK( m_susyObj[m_eleIDDefault]->GetJets(m_xaodJets, m_xaodJetsAux, true) );
            if(dbg()>=5) cout << "XaodAnalysis::xaodJets    Jets (syst affected) retrieved (size: " << m_xaodJets->size() << ")" << endl;
        }
        return m_xaodJets;
    }
    else {
        if(!m_xaodJets_nom) {
            CHECK( m_susyObj[m_eleIDDefault]->GetJets(m_xaodJets_nom, m_xaodJetsAux_nom, true) );
            if(dbg()>=5) cout << "XaodAnalysis::xaodJets    Jets (nominal) retrieved (size: " << m_xaodJets_nom->size() << ")" << endl;
        }
        return m_xaodJets_nom;
    }
    return nullptr;
}
//////////////////////////////////////////////////////////////////////////////
xAOD::TauJetContainer* XaodAnalysis::xaodTaus(ST::SystInfo sysInfo, SusyNtSys sys)
{
    bool syst_affects = ST::testAffectsObject(xAOD::Type::Tau, sysInfo.affectsType);
    if(sys!=NtSys::NOM && syst_affects) {
        if(!m_xaodTaus) {
            CHECK( m_susyObj[m_eleIDDefault]->GetTaus(m_xaodTaus, m_xaodTausAux, true) );
            if(dbg()>=5) cout << "XaodAnalysis::xaodTaus    Taus (syst affected) retrieved (size: " << m_xaodTaus->size() << ")" << endl;
        }
        return m_xaodTaus;
    }
    else {
        if(!m_xaodTaus_nom) {
            CHECK( m_susyObj[m_eleIDDefault]->GetTaus(m_xaodTaus_nom, m_xaodTausAux_nom, true) );
            if(dbg()>=5) cout << "XaodAnalysis::xaodTaus    Taus (nominal) retrieved (size: " << m_xaodTaus_nom->size() << ")" << endl;
        }
        return m_xaodTaus_nom;
    }
    return nullptr;
}
//////////////////////////////////////////////////////////////////////////////
xAOD::PhotonContainer* XaodAnalysis::xaodPhotons(ST::SystInfo sysInfo, SusyNtSys sys)
{
    bool syst_affects = ST::testAffectsObject(xAOD::Type::Photon, sysInfo.affectsType);
    if(sys!=NtSys::NOM && syst_affects) {
        if(!m_xaodPhotons) {
            CHECK( m_susyObj[m_eleIDDefault]->GetPhotons(m_xaodPhotons, m_xaodPhotonsAux, true) );
            if(dbg()>=5) cout << "XaodAnalysis::xaodPhotons    Photons (syst affected) retrieved (size: " << m_xaodPhotons->size() << ")" << endl;
        }
        return m_xaodPhotons;
    }
    else {
        if(!m_xaodPhotons_nom) {
            CHECK( m_susyObj[m_eleIDDefault]->GetPhotons(m_xaodPhotons_nom, m_xaodPhotonsAux_nom, true) );
            if(dbg()>=5) cout << "XaodAnalysis::xaodPhotons    Photons (nominal) retrieved (size: " << m_xaodPhotons_nom->size() << ")" << endl;
        }
        return m_xaodPhotons_nom;
    }
    return nullptr;
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::retrieveXaodMet(ST::SystInfo sysInfo, SusyNtSys sys)
{
    if(dbg()>=5) cout <<"XaodAnalysis::retrieveXaodMet    Building met for sys " << SusyNtSysNames.at(sys) << endl;

    m_metContainer = new xAOD::MissingETContainer;
    m_metAuxContainer = new xAOD::MissingETAuxContainer;
    m_metContainer->setStore(m_metAuxContainer);
    m_metContainer->reserve(10);

    xAOD::ElectronContainer* electrons = xaodElectrons(sysInfo, sys);
    xAOD::MuonContainer* muons = xaodMuons(sysInfo, sys);
    xAOD::JetContainer* jets = xaodJets(sysInfo, sys);
    xAOD::PhotonContainer* photons = xaodPhotons(sysInfo, sys);

    // GetMet(met, jet, elec, muon, gamma, tau, doTST = true, doJVT=true, invis = 0)
    m_susyObj[m_eleIDDefault]->GetMET(*m_metContainer,
                                        jets,
                                        electrons,
                                        muons,
                                        photons,
                                        0);

    if(dbg()>=5) cout << "XaodAnalysis::retrieveXaodMet    Built MET with "
            << electrons->size() << " electrons "
            << muons->size() << " muons "
            << jets->size() << " jets "
            << photons->size() << " photons" << endl;
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::retrieveXaodTrackMet(ST::SystInfo sysInfo, SusyNtSys sys)
{
    if(dbg()>=5) cout << "XaodAnalysis::retrieveXaodTrackMet    Building track met for sys " << SusyNtSysNames.at(sys) << endl;

    m_trackMetContainer = new xAOD::MissingETContainer;
    m_trackMetAuxContainer = new xAOD::MissingETAuxContainer;

    m_trackMetContainer->setStore(m_trackMetAuxContainer);
    m_trackMetContainer->reserve(10);

    xAOD::ElectronContainer* electrons = xaodElectrons(sysInfo, sys);
    xAOD::MuonContainer* muons = xaodMuons(sysInfo, sys);
    xAOD::JetContainer* jets = xaodJets(sysInfo, sys);

    // GetTrackMET(met, jet, elec, muon)
    m_susyObj[m_eleIDDefault]->GetTrackMET(*m_trackMetContainer,
                                            jets,
                                            electrons,
                                            muons);

    if(dbg()>=5) cout << "XaodAnalysis::retriveXaodTrackMet   Built track MET with "
                    << electrons->size() << " electrons "
                    << muons->size() << " muons "
                    << jets->size() << " jets" << endl;

}
//////////////////////////////////////////////////////////////////////////////
const xAOD::TruthParticleContainer* XaodAnalysis::retrieveTruthParticles(xAOD::TEvent& e, bool dbg)
{
    const xAOD::TruthParticleContainer* truth = nullptr;
    e.retrieve(truth, "TruthParticles");
    if(dbg) {
        if(truth) cout << "XaodAnalysis::retrieveTruthParticles    Retrieved truth particles (size: " << truth->size() << ")" << endl;
        else    cout << "XaodAnalysis::retrieveTruthParticles    WARNING Failed to retrieve truth particles" << endl;
    }
    return truth;
}
//////////////////////////////////////////////////////////////////////////////
const xAOD::TruthParticleContainer* XaodAnalysis::xaodTruthParticles()
{
    if(mc() && m_xaodTruthParticles==nullptr) {
        m_xaodTruthParticles = retrieveTruthParticles(m_event, dbg());
    }
    return m_xaodTruthParticles;
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::fill_objects(SusyNtSys sys, ST::SystInfo sysInfo)
{
    // fill the containers for:
    //  1) "pre" susyNt objects --> THESE ARE THE ONES THAT GET STORED IN THE OUTPUT NTUPLE
    //  2) "baseline" objects, those that base the default SUSYTools' baseline defintion
    fill_baseline_objects(sys, sysInfo);

    // fill the containers for:
    // 3) "signal" objects, those that base the default SUSYTools' signal definition
    fill_signal_objects(sys, sysInfo);
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::fill_baseline_objects(SusyNtSys sys, ST::SystInfo sysInfo)
{
    if(dbg()>=10) cout << "XaodAnalysis::fill_baseline_objects    Filling baseline objects (sys=" << SusyNtSysNames.at(sys) << ")" << endl;

    ///////////////////////////////////////////
    // containers
    ///////////////////////////////////////////
    xAOD::ElectronContainer* electrons = xaodElectrons(sysInfo, sys);
    xAOD::MuonContainer* muons = xaodMuons(sysInfo, sys);
    xAOD::JetContainer* jets = xaodJets(sysInfo, sys);
    xAOD::TauJetContainer* taus = xaodTaus(sysInfo, sys);

    ////////////////////////////////////////////
    // electrons
    ////////////////////////////////////////////
    if(electrons) {
        int iEl = -1;
        for(const auto& el : *electrons) {
            iEl++;
            if(dbg()>=10) cout << "XaodAnalysis::fill_baseline_objects     Electron[" << iEl << "]    (pt,eta,phi)=(" << el->pt() <<","<<el->eta() << "," << el->phi() << ")" << endl;

            // here we assume that the ST baseline definition is fine
            if( (bool)el->auxdata<char>("baseline")==1) m_baseElectrons.push_back(iEl);

            //if(dbg()>=10) cout << "XaodAnalysis::fill_baseline_objects    \t> passing baseline ? " << (bool)(el->auxdata<char>("baseline") << endl;

            /////////////////////////////////////////////////////////
            // STORING THIS ELECTRON
            const xAOD::CaloCluster* cluster = el->caloCluster();
            double et = cluster->e()/cosh(cluster->eta());
            // here we put a cut on the electron pT to store in the output susyNt
            if( et * MeV2GeV > 5 )
                m_preElectrons.push_back(iEl);
            /////////////////////////////////////////////////////////
        } // el
        if(dbg()>=10) cout << "XaodAnalysis::fill_baseline_objects    preElectrons size = " << m_preElectrons.size() << endl;
    }
    else {
        cout << "XaodAnalysis::fill_baseline_objects    WARNING Electrons container is null!" << endl;
    }
    

    ////////////////////////////////////////////
    // muons
    ////////////////////////////////////////////
    if(muons) {
        int iMu = -1;
        for(const auto& mu : *muons) {
            iMu++;

            /////////////////////////////////////////////////////////
            // STORING THIS MUON
            // here do not place any selection on muons that we store in the output susyNt
            m_preMuons.push_back(iMu);
            /////////////////////////////////////////////////////////

            if(dbg()>=10) cout << "XaodAnalysis::fill_baseline_objects    Muon[" << iMu << "]    (pt,eta,phi)=(" << mu->pt() << "," << mu->eta() << "," << mu->phi() << ")" << endl;

            if( (bool)mu->auxdata<char>("baseline")==1 ) m_baseMuons.push_back(iMu);
        } // mu
        if(dbg()>=10) cout << "XaodAnalysis::fill_baseline_objects    preMuons size = " << m_preMuons.size() << endl;
    }
    else {
        cout << "XaodAnalysis::fill_baseline_objects    WARNING Muons container is null!" << endl;
    }

    //////////////////////////////////
    // For updated jet selection (as of SUSY,2.3.15a)
    // we need OR flags for Jets in order to check for b-tagging
    // and "bad" jets
    //////////////////////////////////
    m_susyObj[m_eleIDDefault]->OverlapRemoval(electrons, muons, jets);

    ////////////////////////////////////////////
    // jets
    ////////////////////////////////////////////
    if(jets) {
        int iJet=-1;
        for(const auto& jet : *jets) {
            iJet++;

            /////////////////////////////////////////////////////////
            // STORING THIS JET
            // place a 20 GeV cut on jets that we store in the output susyNt
            if(jet->pt()*MeV2GeV > 20.0) m_preJets.push_back(iJet);
            /////////////////////////////////////////////////////////

            if(dbg()>=10) cout << "XaodAnalysis::fill_baseline_objects    Jet[" << iJet << "]    (pt,eta,phi)=(" << jet->pt() << "," << jet->eta() << "," << jet->phi() << ")" << endl;

            if( (bool)jet->auxdata<char>("baseline")==1 ) m_baseJets.push_back(iJet);
        } // jet
        if(dbg()>=10) cout << "XaodAnalysis::fill_baseline_objects    preJets size = " << m_preJets.size() << endl;
    }
    else {
        cout << "XaodAnalysis::fill_baseline_objects    WARNING Jets container is null!" << endl;
    }
    
    ////////////////////////////////////////////
    // taus
    ////////////////////////////////////////////
    if(taus) {
        int iTau=-1;
        for(const auto& tau : *taus) {
            iTau++;

            int nTracks = tau->nTracks();
            if( std::abs(tau->charge()) == 1 && (nTracks==1 || nTracks==3 || nTracks==5)) {
                m_contTaus.push_back(iTau);
            }
            if((bool)tau->auxdata<char>("baseline")==1) {
                m_preTaus.push_back(iTau);
                m_baseTaus.push_back(iTau);
            }
        } // tau
    }
    else {
        cout << "XaodAnalysis::fill_baseline_objects    WARNING Taus container is null!" << endl;
    }

    //////////////////////////////////
    // if nominal, keep track of idx
    //////////////////////////////////
    if(sys==NtSys::NOM) {
        m_preElectrons_nom = m_preElectrons;
        m_preMuons_nom = m_preMuons;
        m_preJets_nom = m_preJets;
        m_contTaus_nom = m_contTaus;
        m_preTaus_nom = m_preTaus;
    }
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::fill_signal_objects(SusyNtSys sys, ST::SystInfo sysInfo)
{
    if(dbg()>=10) cout << "XaodAnalysis::fill_signal_objects    Filling signal objects (sys=" << SusyNtSysNames.at(sys) << ")" << endl;

    /////////////////////////////////
    // grab the containers 
    /////////////////////////////////
    xAOD::ElectronContainer* electrons = xaodElectrons(sysInfo, sys);
    xAOD::MuonContainer* muons = xaodMuons(sysInfo, sys);
    xAOD::JetContainer* jets = xaodJets(sysInfo, sys);
    xAOD::TauJetContainer* taus = xaodTaus(sysInfo, sys);
    xAOD::PhotonContainer* photons = xaodPhotons(sysInfo, sys);

    ////////////////////////////////////////////
    // electrons
    ////////////////////////////////////////////
    if(electrons) {
        int iEl = 0;
        for(const auto& el : *electrons) {
            bool pass_signal_def = (bool)el->auxdata<char>("signal")==1;
            bool pass_OR = (bool)el->auxdata<char>("passOR")==1;
            if(dbg()>=10) cout << "XaodAnalysis::fill_signal_objects    Electron[" << iEl << "]    (pt,eta,phi)=(" << el->pt() << "," << el->eta() << "," << el->phi() << ")" << endl;
            if(pass_signal_def && pass_OR) m_sigElectrons.push_back(iEl);
            iEl++;
        } 
        if(dbg()>=10) cout << "XaodAnalysis::fill_signal_objects    sigElectrons size = " << m_sigElectrons.size() << endl; 
    }
    else {
        cout << "XaodAnalysis::fill_signal_objects    WARNING Electron container null" << endl;
    }
    
    ////////////////////////////////////////////
    // muons
    ////////////////////////////////////////////
    if(muons) {
        int iMu = 0;
        for(const auto& mu : *muons) {
            bool pass_signal_def = (bool)mu->auxdata<char>("signal")==1;
            bool pass_OR = (bool)mu->auxdata<char>("passOR")==1;
            if(dbg()>=10) cout << "XaodAnalysis::fill_signal_objects    Muon[" << iMu << "]    (pt,eta,phi)=(" << mu->pt() << "," << mu->eta() << "," << mu->phi() << ")" << endl;
            if(pass_signal_def && pass_OR) m_sigMuons.push_back(iMu);
            iMu++;
        }
        if(dbg()>=10) cout << "XaodAnalysis::fill_signal_objects    sigMuons size = " << m_sigMuons.size() << endl;
    }
    else {
        cout << "XaodAnalysis::fill_signal_objects    WARNING Muon container null" << endl;
    }

    ////////////////////////////////////////////
    // jets
    ////////////////////////////////////////////
    if(jets) {
        int iJet=0;
        for(const auto& jet : *jets) {
            bool pass_pt = jet->pt()*MeV2GeV > 20.;
            bool pass_OR = (bool)jet->auxdata<char>("passOR")==1;
            bool pass_bad = (!(bool)jet->auxdata<char>("bad")==1);
            if(dbg()>=10) cout << "XaodAnalysis::fill_signal_objects    Jet[" << iJet << "]    (pt,eta,phi)=(" << jet->pt() << "," << jet->eta() << "," << jet->phi() << ")" << endl;
            if(pass_pt && pass_OR && pass_bad) m_sigJets.push_back(iJet);
            iJet++;
        }
        if(dbg()>=10) cout << "XaodAnalysis::fill_signal_objects    sigJets size = " << m_sigJets.size() << endl;
    }
    else {
        cout << "XaodAnalysis::fill_signal_objects    WARNING Jet container null" << endl;
    }

    ////////////////////////////////////////////
    // taus
    ////////////////////////////////////////////
    if(taus) {
        int iTau = 0;
        for(const auto& tau : *taus) {
            bool pass_pt = tau->pt() * MeV2GeV > 20.;
            bool pass_signal = (bool)tau->auxdata<char>("signal")==1;
            if(dbg()>=10) cout << "XaodAnalysis::fill_signal_objects    Tau[" << iTau << "]    (pt,eta,phi)=(" << tau->pt() << "," << tau->eta() << "," << tau->phi() << ")" << endl;
            if(pass_pt & pass_signal) m_sigTaus.push_back(iTau);
            iTau++;
        }
        if(dbg()>=10) cout <<"XaodAnalysis::fill_signal_objects    sigTuas size = " << m_sigTaus.size() << endl;
    } 
    else {
        cout << "XaodAnalysis::fill_signal_objects    WARNING Tau container null" << endl;
    }

    ////////////////////////////////////////////
    // photons
    ////////////////////////////////////////////
    if(photons) {
    int iPh = 0;
    for(const auto& ph : *photons) {
        #warning need to set photon decorators in SUSYTools
        bool pass_cleaning = (bool)ph->auxdata<char>("passCleaning");
        bool pass_ambiguity = (bool)ph->auxdata<char>("passAmbiguity");
        if(pass_cleaning && pass_ambiguity) { m_prePhotons.push_back(iPh); m_basePhotons.push_back(iPh); m_sigPhotons.push_back(iPh); }
        iPh++;
    }
    if(dbg()>=10) cout << "XaodAnalysis::fill_signal_objects    sigPhotons size = " << m_sigPhotons.size() << endl;
    }
    else {
        if(dbg()) cout << "XaodAnalysis::fill_signal_objects    WARNING Photon container null" << endl;
    }
}
//////////////////////////////////////////////////////////////////////////////
void XaodAnalysis::sample_event_triggers()
{
    if(dbg()>=5) cout << "XaodAnalysis::sample_event_triggers    Checking if any of our triggers fired" << endl;
    m_evtTrigBits.ResetAllBits();
    vector<string> trigs = xaodTriggers();
    for(unsigned int itrig = 0; itrig < trigs.size(); itrig++) {
        if(m_susyObj[m_eleIDDefault]->IsTrigPassed(trigs[itrig])) m_evtTrigBits.SetBitNumber(itrig, true);
    }
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::eleIsOfType(const xAOD::Electron& in, ElectronId id)
{
    if     (id==ElectronId::VeryLooseLLH  && m_elecSelLikelihoodVeryLoose->accept(in))  return true;
    else if(id==ElectronId::LooseLLH  && m_elecSelLikelihoodLoose->accept(in))  return true;
    else if(id==ElectronId::LooseLLHBLayer && m_elecSelLikelihoodLooseBLayer->accept(in)) return true;
    else if(id==ElectronId::MediumLLH && m_elecSelLikelihoodMedium->accept(in)) return true;
    else if(id==ElectronId::TightLLH  && m_elecSelLikelihoodTight->accept(in))  return true;
    
    return false;
}
//////////////////////////////////////////////////////////////////////////////
TBits XaodAnalysis::matchElectronTriggers(const xAOD::Electron& in)
{
    // DA: maybe split up trigger list in terms of electron/muon/etc triggers? but then bit numbers
    // may be out of sync w.r.t. the stored histogram... in any case, the non-passed triggers are always
    // false and on SusyNtuple side the user will provide specific string for the ele/muo trigger that
    // is checked against the trig histo
    if(m_dbg>=15) cout << "XaodAnalysis::matchElectronTriggers" << endl;
    TBits eleTrigBits(m_nTriggerBits);
    eleTrigBits.ResetAllBits();
    std::vector<std::string> trigs = XaodAnalysis::xaodTriggers();
    for(unsigned int iTrig = 0; iTrig < trigs.size(); iTrig++) {
    
        std::size_t elfound = trigs[iTrig].find("HLT_e");
        std::size_t elfound2 = trigs[iTrig].find("HLT_2e");
        std::size_t mufound = trigs[iTrig].find("mu");
        std::size_t mufound2 = trigs[iTrig].find("MU");
        std::size_t metfound = trigs[iTrig].find("xe");
    
        bool mu_chain = ((mufound != std::string::npos) || (mufound2 != std::string::npos));
        bool el_chain = ((elfound != std::string::npos) || (elfound2 != std::string::npos));
        bool met_chain = (metfound != std::string::npos);
        bool dilepton = mu_chain && el_chain;
    
        bool ismatch = false;
        if( (el_chain || dilepton) && !met_chain)
            ismatch = m_susyObj[m_eleIDDefault]->IsTrigMatched(&in, trigs[iTrig]);
        if(ismatch) {
            //cout << "     > ele match to : " << trigs[iTrig] << endl;
            eleTrigBits.SetBitNumber(iTrig,true);
        }
    }
    return eleTrigBits; 
}
//////////////////////////////////////////////////////////////////////////////
TBits XaodAnalysis::matchMuonTriggers(const xAOD::Muon& in)
{
    if(dbg()>=15) cout << "XaodAnalysis::matchMuonTriggers" << endl;
    TBits muoTrigBits(m_nTriggerBits);
    muoTrigBits.ResetAllBits();
    std::vector<std::string> trigs = XaodAnalysis::xaodTriggers();
    for(unsigned int iTrig=0; iTrig<trigs.size(); iTrig++){
    
        std::size_t elfound = trigs[iTrig].find("el");
        std::size_t mufound = trigs[iTrig].find("mu");
        std::size_t mufound2 = trigs[iTrig].find("MU");
        std::size_t metfound = trigs[iTrig].find("xe");
    
        bool mu_chain = ((mufound != std::string::npos) || (mufound2 != std::string::npos));
        bool el_chain = (elfound != std::string::npos);
        bool met_chain = (metfound != std::string::npos);
        bool dilepton = mu_chain && el_chain;
    
        if( (mu_chain || dilepton) && !met_chain)
            if(m_susyObj[m_eleIDDefault]->IsTrigMatched(&in, trigs[iTrig])) muoTrigBits.SetBitNumber(iTrig, true);
    }
    return muoTrigBits;
}
//////////////////////////////////////////////////////////////////////////////
std::map<std::string, std::vector<unsigned int>> XaodAnalysis::getDiMuTrigMap(const xAOD::Muon &in, const xAOD::MuonContainer &muons)
{
    if(dbg()>=15) cout << "XaodAnalysis::getDiMuTrigMap" << endl;

    std::map<std::string, std::vector<unsigned int>> diMuTrigMap;

    std::vector<std::string> trigs = XaodAnalysis::xaodTriggers();
    for(unsigned int iTrig=0; iTrig<trigs.size(); ++iTrig){
        string trig = trigs[iTrig];

        // currently: dimuon trigger iff two instances of *lowercase* 'mu'
        // HLT_mu20_iloose_L1MU15 (e.g.) is a *single-muon* trigger
        int count = 0;
        string token = "mu";
        for (size_t offset = trig.find(token); offset != std::string::npos; offset = trig.find(token, offset + token.length())) {
            ++count;
        }

        if (count >= 2) {
            diMuTrigMap[trig] = {};
            //for (unsigned int i = 0; i < muons.size(); ++i) {
            //    if (m_susyObj[m_eleIDDefault]->IsTrigMatched(&in, muons[i], trig)) {
            //        diMuTrigMap[trig].push_back(i);
            //    }
            //} // i

            // dantrim -- loop over indices of pre muons container, since pre-muons are the ones
            //              that will get stored in the output susyNt file and whose indices
            //              we store
            for(auto &i : m_preMuons) {
                if (m_susyObj[m_eleIDDefault]->IsTrigMatched(&in, muons.at(i), trig)) {
                    diMuTrigMap[trig].push_back(i);
                }
            } // i
        } // >=2
    } // iTrig

    return diMuTrigMap;
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::muIsOfType(const xAOD::Muon& in, MuonId id)
{
    if     (id==MuonId::VeryLoose && m_muonSelectionToolVeryLoose->accept(in))  return true;
    else if(id==MuonId::Loose     && m_muonSelectionToolLoose    ->accept(in))  return true;
    else if(id==MuonId::Medium    && m_muonSelectionToolMedium   ->accept(in))  return true;
    else if(id==MuonId::Tight     && m_muonSelectionToolTight    ->accept(in))  return true;
    return false;
}
//////////////////////////////////////////////////////////////////////////////
bool XaodAnalysis::dilepton_trigger_matched(const xAOD::IParticle* part1, const xAOD::IParticle* part2,
            string trigger_chain)
{

    if(trigger_chain=="") {
        cout << "XaodAnalysis::dilepton_trigger_matched    ERROR Trigger chain is \"\", returning false" << endl;
        return false;
    }
    return m_susyObj[m_eleIDDefault]->IsTrigMatched({part1, part2}, trigger_chain);
}
