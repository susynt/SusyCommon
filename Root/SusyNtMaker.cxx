#include "SusyCommon/SusyNtMaker.h"

//SusyCommon
#include "SusyCommon/SusyObjId.h"
#include "SusyCommon/ss3l_chargeflip.h"

//SusyNtuple
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/TriggerTools.h"

//xAOD
//#include "EventPrimitives/EventPrimitivesHelper.h"
#include "AthContainers/AuxElement.h"
#include "xAODPrimitives/IsolationType.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackParticlexAODHelpers.h"
#include "xAODEgamma/EgammaxAODHelpers.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/MuonAuxContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODTau/TauxAODHelpers.h"

//Tools
#include "ElectronPhotonSelectorTools/AsgElectronChargeIDSelectorTool.h"

//SUSY
//#include "SUSYTools/SUSYCrossSection.h"

//std/stl
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <iostream>
using namespace std;

using Susy::SusyNtMaker;

using GhostList_t = std::vector< ElementLink<xAOD::IParticleContainer> >;
static SG::AuxElement::ConstAccessor<GhostList_t> ghostAcc("GhostTrack");

//////////////////////////////////////////////////////////////////////////////
SusyNtMaker::SusyNtMaker() :
    m_flags_checked(false),
    m_file_outtree(0),
    m_outtree(0),
    m_susyNt(0),
    m_susyFinalState(0)
{
    n_pre_ele = 0;
    n_pre_muo = 0;
    n_pre_tau = 0;
    n_pre_jet = 0;
    n_base_ele = 0;
    n_base_muo = 0;
    n_base_tau = 0;
    n_base_jet = 0;
    n_sig_ele = 0;
    n_sig_muo = 0;
    n_sig_tau = 0;
    n_sig_jet = 0;

}
//////////////////////////////////////////////////////////////////////////////
SusyNtMaker::~SusyNtMaker()
{

}
//////////////////////////////////////////////////////////////////////////////
struct FillCutFlow {
    int iCut; // index of sequential cut
    bool passAll; // whether we've survived all cuts so far
    bool includeThisCut; // whether this cut should be used when computed passAll
    vector<size_t> *counters;
    FillCutFlow(vector<size_t> *cs) :
        iCut(0), passAll(true), includeThisCut(true), counters(cs) {}
    FillCutFlow& operator()(bool thisEventDoesPassThisCut, float weight) {
        if(thisEventDoesPassThisCut && passAll) {
            counters->at(iCut) += 1;
        } else {
            if(includeThisCut) passAll = false;
        }
        iCut++;
        return *this;
    }
    FillCutFlow& disableFilterNextCuts() { includeThisCut = false; return *this; }
    FillCutFlow& enableFilterNextCuts() { includeThisCut = true; return *this; }
}; // struct
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::Init(TTree *tree)
{
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::SlaveBegin(TTree* tree)
{
    cout << "SusyNtMaker::SlaveBegin" << endl;
    XaodAnalysis::SlaveBegin(tree);
    XaodAnalysis::Init(tree);

    if(fill_nt()) {
        initialize_output_tree();
    }
    // dantrim 2017 June 11 - this is where we had the cutflow histograms, but they are useless these days...
    //initializeCutflowHistograms();
    initialize_counters();

    // start the analysis timer
    m_timer.Start();
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::initialize_counters()
{
    // event level trigger histo
    vector<string> trigs = XaodAnalysis::xaodTriggers();
    h_passTrigLevel = new TH1D("trig", "Event Level Triggers Fired", trigs.size()+1, 0.0, trigs.size()+1);
    h_passTrigLevel->GetXaxis()->SetLabelSize(0.8*h_passTrigLevel->GetLabelSize());
    for(uint itrig = 0; itrig < trigs.size(); itrig++)
        h_passTrigLevel->GetXaxis()->SetBinLabel(itrig+1, trigs[itrig].c_str());

    // cutflow
    m_cutstageCounters.clear();
    for(uint icut = 0; icut < cutflow_labels().size(); icut++)
        m_cutstageCounters.push_back(0);

}
//////////////////////////////////////////////////////////////////////////////
const vector<string> SusyNtMaker::cutflow_labels()
{
    vector<string> labels;
    labels.push_back("Initial");
    labels.push_back("GRL");
    labels.push_back("error flags");
    labels.push_back("good pvs");
    labels.push_back("bad muon");
    labels.push_back("cosmic muon");
    labels.push_back("jet cleaning");
    labels.push_back(">=1 base lepton");
    labels.push_back(">=1 signal lepton");
    return labels;
}
//////////////////////////////////////////////////////////////////////////////
Bool_t SusyNtMaker::Process(Long64_t entry)
{
    static Long64_t chainEntry = -1;
    chainEntry++;
    m_event.getEntry(chainEntry);

    // clear the output storage indices
    clear_output_objects();

    // start the PRW tool since all tools depend on it downstream
    for(int susyObjId : Susy::leptonIds()) {
        m_susyObj[susyObjId]->ApplyPRWTool();
        if(m_run_oneST) break;
    }

    const xAOD::EventInfo* eventinfo = XaodAnalysis::xaodEventInfo();

    if(dbg() || chainEntry % 500 == 0) {
        cout << " *** Processing entry " << setw(6) << chainEntry
                << "  run " << setw(6) << eventinfo->runNumber()
                << "  event " << setw(7) << eventinfo->eventNumber() << " *** " << endl;
    }


    // before filling check that things are consistent
    if(!m_flags_checked) {
        m_flags_checked = true;
        if(!running_options_are_valid()) {
            cout << "SusyNtMaker    Running options are invalid/inconsistent, exiting" << endl;
            abort();
        }
    }

    // fill our object containers with this event's objects
    retrieve_xaod_collections();

    // fill the event level trigger histo
    fill_event_trigger_histo();

    // SUSY final state
    susy_finalstate();

    ///////////////////////////////////////////////////////////
    // object selections
    ///////////////////////////////////////////////////////////

    // clear the SusyNtObject for the new event
    m_susyNt.clear();

    bool pass_event_selection = pass_event_level_selection();
    if(!pass_event_selection) return kTRUE;

    // get the nominal objects
    fill_nominal_objects();
    


    // clear the object containers before the next event
    delete_shallow_copies();
    //clear_containers();
    clear_output_objects();

    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
bool SusyNtMaker::running_options_are_valid()
{
    bool ok = true;
    bool is_sim = xaodEventInfo()->eventType( xAOD::EventInfo::IS_SIMULATION );
    bool is_data = !is_sim;
    if(mc() != is_sim) {
        ok = false;
        if(dbg()) {
            cout << "SusyNtMaker::running_options_are_valid    Invalid run options: "
                << "(SusyNtMaker isMC = " << mc() << ", xAOD::EventInfo::IS_SIMULATION = " << is_sim << ")"
                << endl;
        }
    }
    if(is_data) {
        const std::vector<xAOD::EventInfo::StreamTag> &streams = xaodEventInfo()->streamTags();
        vector<string> streamnames(streams.size());
        std::transform(streams.begin(), streams.end(), streamnames.begin(),
                [](const xAOD::EventInfo::StreamTag &s) { cout << "SusyNtMaker::running_options_are_valid    StreamTag " << s.name() << endl; return s.name(); });
        bool is_physicsMain = (find(streamnames.begin(), streamnames.end(), "Main") != streamnames.end());
        bool consistent_stream = (is_physicsMain ? m_stream==Stream_PhysicsMain : false);

        if(!consistent_stream) {
            ok = false;
            if(dbg()) {
                cout << "SusyNtMaker::running_options_are_valid    Inconsistent DataStream: "
                    << " SusyNtMaker stream: "
                    << (m_stream==Stream_PhysicsMain ? "Stream_PhysicsMain" : "Unknown")
                    << " EventInfo: "
                    << accumulate(streamnames.begin(), streamnames.end(), std::string(),
                        [](const std::string& a, const std::string& b) -> std::string {
                            return a + (a.length() > 0 ? "," : "") + b;
                        })
                    << endl;
            }
        }
    }
    if(dbg()) cout << "SusyNtMaker::running_options_are_valid    Valid options? " << ok << endl;
    return ok;
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::fill_event_trigger_histo()
{
    std::vector<std::string> trigs = XaodAnalysis::xaodTriggers();
    for(uint itrig = 0; itrig < trigs.size(); itrig++) {
        if(m_susyObj[m_eleIDDefault]->IsTrigPassed(trigs[itrig])) h_passTrigLevel->Fill(itrig+0.5);
    }
}
//////////////////////////////////////////////////////////////////////////////
int SusyNtMaker::susy_finalstate()
{
    int pdg1 = 0;
    int pdg2 = 0;
    m_susyFinalState = 0;
    if(mc() && !xaodTruthParticles()->empty() && xaodTruthParticles()!=nullptr) {
        m_susyObj[m_eleIDDefault]->FindSusyHP(xaodTruthParticles(), pdg1, pdg2);
    }
    if(pdg1 != 0 && pdg2 !=0) m_susyFinalState = SUSY::finalState(pdg1, pdg2);
    return m_susyFinalState;
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::Terminate()
{
    cout << "SusyNtMaker::Terminate" << endl;
    XaodAnalysis::Terminate();
    m_timer.Stop();

    // print cutflow
    cout << counter_summary() << endl;
    cout << timer_summary() << endl;

    if(fill_nt()) {
        save_output_tree();
    }

    return;
}
//////////////////////////////////////////////////////////////////////////////
string SusyNtMaker::counter_summary()
{
    ostringstream oss;
    oss << "-------------------------------------------------------------" << endl
        << " NtMaker Counter Summary " << endl
        << endl
        << " Object Counter " << endl
        << "   > pre ele   " << n_pre_ele << endl
        << "   > pre muo   " << n_pre_muo << endl
        << "   > pre tau   " << n_pre_tau << endl
        << "   > pre jet   " << n_pre_jet << endl
        << "   > base ele  " << n_pre_ele << endl
        << "   > base muo  " << n_pre_muo << endl
        << "   > base tau  " << n_pre_tau << endl
        << "   > base jet  " << n_pre_jet << endl
        << "   > sig ele   " << n_pre_ele << endl
        << "   > sig muo   " << n_pre_muo << endl
        << "   > sig tau   " << n_pre_tau << endl
        << "   > sig jet   " << n_pre_jet << endl
        << endl;

    oss << " Event Counter " << endl;
    vector<string> labels = cutflow_labels();
    struct shorter { bool operator()(const string& a, const string& b) { return a.size() < b.size(); }};
    size_t max_label_length = max_element(labels.begin(), labels.end(), shorter())->size();

    for(size_t i = 0; i < m_cutstageCounters.size(); ++i)
        oss << "   " << setw(max_label_length+2) << std::left<<labels[i] << m_cutstageCounters[i] << endl;
    oss << endl;
    oss << "-------------------------------------------------------------" << endl;
    return oss.str();
}
//////////////////////////////////////////////////////////////////////////////
string SusyNtMaker::timer_summary()
{
    double realTime = m_timer.RealTime();
    double cpuTime = m_timer.CpuTime();
    int hours = int(realTime / 3600);
    realTime -= hours * 3600;
    int min = int(realTime / 60);
    realTime -= min * 60;
    int sec = int(realTime);
    int nEventInput = m_cutstageCounters.front();
    int nEventOutput = m_outtree ? m_outtree->GetEntries() : -1;
    float speed = nEventInput / m_timer.RealTime()/1000;
    TString line1; line1.Form("Real %d:%02d:%02d, CPU %.3f", hours, min, sec, cpuTime);
    TString line2; line2.Form("%2.3f",speed);
    ostringstream oss;
    oss << "-------------------------------------------------------------" << endl;
    oss << " Number of events processed : " << nEventInput << endl
        << " Number of events saved     : " << nEventOutput << endl
        << " Analysis time              : " << line1 << endl
        << " Analysis speed [kHz]       : " << line2 << endl;
    oss << "-------------------------------------------------------------" << endl;
    return oss.str();

}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::initialize_output_tree()
{
    m_file_outtree = new TFile(output_name().c_str(), "recreate");
    m_outtree = new TTree("susyNt", "susyNt");
    m_outtree->SetAutoSave(10000000);
    //m_outtree->SetMaxTreeSize(3000000000u); // dantrim June 11 2017 - Steve's/Davide's value is significantly smaller than the default? Comment this out... perhaps this removes the multile output files
    m_susyNt.SetActive();
    m_susyNt.WriteTo(m_outtree);

    if(dbg())
        cout << "SusyNtMaker::initialize_output_tree    " << m_outtree << endl;

    return;
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::save_output_tree()
{
    cout << "SusyNtMaker::save_output_tree    " << m_file_outtree << "  " << m_outtree << endl;
    m_file_outtree = m_outtree->GetCurrentFile();
    m_file_outtree->Write(0, TObject::kOverwrite);
    cout << "SusyNtMaker::save_output_tree    susyNt saved to " << m_file_outtree->GetName() << endl;
    write_metadata();
    m_file_outtree->Close();
    return;
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::write_metadata()
{

    cout << "SusyNtMaker::write_metadata" << endl;
    struct {
        string operator()(const string &s) { return (s.size()==0 ? " WARNING empty string!" : ""); }
    } warn_if_empty;
    if(dbg()) {
        cout << "SusyNtMaker::write_metadata    Writing the following info to file: " << endl;
        cout << "SusyNtMaker::write_metadata     > input container name   : "
                << input_container() << warn_if_empty(input_container()) << endl;
        cout << "SusyNtMaker::write_metadata     > output container name  : "
                << output_container() << warn_if_empty(output_container()) << endl;
        cout << "SusyNtMaker::write_metadata     > production tag         : "
                << production_tag() << warn_if_empty(production_tag()) << endl;
        cout << "SusyNtMaker::write_metadata     > production command     : "
                << production_command() << warn_if_empty(production_command()) << endl; 
    }

    if(m_file_outtree) {
        TDirectory* current_directory = gROOT->CurrentDirectory();
        m_file_outtree->cd();
        TNamed input_cont_name("inputContainerName", input_container().c_str());
        TNamed output_cont_name("outputContainerName", output_container().c_str());
        TNamed prod_tag("productionTag", production_tag().c_str());
        TNamed production_cmd("productionCommand", production_command().c_str());
        input_cont_name.Write();
        output_cont_name.Write();
        prod_tag.Write();
        production_cmd.Write();
        current_directory->cd();
    }
    else {
        cout << "SusyNtMaker::write_metadata    WARNING Missing output file, cannot write metadata" << endl;
    }
    return;
}
//////////////////////////////////////////////////////////////////////////////
bool SusyNtMaker::pass_event_level_selection()
{
    const xAOD::EventInfo* eventinfo = XaodAnalysis::xaodEventInfo();
    float w = mc() ? eventinfo->mcEventWeight() : 1;

    fill_event_cleaning_flags();

    FillCutFlow fillCutFlow(&m_cutstageCounters);

    bool keep_all_events = !do_event_filter();

    bool pass_grl(m_cutFlags & ECut_GRL);
    bool pass_lar(m_cutFlags & ECut_LarErr);
    bool pass_tile(m_cutFlags & ECut_TileErr);
    bool pass_TTC(m_cutFlags & ECut_TTC);
    bool pass_SCT(m_cutFlags & ECut_SCTErr);
    bool pass_errorFlags(pass_lar && pass_tile && pass_TTC && pass_SCT);

    fillCutFlow(true, w); // initial
    fillCutFlow(pass_grl, w);
    fillCutFlow(pass_errorFlags, w);

    if(dbg()>=5 && !(keep_all_events || fillCutFlow.passAll)) {
        cout << "SusyNtMaker::pass_event_level_selection    "
                << "Event " << eventinfo->eventNumber() << " failed event level selection" << endl;
    }

    //return (keep_all_events || fillCutFlow.passAll);
    #warning bypassing event level selection filtering
    return true;
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::fill_nominal_objects()
{
    SusyNtSys sys = NtSys::NOM;
    ST::SystInfo sysInfo = systInfoList[0];
    fill_objects(sys, sysInfo);
}
