#include "SusyCommon/SusyNtMaker.h"

//SusyCommon
#include "SusyCommon/SusyObjId.h"
#include "SusyCommon/ss3l_chargeflip.h"

//SusyNtuple
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/TriggerTools.h"
#include "SusyNtuple/RecoTruthClassification.h" // isFakeLepton

//xAOD
#include "EventPrimitives/EventPrimitivesHelpers.h" // Amg::error
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
    n_pre_pho = 0;
    n_base_ele = 0;
    n_base_muo = 0;
    n_base_tau = 0;
    n_base_jet = 0;
    n_base_pho = 0;
    n_sig_ele = 0;
    n_sig_muo = 0;
    n_sig_tau = 0;
    n_sig_jet = 0;
    n_sig_pho = 0;

    CP::SystematicCode::enableFailure();

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
////////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::SlaveBegin(TTree* tree)
{
    if(dbg()) cout << "SusyNtMaker::SlaveBegin" << endl;
    XaodAnalysis::SlaveBegin(tree);

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
    labels.push_back(">=1 pre lepton");
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

    // start the PRW tool since all tools depend on it downstream
    // (specifically the RandomRunNumber being attached to EventInfo)
    for(int susyObjId : Susy::leptonIds()) {
        m_susyObj[susyObjId]->ApplyPRWTool();
        if(m_run_oneST) break;
    }

    const xAOD::EventInfo* eventinfo = XaodAnalysis::xaodEventInfo();
    if(dbg() || chainEntry % 5000 == 0) {
        cout << "SusyNtMaker::Process     *** Processing entry " << setw(6) << chainEntry
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

    // collect the SUSY final state (if there is one)
    susy_finalstate();

    ///////////////////////////////////////////////////////////
    // object selections
    ///////////////////////////////////////////////////////////

    // clear the SusyNtObject for the new event
    m_susyNt.clear();

    bool pass_event_selection = pass_event_level_selection();
    if(!pass_event_selection) { clear_event(); return kTRUE; }

    // get the nominal objects
    fill_nominal_objects();

    // apply selection on objects
    bool pass_object_selection = pass_object_level_selection();
    if(!pass_object_selection) { clear_event(); return kTRUE; }

    ////////////////////////////////////////////////////////////
    // store the objects in the output SusyNtObject
    ////////////////////////////////////////////////////////////
    if(fill_nt()) {

        if(mc()) m_tauTruthMatchingTool->initializeEvent(); // gives the tool the truth info

        // dantrim June 12 2017 -- TODO update the dilepton trigger matching to be Event
        // store the event-wide trigger bits to Susy::Event (not doing trigger object matching here)
        sample_event_triggers();

        // store objects to the output susyNt
        fill_nt_variables();

        // run systematics
        if(mc() && sys()) {
            run_kinematic_systematics();
        }

        // fill the output tree
        int bytes = m_outtree->Fill();
        if(bytes < 0) {
            cout << "SusyNtMaker::Process    ERROR Unable to fill output tree, abort (TTree::Fill returns " << bytes << ")" << endl;
            abort();
        }
    }

    clear_event();

    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::clear_event()
{
    // clear the object containers before the next event
    delete_shallow_copies();
    //clear_containers();
    clear_output_objects();
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
vector<int> SusyNtMaker::susy_finalstate()
{
    vector<int> out;
    int pdg1 = 0;
    int pdg2 = 0;
    m_susyFinalState = 0;
    if(mc() && !xaodTruthParticles()->empty() && xaodTruthParticles()!=nullptr) {
        m_susyObj[m_eleIDDefault]->FindSusyHP(xaodTruthParticles(), pdg1, pdg2);
    }
    if(pdg1 != 0 && pdg2 !=0) m_susyFinalState = SUSY::finalState(pdg1, pdg2);
    out.push_back(m_susyFinalState);
    out.push_back(pdg1);
    out.push_back(pdg2);
    return out;
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::Terminate()
{
    if(dbg()) cout << "SusyNtMaker::Terminate" << endl;
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
        << "   > pre pho   " << n_pre_pho << endl
        << " - - - - - - - - - - - - - - - - " << endl
        << "   > base ele  " << n_base_ele << endl
        << "   > base muo  " << n_base_muo << endl
        << "   > base tau  " << n_base_tau << endl
        << "   > base jet  " << n_base_jet << endl
        << "   > base pho  " << n_base_pho << endl
        << " - - - - - - - - - - - - - - - - " << endl
        << "   > sig ele   " << n_sig_ele << endl
        << "   > sig muo   " << n_sig_muo << endl
        << "   > sig tau   " << n_sig_tau << endl
        << "   > sig jet   " << n_sig_jet << endl
        << "   > sig pho   " << n_sig_pho << endl
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
//////////////////////////////////////////////////////////////////////////////
bool SusyNtMaker::pass_object_level_selection()
{
    const xAOD::EventInfo* eventinfo = XaodAnalysis::xaodEventInfo();
    double w = mc() ? eventinfo->mcEventWeight() : 1;

    fill_object_cleaning_flags(); // bad jet, bad muon, cosmic muon

    ////////////////////////////////////////
    // update cutflow
    ////////////////////////////////////////
    FillCutFlow fillCutFlow(&m_cutstageCounters);
    // we've filled up the event counters
    fillCutFlow.iCut = 3;

    bool pass_jet_cleaning(m_cutFlags & ECut_BadJet);
    bool pass_good_pv(m_cutFlags & ECut_GoodVtx);
    bool pass_bad_muon(m_cutFlags & ECut_BadMuon);
    bool pass_cosmic(m_cutFlags & ECut_Cosmic);

    bool pass_ge1pl( (m_preElectrons.size() + m_preMuons.size()) >= 1);
    bool pass_ge1bl( (m_baseElectrons.size() + m_baseMuons.size()) >= 1);
    bool pass_ge1sl( (m_sigElectrons.size() + m_sigMuons.size()) >= 1);

    fillCutFlow(pass_good_pv, w);
    fillCutFlow(pass_bad_muon, w);
    fillCutFlow(pass_cosmic, w);
    fillCutFlow(pass_jet_cleaning, w);
    fillCutFlow(pass_ge1pl, w);
    fillCutFlow(pass_ge1bl, w);
    fillCutFlow(pass_ge1sl, w);

    ////////////////////////////////////////
    // apply filtering
    ////////////////////////////////////////
    bool event_passes = true;
    if(do_event_filter()) {
        bool pass_nlep_filter = true;
        bool pass_trig_filter = true;

        // filter on pre lepton objects
        if(nlep_for_filter()>0) {
            int n_lep = (m_preElectrons.size() + m_preMuons.size());
            pass_nlep_filter = (n_lep >= nlep_for_filter()); 
            //cout << "filtering leptons : " << n_lep << "  pass_nlep: " << pass_nlep_filter << "   pass_g21pl? " << pass_ge1pl <<  endl;
            event_passes = (event_passes && pass_nlep_filter);
        }

        // simply filter by requiring any of our event level triggers to have fired (and we have many)
        if(do_trig_filter()) {
            pass_trig_filter = (h_passTrigLevel->Integral() > 0);
            event_passes = (event_passes && pass_trig_filter);
        }
        //cout << "do object filter? " << do_event_filter() << "  nlep_for_filter: " << nlep_for_filter() << "  pass_nlep_filter: " << pass_nlep_filter << "  do trig filter: " << do_trig_filter() << "  pass trig? " << pass_trig_filter << "  EVENT PASSES ? " << event_passes << endl;
    }
    if(dbg()>=5 && !event_passes) {
        cout << "SusyNtMaker::pass_object_level_selection    Event " << eventinfo->eventNumber() << " passes object selection" << endl;
    }

    /////////////////////////////////////
    // update object counters
    /////////////////////////////////////
    if(event_passes) {
        // pre objects
        n_pre_ele   += m_preElectrons.size();
        n_pre_muo   += m_preMuons.size();
        n_pre_tau   += m_preTaus.size();
        n_pre_jet   += m_preJets.size();
        n_pre_pho   += m_prePhotons.size();

        // base objects
        n_base_ele  += m_baseElectrons.size();
        n_base_muo  += m_baseMuons.size();
        n_base_tau  += m_baseTaus.size();
        n_base_jet  += m_baseJets.size();
        n_base_pho  += m_basePhotons.size();

        // signal objects
        n_sig_ele   += m_sigElectrons.size();
        n_sig_muo   += m_sigMuons.size();
        n_sig_tau   += m_sigTaus.size();
        n_sig_jet   += m_sigJets.size();
        n_sig_pho   += m_sigPhotons.size();
    }

    return event_passes;

}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::fill_nt_variables()
{

    // Susy::Event
    fill_event_variables();

    // Susy::Electron
    fill_electron_variables();

    // Susy::Muon
    fill_muon_variables();

    // Susy::Jet
    fill_jet_variables();

    // Susy::Tau
    fill_tau_variables();

    // Susy::Photon
    fill_photon_variables();

    // Susy::Met
    fill_met_variables();

    // Susy::TrackMet
    fill_track_met_variables();

    #warning NEED TO FILL TRUTH VARIABLES

    // perform dilepton trigger matching 
    // WARNING this must come AFTER or AT THE END OF 'fill_nt_variables'
    perform_dilepton_trigger_matching();
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::fill_event_variables()
{
    if(dbg()>=5) cout << "SusyNtMaker::fill_event_variables    Filling Susy::Event" << endl;

    Susy::Event* evt = m_susyNt.evt();
    const xAOD::EventInfo* eventinfo = XaodAnalysis::xaodEventInfo();

    // run coordinates
    evt->run = eventinfo->runNumber();
    evt->eventNumber = eventinfo->eventNumber();
    evt->lb = eventinfo->lumiBlock();
    evt->stream = m_stream;
    int year = 0;
    if(mc()) {
        year = m_susyObj[m_eleIDDefault]->treatAsYear();
    }
    else {
        if(data15()) year = 2015;
        else if(data16()) year = 2016;
    }
    if(year==0) cout << "SusyNtMaker::fill_event_variables    WARNING treatAsYear was not found correctly for event " << eventinfo->eventNumber() << "!" << endl;
    evt->treatAsYear = year;

    // <mu>
    evt->avgMu = m_susyObj[m_eleIDDefault]->GetCorrectedAverageInteractionsPerCrossing(false);
    evt->avgMuDataSF = m_susyObj[m_eleIDDefault]->GetCorrectedAverageInteractionsPerCrossing(true);

    evt->isMC = mc();
    evt->mcChannel = mc() ? eventinfo->mcChannelNumber() : 0;
    evt->w = mc() ? eventinfo->mcEventWeight() : 1;

    //cout << "SusyNtMaker::fill_event_variables    mc channel : " << evt->mcChannel << "  year: " << year << endl;

    evt->initialNumberOfEvents = m_nEventsProcessed;
    evt->sumOfEventWeights = m_sumOfWeights;
    evt->sumOfEventWeightsSquared = m_sumOfWeightsSquared;

    evt->nVtx = xaodVertices()->size();

    vector<int> susyinfo = susy_finalstate();
    evt->susyFinalState = susyinfo.at(0);
    evt->susySpartId1 = susyinfo.at(1);
    evt->susySpartId2 = susyinfo.at(2);

    // dantrim June 12 2017 - TODO items to remove from Susy::Event
    // hDecay, eventWithSusyProp

    evt->trigBits = m_evtTrigBits;

    //////////////////////////////////////////
    // PRW
    //////////////////////////////////////////
    evt->wPileup = mc() ? m_susyObj[m_eleIDDefault]->GetPileupWeight() : 1;
    if(mc() && sys()) {
        for(const auto& sysInfo : systInfoList) {
            if(!(sysInfo.affectsType == ST::SystObjType::EventWeight && sysInfo.affectsWeights)) continue;
            const CP::SystematicSet& sys = sysInfo.systset;
            SusyNtSys ourSys = CPsys2sys((sys.name()).c_str());
            if(!(ourSys==NtSys::PILEUP_UP || ourSys==NtSys::PILEUP_DN)) continue;
            bool do_down = ourSys==NtSys::PILEUP_DN;
            // configure the tools
            if ( m_susyObj[m_eleIDDefault]->applySystematicVariation(sys) != CP::SystematicCode::Ok) {
                cout << "SusyNtMaker::fill_event_variables    cannot configure SUSYTools for systematic " << sys.name() << " (" << SusyNtSysNames.at(ourSys) << ")" << endl;
                continue;
            }
            // get and store
            if(do_down) { evt->wPileup_dn = m_susyObj[m_eleIDDefault]->GetPileupWeight(); }
            else { evt->wPileup_up = m_susyObj[m_eleIDDefault]->GetPileupWeight(); }
        } // sysInfo
        if( m_susyObj[m_eleIDDefault]->resetSystematics() != CP::SystematicCode::Ok) {
            cout << "SusyNtMaker::fill_event_variables    cannot reset SUSYTools systematics. Aborting." << endl;
            abort();
        }
    } // mc && sys

    // sherpa 2.2 V+jets weight
    float vjetweight = 1.0;
    bool is_sherpa22vjet = false;
    if(mc()) {
        int mc_number = eventinfo->mcChannelNumber();
        if(SusyNtTools::isSherpa22Vjet(mc_number)) {
            vjetweight = m_susyObj[m_eleIDDefault]->getSherpaVjetsNjetsWeight();
            is_sherpa22vjet = true;
        }
    } // mc
    evt->isSherpaVjetsSample = is_sherpa22vjet;
    evt->sherpa22VjetsWeight = vjetweight;

    if(mc())
        evt->mcWeights = get_mc_weights(eventinfo);

    m_susyNt.evt()->cutFlags[NtSys::NOM] = m_cutFlags;
}
//////////////////////////////////////////////////////////////////////////////
vector<float> SusyNtMaker::get_mc_weights(const xAOD::EventInfo* ei)
{
    vector<float> weights;
    for(const float w : ei->mcEventWeights()) {
        weights.push_back(w);
    }
    // alternative method:
    //const xAOD::TruthEventContainer *truthE;
    //CHECK( m_event.retrieve( truthE, "TruthEvents" ) );
    //const xAOD::TruthEvent *truthEvent = (*truthE)[0];
    //for (const float weight : truthEvent->weights()) {
    //    truthWeights.push_back(weight);
    //}
    return weights;
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::fill_electron_variables()
{
    if(dbg()>=5) cout << "SusyNtMaker::fill_electron_variables    Filling Susy::Electron" << endl;

    xAOD::ElectronContainer* electrons = xaodElectrons(systInfoList[0]);

    if(electrons) {
        if(mc() && m_derivation.Contains("SUSY")) {
            int dsid = xaodEventInfo()->mcChannelNumber();
            fillElectronChargeFlip(electrons, xaodTruthParticles(), dsid);
        }
        for(auto& i : m_preElectrons) {
            store_electron(*(electrons->at(i)), i);
        }
    }
    else {
        cout << "SusyNtMaker::fill_electron_variables    WARNING Electron container is null" << endl;
    }
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::fill_muon_variables()
{
    if(dbg()>=5) cout << "SusyNtMaker::fill_muon_variables    Filling Susy::Muon" << endl;

    xAOD::MuonContainer* muons = xaodMuons(systInfoList[0]);

    if(muons) {
        for(auto &i : m_preMuons) {
            store_muon(*(muons->at(i)), *muons);
        } // i
    //    cout << "---------------------------------" << endl;
    // xAOD::MuonContainer *sf_muon = new xAOD::MuonContainer;
    // xAOD::MuonAuxContainer *sf_muon_aux = new xAOD::MuonAuxContainer;
    // sf_muon->setStore(sf_muon_aux);
    // xAOD::Muon* sfMu = new xAOD::Muon;
    // sfMu->makePrivateStore(in);
    // sf_muon->push_back(sfMu);
    //    if(m_preMuons.size() == 2 && (m_preMuons.size() + m_preElectrons.size())==2 ) {
    //        xAOD::MuonContainer* test_mu = new xAOD::MuonContainer;
    //        xAOD::MuonAuxContainer* test_mu_aux = new xAOD::MuonAuxContainer;
    //        test_mu->setStore(test_mu_aux);
    //        for(auto &i : m_preMuons) {
    //            xAOD::Muon* mu = new xAOD::Muon;
    //            mu->makePrivateStore(*(muons->at(i)));
    //            test_mu->push_back(mu);
    //            cout << "  >> mu with pt = " << mu->pt() * MeV2GeV << endl;
    //        }
    //        std::string trig = "HLT_mu24_iloose_L1MU15";
    //        cout << "trigger sf " << m_susyObj[m_eleIDDefault]->GetTotalMuonTriggerSF(*test_mu, trig) << endl;
    //    }

    }
    else {
        cout << "SusyNtMaker::fill_muon_variables    WARNING Muon container is null" << endl;
    }
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::fill_jet_variables()
{
    if(dbg()>=15) cout << "SusyNtMaker::fill_jet_variables    Filling Susy::Jet" << endl;

    xAOD::JetContainer* jets = xaodJets(systInfoList[0]);

    if(jets) {
        // b-jet SF
        // dantrim June 13 2017 - TODO check if better to do this on per-jet basis inside (as with muon trigger SF)
        if(mc()) m_susyObj[m_eleIDDefault]->BtagSF(jets); // decorates jets with effscalefact
        for(auto& i : m_preJets_nom) {
            store_jet(*(jets->at(i)));
        }
    }
    else {
        cout << "SusyNtMaker::fill_jet_variables    WARNING Jet container is null" << endl;
    }

}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::fill_tau_variables()
{
    if(dbg()>=15) cout << "SusyNtMaker::fill_tau_variables    Filling Susy::Tau" << endl;

    xAOD::TauJetContainer* taus = xaodTaus(systInfoList[0]);

    if(taus) {
        vector<int>& saveTaus = m_saveContTaus ? m_contTaus : m_preTaus;
        for(auto& i: saveTaus) {
            store_tau(*(taus->at(i)));
        }
    }
    else {
        cout << "SusyNtMaker::fill_tau_variables    WARNING Tau container is null" << endl;
    }
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::fill_photon_variables()
{
    if(dbg()>=15) cout << "SusyNtMaker::fill_photon_variables    Filling Susy::Photon" << endl;

    xAOD::PhotonContainer* photons = xaodPhotons(systInfoList[0]);

    if(photons) {
        for(auto& i : m_prePhotons) {
            store_photon(*(photons->at(i)));
        }
    }
    else {
        cout << "SusyNtMaker::fill_photon_variables    WARNING Photon container is null" << endl;
    }
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::fill_met_variables(SusyNtSys sys)
{

    if(dbg()>=15) cout << "SusyNtMaker::fill_met_variables    Filling MET (sys=" << SusyNtSysNames.at(sys) << ")" << endl;

    xAOD::MissingETContainer::const_iterator met_it = xaodMET()->find("Final");

    if(dbg()>=15) {
        cout << "SusyNtMaker::fill_met_variables    Dumping xAOD MET container " << endl;
        for(auto it = xaodMET()->begin(), end = xaodMET()->end(); it!=end; ++it) {
            cout << "SusyNtMaker::fill_met_variables      > " << (*it)->name() << endl;
        }
    }

    if(met_it == xaodMET()->end()) {
        cout << "SusyNtMaker::fill_met_variables    WARNING No RefFinal inside MET container, not filling MET variables" << endl;
        return;
    }

    Susy::Met* met = new Susy::Met();
    met->Et = (*met_it)->met()*MeV2GeV;
    met->phi = (*met_it)->phi();
    met->sumet = (*met_it)->sumet()*MeV2GeV;
    met->sys = sys;

    if(dbg()>=15) cout << "SusyNtMaker::fill_met_variables    (sys=" << SusyNtSysNames.at(sys) << ") (Et,phi,pt)=("<<met->Et<<","<<met->phi<<","<<met->lv().Pt()<<")"<<endl;

    ////////////////////////////////////////////
    // Electron Term (RefEle)
    ////////////////////////////////////////////
    xAOD::MissingETContainer::const_iterator met_find = xaodMET()->find("RefEle");
    if(met_find != xaodMET()->end()) {
        met->refEle_et = (*met_find)->met()*MeV2GeV;
        met->refEle_phi = (*met_find)->phi();
        met->refEle_sumet = (*met_find)->sumet()*MeV2GeV;
    }

    ////////////////////////////////////////////
    // Photon Term (RefGamma)
    ////////////////////////////////////////////
    met_find = xaodMET()->find("RefGamma");
    if(met_find != xaodMET()->end()) {
        met->refGamma_et = (*met_find)->met()*MeV2GeV;
        met->refGamma_phi = (*met_find)->phi();
        met->refGamma_sumet = (*met_find)->sumet()*MeV2GeV;
    }


    ////////////////////////////////////////////
    // Tau Term (RefTau)
    ////////////////////////////////////////////
    met_find = xaodMET()->find("RefTau");
    if(met_find != xaodMET()->end()) {
        met->refTau_et = (*met_find)->met()*MeV2GeV;
        met->refTau_phi = (*met_find)->phi();
        met->refTau_sumet = (*met_find)->sumet()*MeV2GeV;

    }

    ////////////////////////////////////////////
    // Jet Term (RefJet)
    ////////////////////////////////////////////
    met_find = xaodMET()->find("RefJet");
    if (met_find != xaodMET()->end()) {
        met->refJet_et = (*met_find)->met()*MeV2GeV;
        met->refJet_phi = (*met_find)->phi();
        met->refJet_sumet = (*met_find)->sumet()*MeV2GeV;
    }

    ////////////////////////////////////////////
    // SoftTerm
    ////////////////////////////////////////////
    met_find = xaodMET()->find("PVSoftTrk"); // using TST (track soft term, not PVSoftClus)
    if (met_find != xaodMET()->end()) {
        met->softTerm_et = (*met_find)->met()*MeV2GeV;
        met->softTerm_phi = (*met_find)->phi();
        met->softTerm_sumet = (*met_find)->sumet()*MeV2GeV;
    }

    ////////////////////////////////////////////
    // Muon Term (RefMuons)
    ////////////////////////////////////////////
    met_find = xaodMET()->find("Muons");
    if (met_find != xaodMET()->end()) {
        met->refMuo_et = (*met_find)->met()*MeV2GeV;
        met->refMuo_phi = (*met_find)->phi();
        met->refMuo_sumet = (*met_find)->sumet()*MeV2GeV;
    }

    //__________________ DONE WITH MET _____________________//
    m_susyNt.met()->push_back(*met);

    if( dbg()>20 && (sys!=NtSys::NOM)) {
        cout << "SusyNtMaker::fill_met_variables       > Nominal MET: " << m_susyNt.met()->at(0).Et << ",  Sys varied MET: " << m_susyNt.met()->back().Et << endl;
    }
    delete met;
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::fill_track_met_variables(SusyNtSys sys)
{

    if(dbg()>=15) cout << "SusyNtMaker::fill_track_met_variables    Filling TrackMet (sys="<<SusyNtSysNames.at(sys)<<")" <<endl;

    xAOD::MissingETContainer::const_iterator trackMet_it = xaodTrackMET()->find("Track");

    if(trackMet_it == xaodTrackMET()->end()) {
        cout << "SusyNtMaker::fill_track_met_variables    WARNING Cannot find 'Track' inside MET_Track container, unable to fill TrackMet" << endl;
        return;
    }

    Susy::TrackMet* tmet = new Susy::TrackMet();

    tmet->Et =  (*trackMet_it)->met()*MeV2GeV;// m_met.Et();
    tmet->phi = (*trackMet_it)->phi();// m_met.Phi();
    tmet->sys = sys;
    tmet->sumet = (*trackMet_it)->sumet()*MeV2GeV;

    //////////////////////////////////////////////
    // Electron Term
    //////////////////////////////////////////////
    xAOD::MissingETContainer::const_iterator met_find = xaodTrackMET()->find("RefEle");
    if(met_find != xaodTrackMET()->end()) {
        tmet->refEle_et = (*met_find)->met()*MeV2GeV;
        tmet->refEle_phi = (*met_find)->phi();
        tmet->refEle_sumet = (*met_find)->sumet()*MeV2GeV;
    }

    //////////////////////////////////////////////
    // Muon Term
    //////////////////////////////////////////////
    met_find = xaodTrackMET()->find("Muons");
    if(met_find != xaodTrackMET()->end()) {
        tmet->refMuo_et = (*met_find)->met()*MeV2GeV;
        tmet->refMuo_phi = (*met_find)->phi();
        tmet->refMuo_sumet = (*met_find)->sumet()*MeV2GeV;
    }

    //////////////////////////////////////////////
    // Jet Term
    //////////////////////////////////////////////
    met_find = xaodTrackMET()->find("RefJet");
    if(met_find != xaodTrackMET()->end()) {
        tmet->refJet_et = (*met_find)->met()*MeV2GeV;
        tmet->refJet_phi = (*met_find)->phi();
        tmet->refJet_sumet = (*met_find)->sumet()*MeV2GeV;
    }

    //////////////////////////////////////////////
    // Soft Term
    //////////////////////////////////////////////
    met_find = xaodTrackMET()->find("PVSoftTrk");
    if(met_find != xaodTrackMET()->end()) {
        tmet->softTerm_et = (*met_find)->met()*MeV2GeV;
        tmet->softTerm_phi = (*met_find)->phi();
        tmet->softTerm_sumet = (*met_find)->sumet()*MeV2GeV;
    } 

    //__________________ DONE WITH TRACK MET _____________________//
    m_susyNt.tkm()->push_back(*tmet);

    delete tmet;
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::store_electron(const xAOD::Electron& in, int ele_idx)
{
    if(dbg()>=15) cout << "SusyNtMaker::store_electron    Electron[" << ele_idx << "], pt=" << in.pt()*MeV2GeV << endl;

    const xAOD::EventInfo* eventinfo = XaodAnalysis::xaodEventInfo();

    Susy::Electron out;

    //////////////////////////////////
    // 4 vector
    //////////////////////////////////
    double pt = ( (in.pt()*MeV2GeV < 0) ? 0 : in.pt()*MeV2GeV );
    double m = ( (in.m()*MeV2GeV < 0) ? 0 : in.m()*MeV2GeV );
    double eta = in.eta();
    double phi = in.phi();
    out.SetPtEtaPhiM(pt, eta, phi, m);
    out.pt = pt;
    out.eta = eta;
    out.phi = phi;
    out.m = m;
    out.q = in.charge();

    //////////////////////////////////
    // electron author
    //////////////////////////////////
    out.author = static_cast<int>(in.author());
    out.authorElectron = in.author() & xAOD::EgammaParameters::AuthorElectron;
    out.authorAmbiguous = in.author() & xAOD::EgammaParameters::AuthorAmbiguous;

    //////////////////////////////////
    // SUSYTools paramters
    //////////////////////////////////
    out.isBaseline = in.auxdata<char>("baseline");
    out.isSignal = in.auxdata<char>("signal");

    //////////////////////////////////
    // ID
    //////////////////////////////////
    out.veryLooseLLH = eleIsOfType(in, ElectronId::VeryLooseLLH);
    out.looseLLH = eleIsOfType(in, ElectronId::LooseLLH);
    out.mediumLLH = eleIsOfType(in, ElectronId::MediumLLH);
    out.looseLLHBLayer = eleIsOfType(in, ElectronId::LooseLLHBLayer);
    out.tightLLH = eleIsOfType(in, ElectronId::TightLLH);
    
    //////////////////////////////////
    // OQ
    //////////////////////////////////
    out.passOQBadClusElectron = (in.isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON) ? true : false);

    //////////////////////////////////
    // Isolation Selection
    //////////////////////////////////
    out.isoGradientLoose          = m_isoToolGradientLooseTight->accept(in) ? true : false;
    out.isoGradient               = m_isoToolGradientTightCalo->accept(in) ? true : false;
    out.isoLooseTrackOnly         = m_isoToolLooseTrackOnlyLoose->accept(in) ? true : false;
    out.isoLoose                  = m_isoToolLoose->accept(in) ? true : false;
    out.isoFixedCutTightTrackOnly = m_isoToolTight->accept(in) ? true : false;
    
    //////////////////////////////////
    // Isolation Variables
    //////////////////////////////////
    out.etconetopo20 = in.isolationValue(xAOD::Iso::topoetcone20) * MeV2GeV;
    out.etconetopo30 = in.isolationValue(xAOD::Iso::topoetcone30) * MeV2GeV;
    out.ptcone20 = in.isolationValue(xAOD::Iso::ptcone20) * MeV2GeV;
    out.ptcone30 = in.isolationValue(xAOD::Iso::ptcone30) * MeV2GeV;
    out.ptvarcone20 = in.auxdataConst<float>("ptvarcone20") * MeV2GeV;
    out.ptvarcone30 = in.auxdataConst<float>("ptvarcone30") * MeV2GeV;

    if(dbg()>=15) {
        cout << "SusyNtMaker::store_electron    Electron[" << ele_idx << "]    pt=" << out.pt
                << "  LLH: (veryLooseLLH,looseLLH,mediumLLH,looseLLHBLayer,tightLLH)=("
                << out.veryLooseLLH<<","<<out.looseLLH<<","<<out.mediumLLH<<","<<out.looseLLHBLayer
                <<","<<out.tightLLH <<")" << endl;
    }

    //////////////////////////////////
    // EleEle shared track
    //////////////////////////////////
    int max_idx = 10;
    for(const auto iel : m_preElectrons) {
        if(iel >= max_idx) max_idx = iel+1;
    } // iel
    out.sharedEleEleTrk.resize(max_idx, 0);
    out.sharedEleEleTrk.assign(max_idx, 0);

    for(const auto iel : m_preElectrons) {
        if(iel == ele_idx) continue;
        xAOD::ElectronContainer* ele = xaodElectrons(systInfoList[0]);
        if(in.trackParticleLink() == ele->at(iel)->trackParticleLink()) {
            out.sharedEleEleTrk.at(iel) = 1;
        }
        else {
            out.sharedEleEleTrk.at(iel) = 0;
        }
    } // iel

    //////////////////////////////////
    // EleMu shared track
    //////////////////////////////////
    out.sharedMuTrk.resize(m_preMuons.size(), 0);
    out.sharedMuTrk.assign(m_preMuons.size(), 0);
    const xAOD::TrackParticle* elTrk = xAOD::EgammaHelpers::getOriginalTrackParticle(&in);
    if(elTrk) {
        xAOD::MuonContainer* muons = xaodMuons(systInfoList[0]);
        for(int im = 0; im < (int)m_preMuons.size(); im++) {
            const xAOD::TrackParticle* muTrk = muons->at(m_preMuons[im])->trackParticle(xAOD::Muon::InnerDetectorTrackParticle);
            if(muTrk) {
                if(elTrk == muTrk)
                    out.sharedMuTrk.at(im) = 1;
            }
            else { out.sharedMuTrk.at(im) = 0; }
        } // im
    }

    //////////////////////////////////
    // Scale Factors
    //////////////////////////////////
    bool recoSF=true;
    bool idSF=true;
    bool trigSF=false;
    bool isoSF=true;

    string single_ele = "singleLepton";
    string double_ele = "diLepton";
    string mixed_ele = "mixedLepton";

    if(mc()) {
        //////////////////////////////////////
        // Lepton SF
        // - one for each electron LH WP
        // - (only Loose, Medium, Tight for now)
        //////////////////////////////////////
        // signature: (input electron, bool doRecoSF, bool doIDSF, bool doTrigSF, bool doIsoSF, string trigExpr)
        // default trigExpr: "e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose"

        
        if(in.pt()*MeV2GeV >= 25.) {
            if(m_run_oneST) {
                out.eleEffSF[ElectronId::TightLLH] =  m_susyObj[m_eleIDDefault]->GetSignalElecSF (in, recoSF, idSF, trigSF, isoSF);

                out.eleTrigSF_single[ElectronId::TightLLH] = m_susyObj[m_eleIDDefault]->GetSignalElecSF(in, false, false, true, false, single_ele);
                out.eleTrigSF_double[ElectronId::TightLLH] = m_susyObj[m_eleIDDefault]->GetSignalElecSF(in, false, false, true, false, double_ele);
                out.eleTrigSF_mixed[ElectronId::TightLLH] = m_susyObj[m_eleIDDefault]->GetSignalElecSF(in, false, false, true, false, mixed_ele);

                //out.eleTrigSF[ElectronId::TightLLH] = m_susyObj[m_eleIDDefault]->GetSignalElecSF (in, false, false, true, false, ele_trig);
                //out.eleCHFSF[ElectronId::TightLLH] = m_susyObj[m_eleIDDefault]->GetSignalElecSF (in, false, false, false, false, ele_trig, true);
            }
            else {
                out.eleEffSF[ElectronId::TightLLH] =  m_susyObj[SusyObjId::eleTightLLH]->GetSignalElecSF (in, recoSF, idSF, trigSF, isoSF);

                out.eleTrigSF_single[ElectronId::TightLLH] = m_susyObj[SusyObjId::eleTightLLH]->GetSignalElecSF(in, false, false, true, false, single_ele);
                out.eleTrigSF_double[ElectronId::TightLLH] = m_susyObj[SusyObjId::eleTightLLH]->GetSignalElecSF(in, false, false, true, false, double_ele);
                out.eleTrigSF_mixed[ElectronId::TightLLH] = m_susyObj [SusyObjId::eleTightLLH]->GetSignalElecSF(in, false, false, true, false, mixed_ele);
                //out.eleCHFSF[ElectronId::TightLLH] = m_susyObj[SusyObjId::eleTightLLH]->GetSignalElecSF (in, false, false, false, false, ele_trig, true);

        
                out.eleEffSF[ElectronId::MediumLLH] = m_susyObj[SusyObjId::eleMediumLLH]->GetSignalElecSF(in, recoSF, idSF, trigSF, isoSF);

                out.eleTrigSF_single[ElectronId::MediumLLH] = m_susyObj[SusyObjId::eleMediumLLH]->GetSignalElecSF(in, false, false, true, false, single_ele);
                out.eleTrigSF_double[ElectronId::MediumLLH] = m_susyObj[SusyObjId::eleMediumLLH]->GetSignalElecSF(in, false, false, true, false, double_ele);
                out.eleTrigSF_mixed [ElectronId::MediumLLH] = m_susyObj[SusyObjId::eleMediumLLH]->GetSignalElecSF(in, false, false, true, false, mixed_ele);
                //out.eleCHFSF[ElectronId::MediumLLH] = m_susyObj[SusyObjId::eleMediumLLH]->GetSignalElecSF(in, false, false, false, false, ele_trig, true);
                //
                //cout << "SusyNtMaker::store_electron    NOMINAL TriggerSF   single = " << out.eleTrigSF_single[ElectronId::MediumLLH]
                //            << "   double = " << out.eleTrigSF_double[ElectronId::MediumLLH]
                //            << "   mixed  = " << out.eleTrigSF_mixed [ElectronId::MediumLLH] << endl;
            }
        }
    } // mc


    //////////////////////////////////
    // Truth matching info
    //////////////////////////////////
    if(mc()) {
        out.mcType = xAOD::TruthHelpers::getParticleTruthType(in);
        out.mcOrigin = xAOD::TruthHelpers::getParticleTruthOrigin(in);
        const xAOD::TruthParticle* truthEle = xAOD::TruthHelpers::getTruthParticle(in);
        out.matched2TruthLepton   = truthEle ? true : false;
        int matchedPdgId = truthEle ? truthEle->pdgId() : -999;
        out.truthType  = isFakeLepton(out.mcOrigin, out.mcType, matchedPdgId);
        out.truthCharge =  truthEle ? truthEle->charge() : 0;
        if(m_derivation.Contains("SUSY")) {
            out.ss3lChargeFlip = in.auxdataConst<int>("chargeFlip");
        
            bool pass_charge_id = m_electronChargeIDTool->accept(in);
            out.passChargeFlipTagger = pass_charge_id;
            out.chargeFlipBDT = m_electronChargeIDTool->calculate(&in, -99); // mu = -99 means will grab mu from EventInfo
        
            // Electron bkg origins
            out.mcBkgMotherPdgId = in.auxdata<int>("bkgMotherPdgId");
            out.mcBkgTruthOrigin = in.auxdata<int>("bkgTruthOrigin");
        }
    } // mc

    //////////////////////////////////
    // Systematic Variations on SF
    //////////////////////////////////
    if(mc() && sys() && (out.veryLooseLLH || out.looseLLH || out.mediumLLH || out.tightLLH) && in.pt()*MeV2GeV > 25.) {
        for(const auto& sysInfo : systInfoList) {
            if(!(sysInfo.affectsType == ST::SystObjType::Electron && sysInfo.affectsWeights)) continue;
    
            const CP::SystematicSet& sys = sysInfo.systset;
            SusyNtSys ourSys = CPsys2sys((sys.name()).c_str());
            for(int i : Susy::electronIds()){
                int index_to_check = (m_run_oneST==true ? (int)m_eleIDDefault : i);
                if(m_susyObj[index_to_check]->applySystematicVariation(sys) != CP::SystematicCode::Ok) {
                    cout << "SusyNtMaker::storeElectron    cannot configure SUSYTools for systematic " << sys.name() << endl;
                    continue;
                }
                if(m_run_oneST) break;
            } // i
            vector<float> sf;
            sf.assign(ElectronId::ElectronIdInvalid, 1);
            vector<float> sf_chf;
            sf_chf.assign(ElectronId::ElectronIdInvalid, 1);
    
            std::string nan_trig = "";
    
            if(m_run_oneST) {
                sf[ElectronId::TightLLH]  = m_susyObj[m_eleIDDefault] ->GetSignalElecSF(in, recoSF, idSF, trigSF, isoSF);
                //sf_chf[ElectronId::TightLLH] = m_susyObj[m_eleIDDefault] ->GetSignalElecSF(in, false, false, false, false, nan_trig, true);
            }
            else {
                sf[ElectronId::TightLLH]  = m_susyObj[SusyObjId::eleTightLLH] ->GetSignalElecSF(in, recoSF, idSF, trigSF, isoSF);
                sf[ElectronId::MediumLLH] = m_susyObj[SusyObjId::eleMediumLLH]->GetSignalElecSF(in, recoSF, idSF, trigSF, isoSF);
    
                //sf_chf[ElectronId::TightLLH] = m_susyObj[m_eleIDDefault] ->GetSignalElecSF(in, false, false, false, false,  nan_trig, true);
                //sf_chf[ElectronId::MediumLLH] = m_susyObj[m_eleIDDefault] ->GetSignalElecSF(in, false, false, false, false, nan_trig, true);
            }

            vector<float> sf_trig_single;
            vector<float> sf_trig_double;
            vector<float> sf_trig_mixed;
            sf_trig_single.assign(ElectronId::ElectronIdInvalid, 1);
            sf_trig_double.assign(ElectronId::ElectronIdInvalid, 1);
            sf_trig_mixed.assign(ElectronId::ElectronIdInvalid, 1);
            if(m_run_oneST) {
                sf_trig_single[ElectronId::TightLLH]  = m_susyObj[m_eleIDDefault]->GetSignalElecSF(in, false, false, true, false, single_ele);
                sf_trig_single[ElectronId::MediumLLH] = m_susyObj[m_eleIDDefault]->GetSignalElecSF(in, false, false, true, false, single_ele);

                sf_trig_double[ElectronId::TightLLH]  = m_susyObj[m_eleIDDefault]->GetSignalElecSF(in, false, false, true, false, double_ele);
                sf_trig_double[ElectronId::MediumLLH] = m_susyObj[m_eleIDDefault]->GetSignalElecSF(in, false, false, true, false, double_ele);

                sf_trig_mixed[ElectronId::TightLLH]  = m_susyObj[m_eleIDDefault]->GetSignalElecSF(in, false, false, true, false, mixed_ele);
                sf_trig_mixed[ElectronId::MediumLLH] = m_susyObj[m_eleIDDefault]->GetSignalElecSF(in, false, false, true, false, mixed_ele);
            }
            else {
                sf_trig_single[ElectronId::TightLLH]  = m_susyObj[SusyObjId::eleTightLLH] ->GetSignalElecSF(in, false, false, true, false, single_ele);
                sf_trig_single[ElectronId::MediumLLH] = m_susyObj[SusyObjId::eleMediumLLH]->GetSignalElecSF(in, false, false, true, false, single_ele);

                sf_trig_double[ElectronId::TightLLH]  = m_susyObj[SusyObjId::eleTightLLH] ->GetSignalElecSF(in, false, false, true, false, double_ele);
                sf_trig_double[ElectronId::MediumLLH] = m_susyObj[SusyObjId::eleMediumLLH]->GetSignalElecSF(in, false, false, true, false, double_ele);

                sf_trig_mixed[ElectronId::TightLLH]  = m_susyObj[SusyObjId::eleTightLLH] ->GetSignalElecSF(in, false, false, true, false, mixed_ele);
                sf_trig_mixed[ElectronId::MediumLLH] = m_susyObj[SusyObjId::eleMediumLLH]->GetSignalElecSF(in, false, false, true, false, mixed_ele);
            }

            for(int i=ElectronId::TightLLH; i<ElectronIdInvalid; i++){
                if     (ourSys == NtSys::EL_EFF_ID_TOTAL_Uncorr_UP)      out.errEffSF_id_up[i]   = sf[i] - out.eleEffSF[i];
                else if(ourSys == NtSys::EL_EFF_ID_TOTAL_Uncorr_DN)      out.errEffSF_id_dn[i]   = sf[i] - out.eleEffSF[i];
                else if(ourSys == NtSys::EL_EFF_Reco_TOTAL_Uncorr_UP)    out.errEffSF_reco_up[i] = sf[i] - out.eleEffSF[i];
                else if(ourSys == NtSys::EL_EFF_Reco_TOTAL_Uncorr_DN)    out.errEffSF_reco_dn[i] = sf[i] - out.eleEffSF[i];
                else if(ourSys == NtSys::EL_EFF_Iso_TOTAL_Uncorr_UP)     out.errEffSF_iso_up[i]  = sf[i] - out.eleEffSF[i];
                else if(ourSys == NtSys::EL_EFF_Iso_TOTAL_Uncorr_DN)     out.errEffSF_iso_dn[i]  = sf[i] - out.eleEffSF[i];

                else if(ourSys == NtSys::EL_EFF_Trigger_TOTAL_UP) {
                    out.errEffSF_trig_up_single[i] = sf_trig_single[i] - out.eleTrigSF_single[i];
                    out.errEffSF_trig_up_double[i] = sf_trig_double[i] - out.eleTrigSF_double[i];
                    out.errEffSF_trig_up_mixed[i] = sf_trig_mixed[i] - out.eleTrigSF_mixed[i];
                    //if(i==ElectronId::MediumLLH) {
                    //    cout << "SusyNtMaker::store_electron    UP     TriggerSF   single = " << out.eleTrigSF_single[ElectronId::MediumLLH]
                    //                << "   double = " << out.eleTrigSF_double[ElectronId::MediumLLH]
                    //                << "   mixed  = " << out.eleTrigSF_mixed [ElectronId::MediumLLH] << endl;
                    //}
                }
                else if(ourSys == NtSys::EL_EFF_Trigger_TOTAL_DN) {
                    out.errEffSF_trig_dn_single[i] = sf_trig_single[i] - out.eleTrigSF_single[i];
                    out.errEffSF_trig_dn_double[i] = sf_trig_double[i] - out.eleTrigSF_double[i];
                    out.errEffSF_trig_dn_mixed[i] = sf_trig_mixed[i] - out.eleTrigSF_mixed[i];
                    //if(i==ElectronId::MediumLLH) {
                    //    cout << "SusyNtMaker::store_electron    DN     TriggerSF   single = " << sf_trig_single[i] 
                    //                << "   double = " << sf_trig_double[i] 
                    //                << "   mixed  = " << sf_trig_mixed[i] << endl; 
                    //}
                }
        
//                if(i==0 || i==1 || i==2){
//                if(i==0)
//                    cout << "ElectronId : TightLH " <<  endl;
//                else if(i==1)
//                    cout << "ElectronId : MediumLH " <<  endl;
//                else if(i==2)
//                    cout << "ElectronId : LooseLH " <<  endl;
//        
//                cout << "   effId           : " << out.errEffSF_id_up[i] << "  " << out.errEffSF_id_dn[i] << endl;
//                cout << "   effReco         : " << out.errEffSF_reco_up[i] << "  " << out.errEffSF_reco_dn[i] << endl;
//                cout << "   effIso          : " << out.errEffSF_iso_up[i] << "  " << out.errEffSF_iso_dn[i] << endl;
//                cout << "   effTrig         : " << out.errEffSF_trig_up[i] << "  " << out.errEffSF_trig_dn[i] << endl;
//                }
        
            }
        } // sysInfo

        for(int i : Susy::electronIds()){
            int index_to_check = (m_run_oneST==true ? (int)m_eleIDDefault : i);
            if(m_susyObj[index_to_check]->resetSystematics() != CP::SystematicCode::Ok){
                cout << "SusyNtMaker::storeElectron    cannot reset SUSYTools systematics. Aborting." << endl;
                abort();
            }
            if(m_run_oneST) break;
        }
    } // if isMC
    else {
        for(int i=ElectronId::TightLLH; i<ElectronIdInvalid; i++){
            out.errEffSF_id_up[i] = out.errEffSF_id_dn[i] = 0;
            out.errEffSF_reco_up[i] = out.errEffSF_reco_dn[i] = 0;
        }
    }

    //////////////////////////////////////
    // Electron cluster information
    //////////////////////////////////////
    if(const xAOD::CaloCluster* c = in.caloCluster()) {
        out.clusE   = c->e()*MeV2GeV;
        out.clusEtaBE = c->etaBE(2);
        out.clusPhiBE = c->phiBE(2);
        out.clusEta = c->eta();
        out.clusPhi = c->phi();
    }

    //////////////////////////////////////
    // Electron track information
    //////////////////////////////////////
    if(const xAOD::TrackParticle* t = in.trackParticle()){
        out.trackPt = t->pt()*MeV2GeV;
        out.trackEta = t->eta();
        out.d0      = t->d0();
        out.d0sigBSCorr = xAOD::TrackingHelpers::d0significance( t, eventinfo->beamPosSigmaX(),
                                        eventinfo->beamPosSigmaY(), eventinfo->beamPosSigmaXY() );
    
        const xAOD::Vertex* PV = m_susyObj[m_eleIDDefault]->GetPrimVtx();
        double  primvertex_z = (PV) ? PV->z() : -999;
        out.z0 = t->z0() + t->vz() - primvertex_z;
    
        out.errD0         = Amg::error(t->definingParametersCovMatrix(),0);
        out.errZ0         = Amg::error(t->definingParametersCovMatrix(),1);
    }

    //////////////////////////////////////
    // Trigger Matching
    //////////////////////////////////////
    out.trigBits = matchElectronTriggers(in);
//    cout << "testing electron trigBits" << endl;
//    int nbins = h_passTrigLevel->GetXaxis()->GetNbins();
//    for(int iTrig=0; iTrig<46; iTrig++){
//        bool bit = out.trigBits.TestBitNumber(iTrig);
//        string trigger = h_passTrigLevel->GetXaxis()->GetBinLabel(iTrig+1);
//        cout << "\t passed trigger [" << iTrig << "] " << trigger << "? " << (bit ? "yes" : "no") << endl;
//    }
//    cout << endl;


    //______________ ALL DONE WITH THE ELECTRON ______________ //
    out.idx = (m_susyNt.ele()->size());
    m_susyNt.ele()->push_back(out);
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::store_muon(const xAOD::Muon& in, const xAOD::MuonContainer& muons)
{
    if(dbg()>=15) cout << "SusyNtMaker::store_muon   Muon (pt=" << in.pt()*MeV2GeV << ")" << endl; 

    const xAOD::EventInfo* eventinfo = xaodEventInfo();

    // testing
    //std::string test_chain = "HLT_mu8noL1";
    //float test_eff_data = m_susyObj[m_eleIDDefault]->GetMuonTriggerEfficiency(in, test_chain, true);
    //float test_eff_mc   = m_susyObj[m_eleIDDefault]->GetMuonTriggerEfficiency(in, test_chain, false);
    //cout << "SusyNtMaker::store_muon    test muon eff for " << test_chain << "  data = " << test_eff_data << "  mc = " << test_eff_mc << endl;

    //std::string test_mixed = "HLT_e17_lhloose_mu14";
    //std::string test_single = "HLT_mu14";

    //float test_mixed_eff = m_susyObj[m_eleIDDefault]->GetMuonTriggerEfficiency(in, test_mixed, false);
    //float test_single_eff = m_susyObj[m_eleIDDefault]->GetMuonTriggerEfficiency(in, test_single, false);
    ////cout << "SusyNtMaker::store_muon   test muon eff for mixed " << test_mixed << " = " << test_mixed_eff << ", for single " << test_single << " = " << test_single_eff << "  ---> equal? " << ( (test_mixed_eff==test_single_eff) ? "YES" : "NO") << endl;
    // xAOD::MuonContainer *sf_muon = new xAOD::MuonContainer;
    // xAOD::MuonAuxContainer *sf_muon_aux = new xAOD::MuonAuxContainer;
    // sf_muon->setStore(sf_muon_aux);
    // xAOD::Muon* sfMu = new xAOD::Muon;
    // sfMu->makePrivateStore(in);
    // sf_muon->push_back(sfMu);
    //float st_mixed = m_susyObj[m_eleIDDefault]->GetTotalMuonTriggerSF(*sf_muon, test_mixed);
    //cout << "SusyNtMaker::store_muon sf returned by SUSYTools for " << test_mixed << " = " << st_mixed << endl;
    //delete sf_muon;
    //delete sf_muon_aux;

    

    Susy::Muon out;

    //////////////////////////////////////
    // 4-vector
    //////////////////////////////////////
    double pt = ( (in.pt()*MeV2GeV < 0) ? 0 : in.pt()*MeV2GeV);
    double m = ( (in.m()*MeV2GeV < 0) ? 0 : in.m()*MeV2GeV);
    double eta(in.eta()), phi(in.phi());
    out.SetPtEtaPhiM(pt, eta, phi, m);
    out.pt  = pt;
    out.eta = eta;
    out.phi = phi;
    out.m   = m;
    out.q   = in.charge();
    
    //////////////////////////////////////
    // SUSYTools flags
    //////////////////////////////////////
    out.isBaseline = (bool)in.auxdata< char >("baseline");
    out.isSignal   = (bool)in.auxdata< char >("signal");
    out.isCaloTagged = (bool)in.muonType()==xAOD::Muon::CaloTagged;
    out.isSiForward = (bool)in.muonType()==xAOD::Muon::SiliconAssociatedForwardMuon;
    out.isCombined = in.muonType()==xAOD::Muon::Combined;
    out.isCosmic   = (bool)in.auxdata< char >("cosmic");  // note: this depends on definition of baseline and OR!
    out.isBadMuon  = (bool)in.auxdata<char>("bad");       // note: independent of definition of baseline/OR

    //////////////////////////////////////
    // Muon ID
    //////////////////////////////////////
    static SG::AuxElement::Accessor<float> mePt_acc("MuonSpectrometerPt");
    static SG::AuxElement::Accessor<float> idPt_acc("InnerDetectorPt");
    bool mu_has_decorations =  mePt_acc.isAvailable(in) && idPt_acc.isAvailable(in);
    if(mu_has_decorations) {
        out.veryLoose   = muIsOfType(in, MuonId::VeryLoose);
        out.loose       = muIsOfType(in, MuonId::Loose);
        out.medium      = muIsOfType(in, MuonId::Medium);
        out.tight       = muIsOfType(in, MuonId::Tight);
    }

    //////////////////////////////////////
    // Isolation selection
    //////////////////////////////////////
    out.isoGradientLoose          = m_isoToolGradientLooseTight->accept(in) ? true : false;
    out.isoGradient               = m_isoToolGradientTightCalo->accept(in) ? true : false;
    out.isoLooseTrackOnly         = m_isoToolLooseTrackOnlyLoose->accept(in) ? true : false;
    out.isoLoose                  = m_isoToolLoose->accept(in) ? true : false;
    out.isoFixedCutTightTrackOnly = m_isoToolTight->accept(in) ? true : false;

    bool all_available = true;
    //////////////////////////////////////
    // Isolation variables
    //////////////////////////////////////
    all_available &= in.isolation(out.ptcone20, xAOD::Iso::ptcone20); out.ptcone20 *= MeV2GeV;
    all_available &= in.isolation(out.ptcone30, xAOD::Iso::ptcone30); out.ptcone30 *= MeV2GeV;
    out.ptvarcone20 = in.auxdataConst<float>("ptvarcone20") * MeV2GeV;
    out.ptvarcone30 = in.auxdataConst<float>("ptvarcone30") * MeV2GeV;
    out.etconetopo20 = in.isolation(xAOD::Iso::topoetcone20) * MeV2GeV;
    out.etconetopo30 = in.isolation(xAOD::Iso::topoetcone30) * MeV2GeV;

    //////////////////////////////////////
    // Muon Track 
    //////////////////////////////////////
    const xAOD::TrackParticle* track;
    if(in.muonType()==xAOD::Muon::SiliconAssociatedForwardMuon) {
        track = in.trackParticle(xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle);
        if(!track) { track = nullptr; track = 0; }
    } // SAF
    else {
        track = in.primaryTrackParticle();
    }
    if(track) {
        const xAOD::Vertex* PV = m_susyObj[m_eleIDDefault]->GetPrimVtx();
        double  primvertex_z = (PV) ? PV->z() : 0.;
        out.d0             = track->d0();
        out.errD0          = Amg::error(track->definingParametersCovMatrix(),0);
        // add protection against missing muon track covariance matrix
        try {
            out.d0sigBSCorr = xAOD::TrackingHelpers::d0significance( track, eventinfo->beamPosSigmaX(),
                                        eventinfo->beamPosSigmaY(), eventinfo->beamPosSigmaXY() );
        }
        catch (...) {
            out.d0sigBSCorr = -99.;
            cout << "SusyNtMaker::store_muon    WARNING Exception caught from d0significance calculation (event #" << eventinfo->eventNumber() << ", muon pT=" << in.pt()*MeV2GeV << "), setting d0sigBSCorr to -99" << endl;
        }
        out.z0             = track->z0() + track->vz() - primvertex_z;
        out.errZ0          = Amg::error(track->definingParametersCovMatrix(),1);
    }

    //////////////////////////////////////
    // Muon ID Track
    //////////////////////////////////////
    if(const xAOD::TrackParticle* idtrack = in.trackParticle( xAOD::Muon::InnerDetectorTrackParticle )){
        out.idTrackPt      = idtrack->pt()*MeV2GeV;
        out.idTrackEta     = idtrack->eta();
        out.idTrackPhi     = idtrack->phi();
        out.idTrackQ       = idtrack->qOverP() < 0 ? -1 : 1;
        out.idTrackQoverP  = idtrack->qOverP()*MeV2GeV;
        out.idTrackTheta   = idtrack->theta();
    }

    //////////////////////////////////////
    // Muon MS Track
    //////////////////////////////////////
    if(const xAOD::TrackParticle* mstrack = in.trackParticle( xAOD::Muon::MuonSpectrometerTrackParticle )){
        out.msTrackPt      = mstrack->pt()*MeV2GeV;
        out.msTrackEta     = mstrack->eta();
        out.msTrackPhi     = mstrack->phi();
        out.msTrackQ       = mstrack->qOverP() < 0 ? -1 : 1;
        out.msTrackQoverP  = mstrack->qOverP()*MeV2GeV;
        out.msTrackTheta   = mstrack->theta();
    }

    //////////////////////////////////////
    // Muon Truth Matching/Info
    //////////////////////////////////////
    if(mc()) {
        const xAOD::TrackParticle* trackParticle = in.primaryTrackParticle();
        if(trackParticle){
            // mcType <==> "truthType" of input xAOD (MCTruthClassifier)
            out.mcType = xAOD::TruthHelpers::getParticleTruthType(*trackParticle);
            // mcOrigin <==> "truthOrigin" of input xAOD (MCTruthClassifier)
            out.mcOrigin = xAOD::TruthHelpers::getParticleTruthOrigin(*trackParticle);
            const xAOD::TruthParticle* truthMu = xAOD::TruthHelpers::getTruthParticle(*trackParticle);
            out.matched2TruthLepton = truthMu ? true : false;
            int matchedPdgId = truthMu ? truthMu->pdgId() : -999;
            out.truthType = isFakeLepton(out.mcOrigin, out.mcType, matchedPdgId);
        }
    }

    //////////////////////////////////////
    // Ghost Association of ID Track
    //////////////////////////////////////
    const xAOD::TrackParticle* idtrack = in.trackParticle( xAOD::Muon::InnerDetectorTrackParticle );
    out.ghostTrack.resize(m_preJets.size(), 0);
    out.ghostTrack.assign(m_preJets.size(), 0);
    if(idtrack) {
        xAOD::JetContainer* jets = XaodAnalysis::xaodJets(systInfoList[0]);
        for(int ij = 0; ij < (int)m_preJets.size(); ij++) {
            for(const auto& ghostLink : ghostAcc(*(jets->at(m_preJets[ij]))) ) {
                if(ghostLink.isValid() && (idtrack == *ghostLink)) {
                    out.ghostTrack.at(ij) = 1;
                    break; // move to next jet
                }
            } // ghostLink loop
        } // preJets loop
    } else {
        out.ghostTrack.assign(m_preJets.size(), 0);
    }

    //////////////////////////////////////
    // Trigger Matching
    //////////////////////////////////////
    out.trigBits = matchMuonTriggers(in);
//    cout << "testing muon trigBits" << endl;
//    int nbins = h_passTrigLevel->GetXaxis()->GetNbins();
//    for(int iTrig=0; iTrig<46; iTrig++){
//        bool bit = out.trigBits.TestBitNumber(iTrig);
//        string trigger = h_passTrigLevel->GetXaxis()->GetBinLabel(iTrig+1);
//        cout << "\t passed trigger " << trigger << "? " << (bit ? "yes" : "no") << endl;
//    }
//    cout << endl;
//    out.diMuTrigMap = getDiMuTrigMap(in, muons);

    //////////////////////////////////////
    // Lepton SF
    // - one for each MuonId that we use:
    // - Loose and Medium
    //////////////////////////////////////
    bool recoSF = true;
    bool isoSF = true;
    if(mc() && fabs(out.eta)<2.5 && out.pt>20){ // SF's are not binned for pt < 20 GeV
        if(m_run_oneST) {
            out.muoEffSF[MuonId::Loose] = m_susyObj[m_eleIDDefault]->GetSignalMuonSF(in, recoSF, isoSF);
        }
        else {
            out.muoEffSF[MuonId::Loose]  = m_susyObj[SusyObjId::muoLoose]->GetSignalMuonSF (in, recoSF, isoSF);
            out.muoEffSF[MuonId::Medium] = m_susyObj[SusyObjId::muoMedium]->GetSignalMuonSF(in, recoSF, isoSF);
        }

        // dantrim Jan 5 2015 : trigger SF kludge -- this is not absolutely correct as the SF are meant for the final signal muons
        // going into your selection, which is not the case here as we are forcing the tool to provide us
        // the SF on a per-muon basis
        xAOD::MuonContainer *sf_muon = new xAOD::MuonContainer;
        xAOD::MuonAuxContainer *sf_muon_aux = new xAOD::MuonAuxContainer;
        sf_muon->setStore(sf_muon_aux);
        xAOD::Muon* sfMu = new xAOD::Muon;
        sfMu->makePrivateStore(in);
        sf_muon->push_back(sfMu);

        TString trig_exp_med = "HLT_mu20_iloose_L1MU15_OR_HLT_mu50";
        if(m_run_oneST) {
            if(m_susyObj[m_eleIDDefault]->treatAsYear()==2016)
                trig_exp_med = "HLT_mu24_imedium";
        }
        else {
            if(m_susyObj[SusyObjId::muoMedium]->treatAsYear()==2016)
                trig_exp_med = "HLT_mu24_imedium";
        }
        trig_exp_med = "HLT_mu22_mu8noL1";
        // dantrim Sept 15 2016 -- don't get trigger SF for loose muons (MuonTriggerScaleFactors tool complains... not yet sure if it is a problem
        // from our mangled setup or the tool's issue)
        if(m_run_oneST) {
            out.muoTrigSF[MuonId::Medium] = m_susyObj[m_eleIDDefault]->GetTotalMuonTriggerSF(*sf_muon, static_cast<string>(trig_exp_med.Data()));
        }
        else {
            out.muoTrigSF[MuonId::Medium] = m_susyObj[SusyObjId::muoMedium]->GetTotalMuonTriggerSF(*sf_muon, static_cast<string>(trig_exp_med.Data()));
        }

        delete sf_muon;
        delete sf_muon_aux;
        //delete sfMu;

        /////////////////////////////////////////
        //  new method for trigger SF
        /////////////////////////////////////////
        //2485             const vector<string> di_triggers = TriggerTools::ele_muo_triggers();
        //2486             const vector<string> triggers = TriggerTools::getTrigNames();
        const vector<string> triggers = TriggerTools::getTrigNames();
        vector<string> muon_triggers;
        vector<string> single_muo = TriggerTools::single_muo_triggers();
        vector<string> di_muo = TriggerTools::di_muo_triggers();
        muon_triggers.insert(muon_triggers.end(), single_muo.begin(), single_muo.end());
        muon_triggers.insert(muon_triggers.end(), di_muo.begin(), di_muo.end());

        for(int itrig = 0; itrig < (int)triggers.size(); itrig++) {
            string chain_name = triggers.at(itrig);
            if(std::find(muon_triggers.begin(), muon_triggers.end(), chain_name) == muon_triggers.end()) continue;
            float eff_data = m_susyObj[SusyObjId::muoMedium]->GetMuonTriggerEfficiency(in, chain_name, true); 
            float eff_mc   = m_susyObj[SusyObjId::muoMedium]->GetMuonTriggerEfficiency(in, chain_name, false); 
            out.muoTrigEffData_medium[itrig] = eff_data;
            out.muoTrigEffMC_medium[itrig] = eff_mc;

            eff_data = m_susyObj[SusyObjId::muoLoose]->GetMuonTriggerEfficiency(in, chain_name, true);
            eff_mc = m_susyObj[SusyObjId::muoLoose]->GetMuonTriggerEfficiency(in, chain_name, false);
            out.muoTrigEffData_loose[itrig] = eff_data;
            out.muoTrigEffMC_loose[itrig] = eff_mc;
            //if(chain_name == "HLT_mu24")
            //    cout << "SusyNtMaker::store_muon    NOMINAL MUON SF pt = " << out.Pt() << "   eff_data = " << eff_data << "  eff_mc = " << eff_mc << endl;
        } // itrig
    }

    //////////////////////////////////////
    // Systematic Varation of SF
    //////////////////////////////////////
    if(mc() && sys() && fabs(out.eta)<2.5 && out.pt>20){
        for(const auto& sysInfo : systInfoList) {
            if(!(sysInfo.affectsType == ST::SystObjType::Muon && sysInfo.affectsWeights)) continue;
            const CP::SystematicSet& sys = sysInfo.systset;
            SusyNtSys ourSys = CPsys2sys((sys.name()).c_str());
            for(int i : Susy::muonIds()){
                int index_to_check = (m_run_oneST==true ? (int)m_eleIDDefault : i);
                if(m_susyObj[index_to_check]->applySystematicVariation(sys) != CP::SystematicCode::Ok) {
                    cout << "SusyNtMaker::store_muon    WARNING Cannot configure SUSYTools for systematic " << sys.name() << endl;
                    continue;
                }
                if(m_run_oneST) break;
            }
            vector<float> sf;
            sf.assign(MuonId::MuonIdInvalid, 1);
            if(m_run_oneST) {
                sf[MuonId::Medium] = m_susyObj[m_eleIDDefault]->GetSignalMuonSF(in, recoSF, isoSF);
                sf[MuonId::Loose] =  m_susyObj[m_eleIDDefault]->GetSignalMuonSF(in, recoSF, isoSF);
            }
            else {
                sf[MuonId::Medium] = m_susyObj[SusyObjId::muoMedium]->GetSignalMuonSF(in, recoSF, isoSF);
                sf[MuonId::Loose] = m_susyObj[SusyObjId::muoLoose]->GetSignalMuonSF(in, recoSF, isoSF);
            }

            for(int i=MuonId::VeryLoose; i<MuonId::MuonIdInvalid; i++){
                if     (ourSys == NtSys::MUON_EFF_STAT_UP)        out.errEffSF_stat_up[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_EFF_STAT_DN)        out.errEffSF_stat_dn[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_EFF_SYS_UP)         out.errEffSF_syst_up[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_EFF_SYS_DN)         out.errEffSF_syst_dn[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_EFF_STAT_LOWPT_UP)  out.errEffSF_stat_lowpt_up[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_EFF_STAT_LOWPT_DN)  out.errEffSF_stat_lowpt_dn[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_EFF_SYS_LOWPT_UP)   out.errEffSF_syst_lowpt_up[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_EFF_SYS_LOWPT_DN)   out.errEffSF_syst_lowpt_dn[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_ISO_STAT_UP)        out.errIso_stat_up[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_ISO_STAT_DN)        out.errIso_stat_dn[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_ISO_SYS_UP)         out.errIso_syst_up[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_ISO_SYS_DN)         out.errIso_syst_dn[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_TTVA_STAT_UP)       out.errTTVA_stat_up[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_TTVA_STAT_DN)       out.errTTVA_stat_dn[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_TTVA_SYS_UP)        out.errTTVA_syst_up[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_TTVA_SYS_DN)        out.errTTVA_syst_dn[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_BADMUON_STAT_UP)    out.errBadMu_stat_up[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_BADMUON_STAT_DN)    out.errBadMu_stat_dn[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_BADMUON_SYS_UP)     out.errBadMu_syst_up[i] = sf[i] - out.muoEffSF[i];
                else if(ourSys == NtSys::MUON_BADMUON_SYS_DN)     out.errBadMu_syst_dn[i] = sf[i] - out.muoEffSF[i];

/*
                if(i==1 || i==2) {
                if(i==1)
                    cout << "MuonId: Loose " << endl;
                else if(i==2)
                    cout << "MuonId: Medium" << endl;
                cout << "    effstat        : " << out.errEffSF_stat_up[i] << "  " << out.errEffSF_stat_dn[i] << endl;
                cout << "    effsyst        : " << out.errEffSF_syst_up[i] << "  " << out.errEffSF_syst_dn[i] << endl;
                cout << "    eff_stat_lowpt : " << out.errEffSF_stat_lowpt_up[i] << "  " << out.errEffSF_stat_lowpt_dn[i] << endl;
                cout << "    eff_syst_lowpt : " << out.errEffSF_syst_lowpt_up[i] << "  " << out.errEffSF_syst_lowpt_dn[i] << endl;
                cout << "    eff_iso_stat   : " << out.errIso_stat_up[i] << "  " << out.errIso_stat_dn[i] << endl;
                cout << "    eff_iso_syst   : " << out.errIso_syst_up[i] << "  " << out.errIso_syst_dn[i] << endl;
                cout << "    eff_ttva_stat  : " << out.errTTVA_stat_up[i]<< "  " << out.errTTVA_stat_dn[i] << endl;
                cout << "    eff_ttva_syst  : " << out.errTTVA_syst_up[i]<< "  " << out.errTTVA_syst_dn[i] << endl;
                cout << "    eff_bad_mu_sys : " << out.errBadMu_syst_up[i]<<"  " << out.errBadMu_syst_dn[i] << endl;
                cout << "    eff_bad_mu_stat: " << out.errBadMu_stat_up[i]<<"  " << out.errBadMu_stat_dn[i] << endl;
                }
*/
            }
        } // sysInfo
        for(int i : Susy::muonIds()){
            int index_to_check = (m_run_oneST==true ? (int)m_eleIDDefault : i);
            if(m_susyObj[index_to_check]->resetSystematics() != CP::SystematicCode::Ok){
                cout << "SusyNtMaker::store_muon    ERROR Cannot reset SUSYTools systematics. Aborting." << endl;
                abort();
            }
            if(m_run_oneST) break;
        }
    } // ifMC && sys
    else {
        for(int i=MuonId::VeryLoose; i<MuonId::MuonIdInvalid; i++){
            out.errEffSF_stat_up[i] = out.errEffSF_stat_dn[i] = 0;
            out.errEffSF_syst_up[i] = out.errEffSF_syst_dn[i] = 0;
            out.errEffSF_stat_lowpt_up[i] = out.errEffSF_stat_lowpt_dn[i] = 0;
            out.errEffSF_syst_lowpt_up[i] = out.errEffSF_syst_lowpt_dn[i] = 0;
            out.errIso_stat_up[i] = out.errIso_stat_dn[i] = 0;
            out.errIso_syst_up[i] = out.errIso_syst_dn[i] = 0;
            out.errTTVA_stat_up[i] = out.errTTVA_stat_dn[i] = 0;
            out.errTTVA_syst_up[i] = out.errTTVA_syst_dn[i] = 0;
            out.errBadMu_stat_up[i] = out.errBadMu_stat_dn[i] = 0;
            out.errBadMu_syst_up[i] = out.errBadMu_syst_dn[i] = 0;
        }
    }

    if(mc() && sys() && fabs(out.eta)<2.5 && out.pt>20){
        for(const auto& sysInfo : systInfoList) {
            if(!(sysInfo.affectsType == ST::SystObjType::Muon && sysInfo.affectsWeights)) continue;
            const CP::SystematicSet& sys = sysInfo.systset;
            SusyNtSys ourSys = CPsys2sys((sys.name()).c_str());

            // only want trigger SF variations
            if(!( (ourSys == NtSys::MUON_EFF_TRIG_STAT_UP) || (ourSys == NtSys::MUON_EFF_TRIG_STAT_DN)
                || (ourSys == NtSys::MUON_EFF_TRIG_SYST_UP) || (ourSys == NtSys::MUON_EFF_TRIG_SYST_DN) )) continue;

            // only store muon trigger SF variations for medium muons
            if(m_susyObj[SusyObjId::muoMedium]->applySystematicVariation(sys) != CP::SystematicCode::Ok) {
                cout << "SusyNtMaker::store_muon    WARNING Cannot configure SUSYTools for systematic " << sys.name() << endl;
                continue;
            }

            const vector<string> triggers = TriggerTools::getTrigNames();
            vector<string> muon_triggers;
            vector<string> single_muo = TriggerTools::single_muo_triggers();
            vector<string> di_muo = TriggerTools::di_muo_triggers();
            muon_triggers.insert(muon_triggers.end(), single_muo.begin(), single_muo.end());
            muon_triggers.insert(muon_triggers.end(), di_muo.begin(), di_muo.end());

            // loop over the global trigger indices
            for(int itrig = 0; itrig < (int)triggers.size(); itrig++) {
                string chain_name = triggers.at(itrig);
                if(std::find(muon_triggers.begin(), muon_triggers.end(), chain_name) == muon_triggers.end()) continue;
                float eff_data_sys = m_susyObj[SusyObjId::muoMedium]->GetMuonTriggerEfficiency(in, chain_name, true);
                float eff_mc_sys = m_susyObj[SusyObjId::muoMedium]->GetMuonTriggerEfficiency(in, chain_name, false);
            

                if(ourSys == NtSys::MUON_EFF_TRIG_STAT_UP) {
                    //if(chain_name == "HLT_mu24")
                    //    cout << "SusyNtMaker::store_muon      > TRIG_STAT_UP   eff_data = " << eff_data_sys << "   eff_mc = " << eff_mc_sys << endl;
                    out.muoTrigEffErrData_stat_up_medium[itrig] = eff_data_sys;
                    out.muoTrigEffErrMC_stat_up_medium[itrig] = eff_mc_sys;
                }
                else if(ourSys == NtSys::MUON_EFF_TRIG_STAT_DN) {
                    //if(chain_name == "HLT_mu24")
                    //    cout << "SusyNtMaker::store_muon      > TRIG_STAT_DN   eff_data = " << eff_data_sys << "   eff_mc = " << eff_mc_sys << endl;
                    out.muoTrigEffErrData_stat_dn_medium[itrig] = eff_data_sys;
                    out.muoTrigEffErrMC_stat_dn_medium[itrig] = eff_mc_sys;
                }
                else if(ourSys == NtSys::MUON_EFF_TRIG_SYST_UP) {
                    //if(chain_name == "HLT_mu24")
                    //    cout << "SusyNtMaker::store_muon      > TRIG_SYST_UP   eff_data = " << eff_data_sys << "   eff_mc = " << eff_mc_sys << endl;
                    out.muoTrigEffErrData_syst_up_medium[itrig] = eff_data_sys;
                    out.muoTrigEffErrMC_syst_up_medium[itrig] = eff_mc_sys;
                }
                else if(ourSys == NtSys::MUON_EFF_TRIG_SYST_DN) {
                    //if(chain_name == "HLT_mu24")
                    //    cout << "SusyNtMaker::store_muon      > TRIG_SYST_DN   eff_data = " << eff_data_sys << "   eff_mc = " << eff_mc_sys << endl;
                    out.muoTrigEffErrData_syst_dn_medium[itrig] = eff_data_sys;
                    out.muoTrigEffErrMC_syst_dn_medium[itrig] = eff_mc_sys;
                }
            } // itrig

            if(m_susyObj[SusyObjId::muoMedium]->resetSystematics() != CP::SystematicCode::Ok) {
                cout << "SusyNtMaker::store_muon    ERROR Cannot reset SUSYTools systematics. Aborting." << endl;
                abort();
            }
        }
    } // ifMC && sys

    //______________ ALL DONE WITH THE MUON ______________ //
    out.idx = (m_susyNt.muo()->size());
    m_susyNt.muo()->push_back(out);

}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::store_jet(const xAOD::Jet& in)
{
    if(dbg()>=15) cout << "SusyNtMaker::store_jet    Storing Jet (pt=" << in.pt()*MeV2GeV << ")" << endl;

    const xAOD::EventInfo* eventinfo = xaodEventInfo();

    Susy::Jet out;

    ///////////////////////////////////////////
    // 4-vector
    ///////////////////////////////////////////
    double pt = ( (in.pt()*MeV2GeV < 0) ? 0 : in.pt()*MeV2GeV);
    double m = ( (in.m()*MeV2GeV < 0) ? 0 : in.m()*MeV2GeV);
    double eta(in.eta()), phi(in.phi());
    out.SetPtEtaPhiM(pt, eta, phi, m);
    out.pt  = pt;
    out.eta = eta;
    out.phi = phi;
    out.m   = m;
    
    ///////////////////////////////////////////
    // Associated tracks
    ///////////////////////////////////////////
    vector<int> nTrkVec;
    in.getAttribute(xAOD::JetAttribute::NumTrkPt500, nTrkVec);
    int jet_nTrk = (m_susyObj[m_eleIDDefault]->GetPrimVtx()==0 || nTrkVec.size()==0) ? 0 : nTrkVec[m_susyObj[m_eleIDDefault]->GetPrimVtx()->index()];
    out.nTracks = jet_nTrk;
    float jet_sumTrkPt = (m_susyObj[m_eleIDDefault]->GetPrimVtx()==0 ? 0 :
                            in.getAttribute< std::vector<float> >(xAOD::JetAttribute::SumPtTrkPt500)[m_susyObj[m_eleIDDefault]->GetPrimVtx()->index()]);
    out.sumTrkPt = jet_sumTrkPt * MeV2GeV;

    ///////////////////////////////////////////
    // JVT (Jet Vertex Tagger)
    ///////////////////////////////////////////
    static SG::AuxElement::Accessor<float> acc_jvt("Jvt");
    out.jvt = acc_jvt(in);

    ///////////////////////////////////////////
    // Truth Labeling
    ///////////////////////////////////////////
    if(mc()) {
        in.getAttribute("HadronConeExclTruthLabelID", out.truthLabel);
    }

    ///////////////////////////////////////////
    // b-tagging
    ///////////////////////////////////////////
    double weight_mv2c20(0.);
    if(!in.btagging()->MVx_discriminant("MV2c20", weight_mv2c20)) {
        cout << "SusyNtMaker::store_jet    WARNING Failed to retrieve MV2c20 weight for jet (event:" << eventinfo->eventNumber() << ", pt=" << in.pt()*MeV2GeV << ")" << endl;
    }
    out.mv2c20 = weight_mv2c20;

    double weight_mv2c10(0.);
    if(!in.btagging()->MVx_discriminant("MV2c10", weight_mv2c10)) {
        cout << "SusyNtMaker::store_jet    WARNING Failed to retrieve MV2c10 weight for jet (event: " << eventinfo->eventNumber() << ", pt=" << in.pt()*MeV2GeV << ")" << endl;
    }
    out.mv2c10 = weight_mv2c10;

    // SUSYTools' label -- not to be used at the analysis level
    out.bjet = m_susyObj[m_eleIDDefault]->IsBJet(in) ? 1 : 0;

    //////////////////////////////////////////
    // b-tagging scale factor
    //////////////////////////////////////////

    // dantrim June 13 2017 -- TODO check implementation/storage of these
    out.effscalefact = mc() ? in.auxdata<double>("effscalefact") : 1;

    
    //////////////////////////////////////////
    // flavor tagging systematics
    //////////////////////////////////////////
    if(mc() && sys()) {
        for(const auto& sysInfo : systInfoList) {
            if(!(sysInfo.affectsType == ST::SystObjType::BTag && sysInfo.affectsWeights)) continue;
            const CP::SystematicSet& sys = sysInfo.systset;
            if(m_susyObj[m_eleIDDefault]->applySystematicVariation(sys) != CP::SystematicCode::Ok) {
                cout << "SusyNtMaker::store_jet    WARNING Cannot configure SUSYTools for systematic " << sys.name() << endl;
                continue;
            }
            SusyNtSys ourSys = CPsys2sys((sys.name()).c_str());
            if(!(ourSys==NtSys::FT_EFF_B_systematics_UP       || ourSys==NtSys::FT_EFF_B_systematics_DN
               ||ourSys==NtSys::FT_EFF_C_systematics_UP       || ourSys==NtSys::FT_EFF_C_systematics_DN
               ||ourSys==NtSys::FT_EFF_Light_systematics_UP   || ourSys==NtSys::FT_EFF_Light_systematics_DN
               ||ourSys==NtSys::FT_EFF_extrapolation_UP       || ourSys==NtSys::FT_EFF_extrapolation_DN
               ||ourSys==NtSys::FT_EFF_extrapolation_charm_UP || ourSys==NtSys::FT_EFF_extrapolation_charm_DN) ) continue;

            // redecorate the jets
            m_susyObj[m_eleIDDefault]->BtagSF(XaodAnalysis::xaodJets(systInfoList[0]));

            double sf_out = out.effscalefact - in.auxdata< double >("effscalefact");
            out.setFTSys(ourSys, sf_out);
        } // sysInfo
        // reset systematics
        if(m_susyObj[m_eleIDDefault]->resetSystematics() != CP::SystematicCode::Ok){
            cout << "SusyNtMaker::store_jet    Cannot reset SUSYTools systematics. Aborting." << endl;
            abort();
        }
    } // isMC && sys

    //////////////////////////////////////////
    // JVT SF systematic variation
    //////////////////////////////////////////
    if(mc()) { 
        if(out.Pt()>20. && out.Pt()<60. && fabs(out.Eta()) < 2.4 && fabs(out.Eta())>0.){ 

            xAOD::JetContainer *jvt_jet = new xAOD::JetContainer;
            xAOD::JetAuxContainer *jvt_jet_aux = new xAOD::JetAuxContainer;
            jvt_jet->setStore(jvt_jet_aux);
            xAOD::Jet* jvtJet = new xAOD::Jet;
            jvtJet->makePrivateStore(in);
            jvt_jet->push_back(jvtJet);

            out.jvtEff = m_susyObj[m_eleIDDefault]->JVT_SF(jvt_jet);
            if(sys()) {
                for(const auto& sysInfo : systInfoList) {
                    if(!(sysInfo.affectsType == ST::SystObjType::Jet && sysInfo.affectsWeights)) continue;
                    const CP::SystematicSet& sys = sysInfo.systset;
                    if(m_susyObj[m_eleIDDefault]->applySystematicVariation(sys) != CP::SystematicCode::Ok) {
                        cout << "SusyNtMaker::store_jet    Cannot configure SUSYTools for systematic " << sys.name () << endl;
                        continue;
                    }
                    SusyNtSys ourSys = CPsys2sys((sys.name()).c_str());
                    if(!(ourSys==NtSys::JET_JVTEff_UP || ourSys==NtSys::JET_JVTEff_DN)) continue;

                    if(ourSys==NtSys::JET_JVTEff_UP) {
                        out.jvtEff_up = out.jvtEff - m_susyObj[m_eleIDDefault]->JVT_SFsys(jvt_jet, sys);
                    }
                    else if(ourSys==NtSys::JET_JVTEff_DN) {
                        out.jvtEff_dn = out.jvtEff - m_susyObj[m_eleIDDefault]->JVT_SFsys(jvt_jet, sys);
                    }
                } // sysInfo
                // reset systematics
                if(m_susyObj[m_eleIDDefault]->resetSystematics() != CP::SystematicCode::Ok){
                    cout << "SusyNtMaker::store_jet    Cannot reset SUSYTools systematics. Aborting." << endl;
                    abort();
                }
            } // if sys
            delete jvt_jet;
            delete jvt_jet_aux;
            //delete jvtJet;
        } // pt/eta range
    }//mc

    //////////////////////////////////////////
    // Jet attributes
    //////////////////////////////////////////
    out.detEta = (in.jetP4(xAOD::JetConstitScaleMomentum)).eta();
    in.getAttribute(xAOD::JetAttribute::EMFrac,out.emfrac);
    in.getAttribute(xAOD::JetAttribute::BchCorrJet,out.bch_corr_jet);
    in.getAttribute(xAOD::JetAttribute::BchCorrCell,out.bch_corr_cell);

    //////////////////////////////////////////
    // Bad Jet ID
    //////////////////////////////////////////
    out.isBadVeryLoose = (bool)in.auxdata<char>("bad") ? 1 : 0;


    //////////////////////////////////////////
    // Hot Tile
    //////////////////////////////////////////
    // dantrim June 13 2017 -- TODO confirm that this variable is non-zero and still used
    float fracSamplingMax, samplingMax;
    in.getAttribute(xAOD::JetAttribute::SamplingMax, samplingMax);
    in.getAttribute(xAOD::JetAttribute::FracSamplingMax, fracSamplingMax);

    //______________ ALL DONE WITH THE JET ______________ //
    out.idx = (m_susyNt.jet()->size());
    m_susyNt.jet()->push_back(out);
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::store_tau(const xAOD::TauJet& in)
{
    if(dbg()>=15) cout << "SusyNtMaker::store_tau   Storing tau (pt=" << in.pt()*MeV2GeV << ")" << endl; 

    Susy::Tau out;

    //////////////////////////////////////////
    // 4-vector
    //////////////////////////////////////////
    double pt(in.pt()*MeV2GeV), eta(in.eta()), phi(in.phi()), m(in.m()*MeV2GeV);

    out.SetPtEtaPhiM(pt, eta, phi, m);
    out.pt  = pt;
    out.eta = eta;
    out.phi = phi;
    out.m   = m;
    out.q = int(in.charge());

    //////////////////////////////////////
    // TauSelectionTool WP
    //////////////////////////////////////
    out.loose   = static_cast<bool>( m_tauSelToolLoose->accept(in) );
    out.medium  = static_cast<bool>( m_tauSelToolMedium->accept(in) );
    out.tight   = static_cast<bool>( m_tauSelToolTight->accept(in) );

    //////////////////////////////////////
    // number of associated tracks
    //////////////////////////////////////
    out.nTrack = in.nTracks();


    //////////////////////////////////////
    // Truth info/classification
    //////////////////////////////////////
    if (mc()){
        auto truthTau = m_tauTruthMatchingTool->getTruth(in);
        out.isTruthMatched = (bool)in.auxdata<char>("IsTruthMatched");
        if(out.isTruthMatched) {
            if(truthTau->isTau()) {
                out.truthNProngs = int(truthTau->auxdata<size_t>("numCharged"));
                out.isHadronicTau = (bool)truthTau->auxdata<char>("IsHadronicTau");
            } //isTau
            out.truthCharge = int(truthTau->charge());
            out.truthPdgId = int(truthTau->absPdgId());

            out.truthType = truthTau->auxdata<unsigned int>("classifierParticleType");
            out.truthOrigin = truthTau->auxdata<unsigned int>("classifierParticleOrigin");

        }//isTruthMatched
        // TODO -- dantrim Feb 23 2016 -- what vars are needed from truth jet
        //auto truthJetLink = in.auxdata< ElementLink< xAOD::JetContainer > >("truthJetLink");
        //if(truthJetLink.isValid()) {
        //    const xAOD::Jet* truthJet = *truthJetLink;
        //    cout << "  > tau was matched to truth jet with (pt, eta, phi, m) = (" << truthJet->p4().Pt() << ", " << truthJet->p4().Eta() << ", " << truthJet->p4().Phi() << ", " << truthJet->p4().M() << endl;
        //}

        // MJF: Scale factors are only valid for pre-selected taus
        if ( in.auxdata< char >("baseline") ) {
            //if (in.nTracks() > 0) m_TauEffEleTool->applyEfficiencyScaleFactor(in);
            // these #'s are dummies for testing the various SF's!
            bool idSF = true;
            bool trigSF = false;
            out.looseEffSF  = m_susyObj[m_eleIDDefault]->GetSignalTauSF(in,idSF,trigSF);
            out.mediumEffSF = m_susyObj[m_eleIDDefault]->GetSignalTauSF(in,idSF,trigSF);
            out.tightEffSF  = m_susyObj[m_eleIDDefault]->GetSignalTauSF(in,idSF,trigSF);
        }

        if(dbg()>=15) cout << "SusyNtMaker::store_tau    MCTruthClassifier found Tau with (truthType, origin) = (" << out.truthType << ", " << out.truthOrigin << ")" << endl;

    }// if isMC

    //______________ ALL DONE WITH THE TAU ______________ //
    m_susyNt.tau()->push_back(out);
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::store_photon(const xAOD::Photon& in)
{

    if(dbg()>=15) cout << "SusyNtMaker::store_photon    Filling photon (pt=" << in.pt()*MeV2GeV << ")" << endl;

    Susy::Photon out;

    //////////////////////////////////////
    // 4-vector
    //////////////////////////////////////
    double pt(in.pt()*MeV2GeV), eta(in.eta()), phi(in.phi()), m(in.m()*MeV2GeV);
    out.SetPtEtaPhiM(pt, eta, phi, m);
    out.pt  = pt;
    out.eta = eta;
    out.phi = phi;
    out.m   = m;

    //////////////////////////////////////
    // Author information
    //////////////////////////////////////
    out.author = static_cast<int>(in.author());
    out.authorPhoton = in.author() & xAOD::EgammaParameters::AuthorPhoton;
    out.authorAmbiguous = in.author() & xAOD::EgammaParameters::AuthorAmbiguous;

    //////////////////////////////////////
    // Photon is Converted
    //////////////////////////////////////
    out.isConv = xAOD::EgammaHelpers::isConvertedPhoton(&in);

    //////////////////////////////////////
    // IsEM ID
    //////////////////////////////////////
    out.loose = (bool)m_photonSelLoose->accept(&in);
    out.tight = (bool)m_photonSelTight->accept(&in);

    //////////////////////////////////////
    // CaloCluster
    //////////////////////////////////////
    const xAOD::CaloCluster* c = in.caloCluster();
    if(c) {
        out.clusE   = c->e()*MeV2GeV;
        // use coordinates from 2nd sampling
        out.clusEtaBE = c->etaBE(2);
        out.clusPhiBE = c->phiBE(2);
        out.clusEta = c->eta();
        out.clusPhi = c->phi();
    }

    //////////////////////////////////////
    // OQ
    //////////////////////////////////////
    out.OQ = in.isGoodOQ(xAOD::EgammaParameters::BADCLUSPHOTON);

    //////////////////////////////////////
    // Isolation Variables
    //////////////////////////////////////
    out.topoEtcone40 = in.isolationValue(xAOD::Iso::topoetcone40) * MeV2GeV;

    //////////////////////////////////////
    // Isolation Selection
    //////////////////////////////////////
    out.isoFixedCutTight         = m_isoToolGradientLooseTight->accept(in) ? true : false;
    out.isoFixedCutTightCaloOnly = m_isoToolGradientTightCalo->accept(in) ? true : false;
    out.isoFixedCutLoose         = m_isoToolLooseTrackOnlyLoose->accept(in) ? true : false;

    //__________________ DONE WITH THE PHOTON _______________ //
    m_susyNt.pho()->push_back(out);
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::run_kinematic_systematics()
{
    if(dbg()>=5) cout << "SusyNtMaker::run_kinematic_systematics    Beginning systematic loop" << endl;

/*
//     useful for figuring out what we have and what we expect
      for(const auto& sysInfo : systInfoList) {
        const CP::SystematicSet& sys = sysInfo.systset;
        if(sys.name()=="") continue;
        SusyNtSys ourSys = CPsys2sys((sys.name()).c_str());
        string affects = "";
        if(sysInfo.affectsType == ST::SystObjType::Jet) affects = "JET";
        else if(sysInfo.affectsType == ST::SystObjType::Muon) affects = "MUON";
        else if(sysInfo.affectsType == ST::SystObjType::Egamma) affects = "EGAMMA";
        else if(sysInfo.affectsType == ST::SystObjType::Electron) affects = "ELECTRON";
        else if(sysInfo.affectsType == ST::SystObjType::Photon) affects = "PHOTON";
        else if(sysInfo.affectsType == ST::SystObjType::Tau) affects = "TAU";
        else if(sysInfo.affectsType == ST::SystObjType::BTag) affects = "BTAG";
        else if(sysInfo.affectsType == ST::SystObjType::MET_TST) affects = "MET_TST";
        else if(sysInfo.affectsType == ST::SystObjType::MET_CST) affects = "MET_CST";
        else if(sysInfo.affectsType == ST::SystObjType::MET_Track) affects = "MET_TRACK";
        else if(sysInfo.affectsType == ST::SystObjType::EventWeight) affects = "EVENTWEIGHT";
        else { affects = "UKNOWN"; }
        string kinOrSys = "";
        if(sysInfo.affectsKinematics) kinOrSys = "Kinematics";
        else if(sysInfo.affectsWeights) kinOrSys = "Weights";
        cout << "systematic: " << (sys.name()).c_str() << "                 ours: " << NtSys::SusyNtSysNames.at(ourSys) << "   affects: " << affects << "  " << kinOrSys << endl;
    }
    cout << endl;
*/

    for(const auto& sysInfo : systInfoList) {
        const CP::SystematicSet& sys = sysInfo.systset;
        if(sys.name()=="") continue; // skip nominal
        if(!sysInfo.affectsKinematics) continue;
        if(dbg()>=15) cout << "SusyNtMaker::run_kinematic_systematics     --------------------------------------------------------------" << endl;
        if(dbg()>=15) cout << "SusyNtMaker::run_kinematic_systematics     > Variation: " << sys.name().c_str() << endl;

        SusyNtSys ourSys = CPsys2sys((sys.name()).c_str());
        if(ourSys == NtSys::SYS_UNKNOWN) continue;

        if(dbg()>=15) cout << "SusyNtMaker::run_kinematic_systematics        >> Matches our systematic: " << NtSys::SusyNtSysNames.at(ourSys) << endl;

        if(m_susyObj[m_eleIDDefault]->applySystematicVariation(sys) != CP::SystematicCode::Ok) {
            cout << "SusyNtMaker::run_kinematic_systematics    WARNING Cannot configure SUYSTools for systematic " << sys.name() << ", will not apply this variation" << endl;
            continue;
        }


        /////////////////////////////////////////////
        // save objects with vthe ariations applied
        ////////////////////////////////////////////

        // these clearing steps ARE NECESSARY TO AVOID MEMORY LEAKS
        clear_output_objects(false);
        delete_shallow_copies(false);

        fill_objects(ourSys, sysInfo);

        // retrieve the MET
        retrieveXaodMet(sysInfo, ourSys);
        retrieveXaodTrackMet(sysInfo, ourSys);

        // electrons
        store_electron_kinematic_sys(sysInfo, ourSys);

        // muons
        store_muon_kinematic_sys(sysInfo, ourSys);

        // taus
        store_tau_kinematic_sys(sysInfo, ourSys);

        // jets
        store_jet_kinematic_sys(sysInfo, ourSys);

        // MET
        fill_met_variables(ourSys);

        // Track MET
        fill_track_met_variables(ourSys);

        // Reset the systematics registry, otherwise the TStore will not be able to load new object collections
        if ( m_susyObj[m_eleIDDefault]->resetSystematics() != CP::SystematicCode::Ok){
            cout << "SusyNtMaker::run_kinematic_systematics    ERROR Cannot reset SUSYTools systematics. Aborting." << endl;
            abort();
        }

    } // sysInfo

}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::store_electron_kinematic_sys(ST::SystInfo sysInfo, SusyNtSys sys)
{
    if(!ST::testAffectsObject(xAOD::Type::Electron, sysInfo.affectsType)) return;

    xAOD::ElectronContainer* electrons     = xaodElectrons(sysInfo,sys);
    xAOD::ElectronContainer* electrons_nom = xaodElectrons(sysInfo,NtSys::NOM);

    if(dbg()>=15) cout << "SusyNtMaker::store_electron_kinematic_sys    " << NtSys::SusyNtSysNames.at(sys) << endl;

    for(const auto &iEl : m_preElectrons) {
        const xAOD::Electron* ele = electrons->at(iEl);
        
        const xAOD::Electron* ele_nom = NULL;
        Susy::Electron* ele_susyNt = NULL;
        int idx_susyNt = -1;
        for(uint idx = 0; idx < m_preElectrons_nom.size(); idx++){
            int iEl_nom = m_preElectrons_nom[idx];
            if(iEl == iEl_nom) {
                ele_nom = electrons_nom->at(iEl_nom);
                ele_susyNt = & m_susyNt.ele()->at(idx);
                idx_susyNt = idx;
                if(dbg()>=15){
                    cout << "SusyNtMaker::store_electron_kinematic_sys    Found matching electron for ele " << idx_susyNt << " (sys=" << SusyNtSysNames.at(sys) << ")   (idx_sys, idx_nom) = (" << iEl << "," << iEl_nom << "), (pT_sys, eta_sys) = (" << ele->pt()*MeV2GeV << "," << ele->eta() << ")  (pT_nom, eta_nom) = (" << ele_nom->pt()*MeV2GeV << "," << ele_nom->eta() << ")" << endl;
                    //ele_susyNt->print();
                }
                if( fabs(ele_nom->eta() - ele->eta())>0.001 || fabs(ele_nom->phi() - ele->phi())>0.001)
                    cout << "SusyNtMaker::store_electron_kinematic_sys    WARNING Index mis-match!" << endl;
                break;
            }
        }
    
        //Nominal electron was not found. Add it at its nominal scale to susyNt and m_preElectron_nom 
        //this happens if the systematic varied object passes the (e.g. pT) threshold that the nominal did not (c.f. XaodAnalysis::fill_baseline_objects)
        if(ele_susyNt == NULL){
            ele_nom = electrons_nom->at(iEl);
            if(dbg()>=20) {
                cout << "SusyNtMaker::store_electron_kinematic_sys    Nominal electron not found at sys idx " << iEl << " (there are " << m_preElectrons_nom.size() << " nominal electrons, " << m_preElectrons.size() << " sys varied electrons) : pT_sys=" << ele->pt()*MeV2GeV << "  pT_nom=" << ele_nom->pt()*MeV2GeV << "  -> Adding (new) nominal electron to output susyNt for SF calculation" << endl;
            }
            store_electron(*ele_nom, iEl);
            m_preElectrons_nom.push_back(iEl);
            ele_susyNt = & m_susyNt.ele()->back(); // now get the newly inserted element and use it
        }

        // now calculate the shift in the electron kinematics
        // store as shift/nom
        float sf = ele->e() / ele_nom->e();
        if(dbg()>=20) cout << "SusyNtMaker::store_electron_kinematic_sys      > (sys="<<SusyNtSysNames.at(sys)<<") electron SF " << sf << endl;
        if     ( sys == NtSys::EG_RESOLUTION_ALL_DN ) ele_susyNt->res_all_dn = sf;
        else if( sys == NtSys::EG_RESOLUTION_ALL_UP ) ele_susyNt->res_all_up = sf;
        else if( sys == NtSys::EG_SCALE_ALL_DN ) ele_susyNt->scale_all_dn = sf;
        else if( sys == NtSys::EG_SCALE_ALL_UP ) ele_susyNt->scale_all_up = sf;
    }
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::store_muon_kinematic_sys(ST::SystInfo sysInfo, SusyNtSys sys)
{
    if(!ST::testAffectsObject(xAOD::Type::Muon, sysInfo.affectsType)) return;

    xAOD::MuonContainer* muons     = xaodMuons(sysInfo,sys);
    xAOD::MuonContainer* muons_nom = xaodMuons(sysInfo, NtSys::NOM);

    if(dbg()>=15) cout << "SusyNtMaker::store_muon_kinematic_sys    " << NtSys::SusyNtSysNames.at(sys) << endl;
    for(const auto& iMu : m_preMuons) {
        const xAOD::Muon* mu = muons->at(iMu);
        
        const xAOD::Muon* mu_nom = NULL;
        Susy::Muon* mu_susyNt = NULL;
        int idx_susyNt = -1;
        for(uint idx = 0; idx < m_preMuons_nom.size(); idx++){
            int iMu_nom = m_preMuons_nom[idx];
            if(iMu == iMu_nom){
                mu_nom = muons_nom->at(iMu_nom);
                mu_susyNt = & m_susyNt.muo()->at(idx);
                idx_susyNt=idx;
                if(dbg()>=15){
                    cout << "SusyNtMaker::store_muon_kinematic_sys    Found matching muon for muo " << idx_susyNt << " (sys=" << SusyNtSysNames.at(sys) << ")   (idx_sys, idx_nom) = (" << iMu << "," << iMu_nom << "), (pT_sys, eta_sys) = (" << mu->pt()*MeV2GeV << "," << mu->eta() << ")  (pT_nom, eta_nom) = (" << mu_nom->pt()*MeV2GeV << "," << mu_nom->eta() << ")" << endl;
                    //mu_susyNt->print();
                }
                if( fabs(mu_nom->eta() - mu->eta())>0.001 || fabs(mu_nom->phi() - mu->phi())>0.001)
                    cout << "SusyNtMaker::store_muon_kinematic_sys    WARNING Index mis-match!" << endl; 
                break;
            }
        }
    
        //Nominal muon was not found. Add it at its nominal scale to susyNt
        //this happens if the systematic varied object passes the (e.g. pT) threshold that the nominal did not (c.f. XaodAnalysis::fill_baseline_objects)
        if(mu_susyNt == NULL){
            mu_nom = muons_nom->at(iMu);
            if(dbg()>=20) {
                cout << "SusyNtMaker::store_muon_kinematic_sys    Nominal muon not found at sys idx " << iMu << " (there are " << m_preMuons_nom.size() << " nominal muons, " << m_preMuons.size() << " sys varied muons) : pT_sys=" << mu->pt()*MeV2GeV << "  pT_nom=" << mu_nom->pt()*MeV2GeV << "  -> Adding (new) nominal muon to output susyNt for SF calculation" << endl;
            }
            store_muon(*mu_nom, *muons);
            m_preMuons_nom.push_back(iMu);
            mu_susyNt = & m_susyNt.muo()->back(); // now get the newly inserted element and use it for SF calculation
        }

        float sf = mu->e() / mu_nom->e();
        if(dbg()>=20) cout << "SusyNtMaker::store_muon_kinematic_sys      > (sys="<< SusyNtSysNames.at(sys)<<") muon SF " << sf << endl;
        if(sys == NtSys::MUON_MS_UP)      mu_susyNt->ms_up = sf;
        else if(sys == NtSys::MUON_MS_DN) mu_susyNt->ms_dn = sf;
        else if(sys == NtSys::MUON_ID_UP) mu_susyNt->id_up = sf;
        else if(sys == NtSys::MUON_ID_DN) mu_susyNt->id_dn = sf;
        else if(sys == NtSys::MUON_SCALE_UP) mu_susyNt->scale_up = sf;
        else if(sys == NtSys::MUON_SCALE_DN) mu_susyNt->scale_dn = sf;
        else if(sys == NtSys::MUON_SAGITTA_RESBIAS_UP) mu_susyNt->sagitta_bias_up = sf;
        else if(sys == NtSys::MUON_SAGITTA_RESBIAS_DN) mu_susyNt->sagitta_bias_dn = sf;
        else if(sys == NtSys::MUON_SAGITTA_RHO_UP) mu_susyNt->sagitta_rho_up = sf;
        else if(sys == NtSys::MUON_SAGITTA_RHO_DN) mu_susyNt->sagitta_rho_dn = sf;
    }
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::store_tau_kinematic_sys(ST::SystInfo sysInfo, SusyNtSys sys)
{
    if(!ST::testAffectsObject(xAOD::Type::Tau, sysInfo.affectsType)) return;

    xAOD::TauJetContainer* taus     = xaodTaus(sysInfo,sys);
    xAOD::TauJetContainer* taus_nom = xaodTaus(sysInfo,NtSys::NOM);

    if(dbg()>=15) cout << "SusyNtMaker::store_tau_kinematic_sys    " << NtSys::SusyNtSysNames.at(sys) << endl;

    vector<int>& saveTaus = cont_taus() ? m_contTaus : m_preTaus;
    vector<int>& saveTaus_nom = cont_taus() ? m_contTaus_nom : m_preTaus_nom;

    for(const auto &iTau : saveTaus) {
        const xAOD::TauJet* tau = taus->at(iTau);
    
        const xAOD::TauJet* tau_nom = NULL;
        Susy::Tau* tau_susyNt = NULL;
        int idx_susyNt = -1;
        for(uint idx = 0; idx < saveTaus_nom.size(); idx++){
            int iTau_nom = saveTaus_nom[idx];
            if(iTau == iTau_nom){
                tau_nom = taus_nom->at(iTau_nom);
                tau_susyNt = & m_susyNt.tau()->at(idx);
                idx_susyNt = idx;
                if(dbg()>=15) {
                    cout << "SusyNtMaker::store_tau_kinematic_sys    Found matching tau for tao " << idx_susyNt << "  (sys=" << SusyNtSysNames.at(sys) << ")    (idx_sys, idx_nom) = (" << iTau << "," << iTau_nom << "), (pT_sys, eta_sys) = (" << tau->pt()*MeV2GeV << "," << tau->eta() << ")  (pT_nom, eta_nom) = (" << tau_nom->pt()*MeV2GeV << "," << tau_nom->eta() << ")" << endl;
                }
                if( fabs(tau_nom->eta() - tau->eta())>0.001 || fabs(tau_nom->phi() - tau->phi())>0.001)
                    cout << "SusyNtMaker::store_tau_kinematic_sys    WARNING Index mis-match!" << endl;
                break;
            }
        }

        //Tau was not found. Add it at its nominal scale to susyNt and m_preTau_nom 
        //this happens if the systematic varied object passes the (e.g. pT) threshold that the nominal did not (c.f. XaodAnalysis::fill_baseline_objects)
        if(tau_susyNt == NULL){
            tau_nom = taus_nom->at(iTau);
            if(dbg()>=20) {
                cout << "SusyNtMaker::store_tau_kinematic_sys    Nominal tau not found at sys idx " << iTau << " (there are " << m_preTaus_nom.size() << " nominal taus, " << m_preTaus.size() << " sys varied taus) : pT_sys=" << tau->pt()*MeV2GeV << "  pT_nom=" << tau_nom->pt()*MeV2GeV << "  -> Adding (new) nominal tau to output susyNt for SF calculation" << endl;
            }
            store_tau(*tau_nom);
            saveTaus_nom.push_back(iTau);
            tau_susyNt = & m_susyNt.tau()->back(); //get the newly inserted taument
        }

        //Calculate systematic SF: shift/nom
        float sf = tau->e() / tau_nom->e();
        if(dbg()>=20) cout << "SusyNtMaker::store_tau_kinematic_sys      > (sys=" << SusyNtSysNames.at(sys) << ") tau SF " << sf << endl;

        if(sys == NtSys::TAUS_SME_TOTAL_UP) tau_susyNt->sme_total_up = sf;
        if(sys == NtSys::TAUS_SME_TOTAL_DN) tau_susyNt->sme_total_dn = sf;
    }
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::store_jet_kinematic_sys(ST::SystInfo sysInfo, SusyNtSys sys)
{
    if(!ST::testAffectsObject(xAOD::Type::Jet, sysInfo.affectsType)) return;

    xAOD::JetContainer* jets     = xaodJets(sysInfo,sys);
    xAOD::JetContainer* jets_nom = xaodJets(sysInfo,NtSys::NOM);
  
    if(dbg()>=15) cout << "SusyNtMaker::store_jet_kinematic_sys    " << NtSys::SusyNtSysNames.at(sys) << endl;

    for(const auto &iJ : m_preJets) {
        const xAOD::Jet* jet = jets->at(iJ);

        const xAOD::Jet* jet_nom = NULL;
        Susy::Jet* jet_susyNt = NULL;
        int idx_susyNt = -1;
        for(uint idx = 0; idx < m_preJets_nom.size(); idx++){
            int iJ_nom = m_preJets_nom[idx];
            if(iJ == iJ_nom){
                jet_nom = jets_nom->at(iJ_nom);
                jet_susyNt = & m_susyNt.jet()->at(idx);
                idx_susyNt = idx;
                if(dbg()>=15) {
                    cout << "SusyNtMaker::store_jet_kinematic_sys    Found matching jet for jet " << idx_susyNt << " (sys=" << SusyNtSysNames.at(sys) << ")  (idx_sys, idx_nom) = (" << iJ <<"," << iJ_nom << "), (pT_sys, eta_sys) = (" << jet->pt()*MeV2GeV <<"," << jet->eta() << ")  (pT_nom, eta_nom) = (" << jet_nom->pt()*MeV2GeV << "," << jet_nom->eta() << ")" << endl;
                }
                if( fabs(jet_nom->eta() - jet->eta())>0.001 || fabs(jet_nom->phi() - jet->phi())>0.001)
                    cout << "SusyNtMaker::store_jet_kinematic_sys    WARNING Index mis-match!" << endl;
                break;
            }
        }

        //Jet was not found. Add it at its nominal scale to susyNt and m_preJet_nom 
        //this happens if the systematic varied object passes the (e.g. pT) threshold that the nominal did not (c.f. XaodAnalysis::fill_baseline_objects)
        if(jet_susyNt == NULL){
            jet_nom = jets_nom->at(iJ);
            if(dbg()>=20) {
                cout << "SusyNtMaker::store_jet_kinematic_sys    Nominal jet not found at sys idx " << iJ << " (there are " << m_preJets_nom.size() << " nominal jets, " << m_preJets.size() << " sys varied jets) : pT_sys=" << jet->pt()*MeV2GeV << "  pT_nom=" << jet_nom->pt()*MeV2GeV << "  -> Adding (new) nominal jet to output susyNt for SF calculation" << endl;
            }
            store_jet(*jet_nom);
            m_preJets_nom.push_back(iJ);
            jet_susyNt = & m_susyNt.jet()->back(); //get the newly inserted jet
        }

        //Calculate systematic SF: shift/nom
        float sf = jet->e() / jet_nom->e();
        if(dbg()>=20) cout << "SusyNtMaker::store_jet_kinematic_sys      > (sys=" << SusyNtSysNames.at(sys) <<") jet SF " << sf << endl;

        if     ( sys == NtSys::JER)                jet_susyNt->jer = sf;
        else if( sys == NtSys::JET_GroupedNP_1_UP) jet_susyNt->groupedNP[0] = sf;
        else if( sys == NtSys::JET_GroupedNP_1_DN) jet_susyNt->groupedNP[1] = sf;
        else if( sys == NtSys::JET_GroupedNP_2_UP) jet_susyNt->groupedNP[2] = sf;
        else if( sys == NtSys::JET_GroupedNP_2_DN) jet_susyNt->groupedNP[3] = sf;
        else if( sys == NtSys::JET_GroupedNP_3_UP) jet_susyNt->groupedNP[4] = sf;
        else if( sys == NtSys::JET_GroupedNP_3_DN) jet_susyNt->groupedNP[5] = sf;
        else if( sys == NtSys::JET_EtaIntercalibration_UP) jet_susyNt->eta_intercal_up = sf;
        else if( sys == NtSys::JET_EtaIntercalibration_DN) jet_susyNt->eta_intercal_dn = sf;
    }
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::perform_dilepton_trigger_matching()
{
    bool verbose_trig_match = false;

    if(dbg()>=10) cout << "SusyNtMaker::perform_dilepton_trigger_matching" << endl;
    const xAOD::EventInfo* ei = xaodEventInfo();

    if(verbose_trig_match) {
        cout << "----------------------------------------------------" << endl;
        cout << " trig evt: " << ei->eventNumber() << endl;
    }

    // get nominal objects
    xAOD::ElectronContainer* electrons = xaodElectrons(systInfoList[0]);
    xAOD::MuonContainer* muons = xaodMuons(systInfoList[0]);

    ///////////////////////////////////////////////////////////////////////////
    // E+E Dilepton Trigger Matching
    ///////////////////////////////////////////////////////////////////////////
    size_t n_stored_ele = m_susyNt.ele()->size();
    for(unsigned int ie_store = 0; ie_store < n_stored_ele; ie_store++) {
        for(unsigned int je_store = 0; je_store < n_stored_ele; je_store++) {
            if(!(je_store > ie_store)) continue;
            int idx_i = m_susyNt.ele()->at(ie_store).idx;
            int idx_j = m_susyNt.ele()->at(je_store).idx;

            // only consider leptons with pT > 17
            bool pass_pt = (m_susyNt.ele()->at(ie_store).pt > TriggerTools::ele_match_pt());
            pass_pt = (pass_pt && (m_susyNt.ele()->at(je_store).pt > TriggerTools::ele_match_pt()));
            if(!pass_pt) continue;

            // this assumes that in SusyNtMaker::store_electron we do not do
            // not apply any further selection
            xAOD::Electron* ele_i = electrons->at(m_preElectrons.at(ie_store));
            xAOD::Electron* ele_j = electrons->at(m_preElectrons.at(je_store));

            const vector<string> di_triggers = TriggerTools::di_ele_triggers();
            const vector<string> triggers = TriggerTools::getTrigNames();
            // loop over entire trigger list in order to use global indices for the triggers
            for(int i = 0; i < (int) triggers.size(); i++) {

                // only test dilepton ee triggers
                if((std::find(di_triggers.begin(), di_triggers.end(), triggers.at(i))==di_triggers.end())) continue;

                bool is_match = dilepton_trigger_matched(ele_i, ele_j, triggers.at(i));

                // build the mask
                DileptonTrigTuple tuple = 0;
                tuple |= (i << 8);
                tuple |= (idx_i << 4);
                tuple |= (idx_j);
                m_susyNt.evt()->m_dilepton_trigger_matches[tuple] = (is_match ? 1 : 0);

                if(is_match && verbose_trig_match) {
                    cout << "SusyNtMaker::perform_dilepton_trigger_matching   EE MATCH [" << triggers.at(i) << "]"
                            << "("<< ie_store << ","<<je_store<<")   matched? " << is_match << "  "
                            << "  ele(" << ie_store << ") pt=" << m_susyNt.ele()->at(ie_store).pt
                            << ", eta=" << m_susyNt.ele()->at(ie_store).eta
                            << "    ele(" << je_store << ") pt=" << m_susyNt.ele()->at(je_store).pt
                            << ", eta=" << m_susyNt.ele()->at(je_store).eta << endl;
                }
            }
        } // je_store
    } // ie_store

    ///////////////////////////////////////////////////////////////////////////
    // M+M Dilepton Trigger Matching
    ///////////////////////////////////////////////////////////////////////////
    size_t n_stored_muo = m_susyNt.muo()->size();
    for(unsigned int im_store = 0; im_store < n_stored_muo; im_store++) {
        for(unsigned int jm_store = 0; jm_store < n_stored_muo; jm_store++) {
            if(!(jm_store > im_store)) continue;

            int idx_i = m_susyNt.muo()->at(im_store).idx;
            int idx_j = m_susyNt.muo()->at(jm_store).idx;

            bool pass_pt = (m_susyNt.muo()->at(im_store).pt > TriggerTools::muo_match_pt());
            pass_pt = (pass_pt && (m_susyNt.muo()->at(jm_store).pt > TriggerTools::muo_match_pt()));
            if(!pass_pt) continue;

            xAOD::Muon* muo_i = muons->at(m_preMuons.at(im_store));
            xAOD::Muon* muo_j = muons->at(m_preMuons.at(jm_store));

            const vector<string> di_triggers = TriggerTools::di_muo_triggers();
            const vector<string> triggers = TriggerTools::getTrigNames();
            // loop over entire trigger list in order to use global indices for the triggers
            for(int i = 0; i < (int) triggers.size(); i++) {

                // only test dilepton mm triggers
                if((std::find(di_triggers.begin(), di_triggers.end(), triggers.at(i))==di_triggers.end())) continue;

                bool is_match = dilepton_trigger_matched(muo_i, muo_j, triggers.at(i));

                // build the mask
                DileptonTrigTuple tuple = 0;
                tuple |= (i << 8);
                tuple |= (idx_i << 4);
                tuple |= (idx_j);
                m_susyNt.evt()->m_dilepton_trigger_matches[tuple] = (is_match ? 1 : 0);

                if(is_match && verbose_trig_match) {
                cout << "SusyNtMaker::perform_dilepton_trigger_matching   MM MATCH [" << triggers.at(i)<< "] "
                        << "("<< im_store << ","<<jm_store<<")   matched? " << is_match << "  "
                        << "  muo(" << im_store << ") pt=" << m_susyNt.muo()->at(im_store).pt
                        << ", eta=" << m_susyNt.muo()->at(im_store).eta
                        << "    muo(" << jm_store << ") pt=" << m_susyNt.muo()->at(jm_store).pt
                        << ", eta=" << m_susyNt.muo()->at(jm_store).eta << endl;
                }
            }
        } // jm_store
    } // im-store

    ///////////////////////////////////////////////////////////////////////////
    // E+M Dilepton Trigger Matching
    ///////////////////////////////////////////////////////////////////////////
    for(unsigned int ie_store = 0; ie_store < n_stored_ele; ie_store++) {
        for(unsigned int im_store = 0; im_store < n_stored_muo; im_store++) {

            int idx_e = m_susyNt.ele()->at(ie_store).idx;
            int idx_m = m_susyNt.muo()->at(im_store).idx;

            bool pass_pt = (m_susyNt.ele()->at(ie_store).pt > TriggerTools::ele_match_pt());
            pass_pt = (pass_pt && (m_susyNt.muo()->at(im_store).pt > TriggerTools::muo_match_pt()));

            if(!pass_pt) continue;

            xAOD::Electron* ele_i = electrons->at(m_preElectrons.at(ie_store));
            xAOD::Muon* muo_i = muons->at(m_preMuons.at(im_store));

            const vector<string> di_triggers = TriggerTools::ele_muo_triggers();
            const vector<string> triggers = TriggerTools::getTrigNames();
            // loop over entire trigger list in order to use global indices for the triggers
            for(int i = 0; i < (int) triggers.size(); i++) {

                // only test dilepton em triggers
                if(std::find(di_triggers.begin(), di_triggers.end(), triggers.at(i))==di_triggers.end()) continue;

                bool is_match = dilepton_trigger_matched(ele_i, muo_i, triggers.at(i));

                // build the mask
                DileptonTrigTuple tuple = 0;
                tuple |= (i << 8);
                tuple |= (idx_e << 4);
                tuple |= (idx_m);
                m_susyNt.evt()->m_dilepton_trigger_matches[tuple] = (is_match ? 1 : 0);

                if(is_match && verbose_trig_match) {
                cout << "SusyNtMaker::perform_dilepton_trigger_matching   EM MATCH [" << triggers.at(i)<< "] "
                        << "("<< ie_store << ","<<im_store<<")   matched? " << is_match << "  "
                        << "  ele(" << ie_store << ") pt=" << m_susyNt.ele()->at(ie_store).pt
                        << ", eta=" << m_susyNt.ele()->at(ie_store).eta
                        << "    muo(" << im_store << ") pt=" << m_susyNt.muo()->at(im_store).pt
                        << ", eta=" << m_susyNt.muo()->at(im_store).eta << endl;
                }
            }
        } // im_store
    } // ie_store
}
