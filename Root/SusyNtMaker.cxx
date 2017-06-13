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
//void SusyNtMaker::Init(TTree *tree)
//{
//}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::SlaveBegin(TTree* tree)
{
    cout << "SusyNtMaker::SlaveBegin" << endl;
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

    // clear the output storage indices
    clear_output_objects();

    // start the PRW tool since all tools depend on it downstream
    for(int susyObjId : Susy::leptonIds()) {
        m_susyObj[susyObjId]->ApplyPRWTool();
        if(m_run_oneST) break;
    }

    const xAOD::EventInfo* eventinfo = XaodAnalysis::xaodEventInfo();

    if(dbg() || chainEntry % 500 == 0) {
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

    // SUSY final state
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
        sample_event_triggers(); // check if triggers fired at event level (not matching)

        // store objects to the output susyNt
        fill_nt_variables();

        // fill the output tree
        int bytes = m_outtree->Fill();
        if(bytes < 0) {
            cout << "SusyNtMaker::Process    ERROR Unable to fill output tree, abort (fill returns " << bytes << ")" << endl;
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

    fill_object_cleaning_flags();

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
    int year = -1;
    if(mc()) {
        year = m_susyObj[m_eleIDDefault]->treatAsYear();
    }
    else {
        if(data15()) year = 2015;
        else if(data16()) year = 2016;
    }
    if(year<0) cout << "SusyNtMaker::fill_event_variables    WARNING treatAsYear was not found correctly for event " << eventinfo->eventNumber() << "!" << endl;
    evt->treatAsYear = year;

    evt->isMC = mc();
    evt->mcChannel = mc() ? eventinfo->mcChannelNumber() : 0;
    evt->w = mc() ? eventinfo->mcEventWeight() : 1;

    cout << "SusyNtMaker::fill_event_variables    mc channel : " << evt->mcChannel << "  year: " << evt->mcChannel << endl;

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
                cout << "SusyNtMaker::fill_event_variables    cannot configure SUSYTools for systematic " << sys.name() << " (" << SusyNtSysNames[ourSys] << ")" << endl;
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
        for(auto& i : m_preJets) {
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

    if(dbg()>=15) cout << "SusyNtMaker::fill_met_variables    Filling MET (sys=" << SusyNtSysNames[sys] << ")" << endl;

    xAOD::MissingETContainer::const_iterator met_it = xaodMET()->find("Final");

    if(dbg()>=15) {
        cout << "SusyNtMaker::fill_met_variables    Dumping MET container " << endl;
        for(auto it = xaodMET()->begin(), end = xaodMET()->end(); it!=end; ++it) {
            cout << "SusyNtMaker::fill_met_variables      > MET " << (*it)->name() << endl;
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

    if(dbg()>=15) cout << "SusyNtMaker::fill_met_variables    (sys=" << SusyNtSysNames[sys] << ") (Et,phi,pt)=("<<met->Et<<","<<met->phi<<","<<met->lv().Pt()<<")"<<endl;

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
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::fill_track_met_variables(SusyNtSys sys)
{

    if(dbg()>-15) cout << "SusyNtMaker::fill_track_met_variables    Filling TrackMet (sys="<<SusyNtSysNames[sys]<<")" <<endl;

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
    if(mc()) {
        //////////////////////////////////////
        // Lepton SF
        // - one for each electron LH WP
        // - (only Loose, Medium, Tight for now)
        //////////////////////////////////////
        // signature: (input electron, bool doRecoSF, bool doIDSF, bool doTrigSF, bool doIsoSF, string trigExpr)
        // default trigExpr: "e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose"
        
        if(in.pt()*MeV2GeV >= 25.) {
            std::string ele_trig;
            std::string ele_trig15 = "e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose";
            ele_trig = ele_trig15;
            std::string ele_trig16 = "e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0";
            ele_trig16 = ele_trig15;
            if(m_run_oneST) {
                if(m_susyObj[m_eleIDDefault]->treatAsYear()==2016)
                    ele_trig = ele_trig16;
                out.eleEffSF[ElectronId::TightLLH] =  m_susyObj[m_eleIDDefault]->GetSignalElecSF (in, recoSF, idSF, trigSF, isoSF);
                out.eleTrigSF[ElectronId::TightLLH] = m_susyObj[m_eleIDDefault]->GetSignalElecSF (in, false, false, true, false, ele_trig);
                //out.eleCHFSF[ElectronId::TightLLH] = m_susyObj[m_eleIDDefault]->GetSignalElecSF (in, false, false, false, false, ele_trig, true);
            }
            else {
                if(m_susyObj[SusyObjId::eleTightLLH]->treatAsYear()==2016)
                    ele_trig = ele_trig16;
        
                out.eleEffSF[ElectronId::TightLLH] =  m_susyObj[SusyObjId::eleTightLLH]->GetSignalElecSF (in, recoSF, idSF, trigSF, isoSF);
                out.eleTrigSF[ElectronId::TightLLH] = m_susyObj[SusyObjId::eleTightLLH]->GetSignalElecSF (in, false, false, true, false, ele_trig);
                //out.eleCHFSF[ElectronId::TightLLH] = m_susyObj[SusyObjId::eleTightLLH]->GetSignalElecSF (in, false, false, false, false, ele_trig, true);
                ele_trig = ele_trig15;
        
        
                if(m_susyObj[SusyObjId::eleMediumLLH]->treatAsYear()==2016)
                    ele_trig = ele_trig16;
                out.eleEffSF[ElectronId::MediumLLH] = m_susyObj[SusyObjId::eleMediumLLH]->GetSignalElecSF(in, recoSF, idSF, trigSF, isoSF);
                out.eleTrigSF[ElectronId::MediumLLH] = m_susyObj[SusyObjId::eleMediumLLH]->GetSignalElecSF(in, false, false, true, false, ele_trig);
                //out.eleCHFSF[ElectronId::MediumLLH] = m_susyObj[SusyObjId::eleMediumLLH]->GetSignalElecSF(in, false, false, false, false, ele_trig, true);
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
    if(m_isMC && m_sys && (out.veryLooseLLH || out.looseLLH || out.mediumLLH || out.tightLLH) && in.pt()*MeV2GeV > 25.) {
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

            std::string ele_trig;
            std::string ele_trig15 = "e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose";
            ele_trig = ele_trig15;
            std::string ele_trig16 = "e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0";
            ele_trig16 = ele_trig15;
            
            vector<float> sf_trig;
            sf_trig.assign(ElectronId::ElectronIdInvalid, 1);
            if(m_run_oneST) {
                if(m_susyObj[m_eleIDDefault]->treatAsYear()==2016)
                    ele_trig = ele_trig16;
                sf_trig[ElectronId::TightLLH]  = m_susyObj[m_eleIDDefault] ->GetSignalElecSF(in, false, false, true, false, ele_trig);
                sf_trig[ElectronId::MediumLLH] = m_susyObj[m_eleIDDefault]->GetSignalElecSF(in, false, false, true, false, ele_trig);
            }
            else {
                if(m_susyObj[m_eleIDDefault]->treatAsYear()==2016)
                    ele_trig = ele_trig16;
                sf_trig[ElectronId::TightLLH]  = m_susyObj[SusyObjId::eleTightLLH] ->GetSignalElecSF(in, false, false, true, false, ele_trig);
                sf_trig[ElectronId::MediumLLH] = m_susyObj[SusyObjId::eleMediumLLH]->GetSignalElecSF(in, false, false, true, false, ele_trig);
            }

            // there are no isolation SF's for electrion ID looseLH
            //sf[ElectronId::LooseLH]  = m_susyObj[SusyObjId::eleLooseLH] ->GetSignalElecSF(in, recoSF, idSF, trigSF);
        
            for(int i=ElectronId::TightLLH; i<ElectronIdInvalid; i++){
                if     (ourSys == NtSys::EL_EFF_ID_TOTAL_Uncorr_UP)      out.errEffSF_id_up[i]   = sf[i] - out.eleEffSF[i];
                else if(ourSys == NtSys::EL_EFF_ID_TOTAL_Uncorr_DN)      out.errEffSF_id_dn[i]   = sf[i] - out.eleEffSF[i];
                else if(ourSys == NtSys::EL_EFF_Reco_TOTAL_Uncorr_UP)    out.errEffSF_reco_up[i] = sf[i] - out.eleEffSF[i];
                else if(ourSys == NtSys::EL_EFF_Reco_TOTAL_Uncorr_DN)    out.errEffSF_reco_dn[i] = sf[i] - out.eleEffSF[i];
                else if(ourSys == NtSys::EL_EFF_Iso_TOTAL_Uncorr_UP)     out.errEffSF_iso_up[i]  = sf[i] - out.eleEffSF[i];
                else if(ourSys == NtSys::EL_EFF_Iso_TOTAL_Uncorr_DN)     out.errEffSF_iso_dn[i]  = sf[i] - out.eleEffSF[i];
                else if(ourSys == NtSys::EL_EFF_Trigger_TOTAL_Uncorr_UP) out.errEffSF_trig_up[i] = sf_trig[i] - out.eleTrigSF[i];
                else if(ourSys == NtSys::EL_EFF_Trigger_TOTAL_Uncorr_DN) out.errEffSF_trig_dn[i] = sf_trig[i] - out.eleTrigSF[i];
        
                if(i==0 || i==1 || i==2){
                if(i==0)
                    cout << "ElectronId : TightLH " <<  endl;
                else if(i==1)
                    cout << "ElectronId : MediumLH " <<  endl;
                else if(i==2)
                    cout << "ElectronId : LooseLH " <<  endl;
        
                cout << "   effId           : " << out.errEffSF_id_up[i] << "  " << out.errEffSF_id_dn[i] << endl;
                cout << "   effReco         : " << out.errEffSF_reco_up[i] << "  " << out.errEffSF_reco_dn[i] << endl;
                cout << "   effIso          : " << out.errEffSF_iso_up[i] << "  " << out.errEffSF_iso_dn[i] << endl;
                cout << "   effTrig         : " << out.errEffSF_trig_up[i] << "  " << out.errEffSF_trig_dn[i] << endl;
                }
        
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
    // dantrim 2017 June 12 - TODO change how do do lepton trigger matching
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
    // dantrim 2017 June 12 - TODO change how do do lepton trigger matching
    out.trigBits = matchMuonTriggers(in);
//    cout << "testing muon trigBits" << endl;
//    int nbins = h_passTrigLevel->GetXaxis()->GetNbins();
//    for(int iTrig=0; iTrig<46; iTrig++){
//        bool bit = out.trigBits.TestBitNumber(iTrig);
//        string trigger = h_passTrigLevel->GetXaxis()->GetBinLabel(iTrig+1);
//        cout << "\t passed trigger " << trigger << "? " << (bit ? "yes" : "no") << endl;
//    }
//    cout << endl;
    out.diMuTrigMap = getDiMuTrigMap(in, muons);

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

    #warning check muon trigger SF variations once we get these SFs
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
            vector<float> sf_trig;
            sf_trig.assign(MuonId::MuonIdInvalid, 1);

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
            if(m_run_oneST) {
                sf_trig[MuonId::Medium] = m_susyObj[m_eleIDDefault]->GetTotalMuonTriggerSF(*sf_muon, static_cast<string>(trig_exp_med.Data()));
            }
            else {
                sf_trig[MuonId::Medium] = m_susyObj[SusyObjId::muoMedium]->GetTotalMuonTriggerSF(*sf_muon, static_cast<string>(trig_exp_med.Data()));
            }
            delete sf_muon;
            delete sf_muon_aux;
            //delete sfMu;

            for(int i=MuonId::VeryLoose; i < MuonId::MuonIdInvalid; i++) {
                if      (ourSys == NtSys::MUON_EFF_TRIG_STAT_UP)    out.errTrigSF_stat_up[i] = sf_trig[i] - out.muoTrigSF[i];
                else if (ourSys == NtSys::MUON_EFF_TRIG_STAT_DN)    out.errTrigSF_stat_dn[i] = sf_trig[i] - out.muoTrigSF[i];
                else if (ourSys == NtSys::MUON_EFF_TRIG_SYST_UP)    out.errTrigSF_syst_up[i] = sf_trig[i] - out.muoTrigSF[i];
                else if (ourSys == NtSys::MUON_EFF_TRIG_SYST_DN)    out.errTrigSF_syst_dn[i] = sf_trig[i] - out.muoTrigSF[i];

               // if(i == 1 || i == 2) {
               //     if(i==1)
               //         cout << "MuonId: Loose" << endl;
               //     else if(i==2)
               //         cout << "MuonId: Medium" << endl;
               //     cout << "       eff trig stat        : " << out.errTrigSF_stat_up[i] << "  " << out.errTrigSF_stat_dn[i] << endl;
               //     cout << "       eff trig syst        : " << out.errTrigSF_syst_up[i] << "  " << out.errTrigSF_syst_dn[i] << endl;
               // }
            } // i

        } // sysInfo
        for(int i : Susy::muonIds()){
            int index_to_check = (m_run_oneST==true ? (int)m_eleIDDefault : i);
            if(m_susyObj[index_to_check]->resetSystematics() != CP::SystematicCode::Ok){
                cout << "SusyNtMaker::store_muon    ERROR cannot reset SUSYTools systematics. Aborting." << endl;
                abort();
            }
            if(m_run_oneST) break;
        }
    } // ifMC && sys
    else {
        for(int i=MuonId::VeryLoose; i<MuonId::MuonIdInvalid; i++){
            out.errTrigSF_stat_up[i] = out.errTrigSF_stat_dn[i] = 0;
            out.errTrigSF_syst_up[i] = out.errTrigSF_syst_dn[i] = 0;
        }
    }

    //______________ ALL DONE WITH THE MUON ______________ //
    out.idx = (m_susyNt.muo()->size());
    m_susyNt.muo()->push_back(out);

}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::store_jet(const xAOD::Jet& in)
{
    if(dbg()>=15) cout << "SusyNtMaker::fill_jet    Storing Jet (pt=" << in.pt()*MeV2GeV << ")" << endl;

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
            if(m_sys) {
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
