#ifndef SUSYCOMMON_SUSYNTMAKER_H
#define SUSYCOMMON_SUSYNTMAKER_H

//SusyNtuple
#include "SusyNtuple/SusyNtObject.h"
#include "SusyNtuple/TriggerTools.h"

//SusyCommon
#include "SusyCommon/XaodAnalysis.h"
#include "SusyCommon/SystematicMapping.h"

//ROOT
#include "TStopwatch.h"

//std/stl
#include <iostream>
#include <string>
#include <vector>

//Tools
namespace Root { class TElectronEfficiencyCorrectionTool; }


namespace Susy {
class SusyNtMaker : public XaodAnalysis
{
    public :
        SusyNtMaker();
        virtual ~SusyNtMaker();

        //////////////////////////////////////////////////////////////////////
        // TSelector
        //////////////////////////////////////////////////////////////////////
        //virtual void Init(TTree* tree);
        virtual void SlaveBegin(TTree* tree);
        virtual Bool_t Process(Long64_t entry);
        virtual void Terminate();

        //////////////////////////////////////////////////////////////////////
        // SusyNtMaker
        //////////////////////////////////////////////////////////////////////

        bool running_options_are_valid();

        void fill_event_trigger_histo();

        std::vector<int> susy_finalstate(); // get final state based on produced SUSY particles, c.f. SUSYTools/SUSYCrossSection.h

        void initialize_counters();

        // output ntuple
        void initialize_output_tree();
        void save_output_tree();
        void write_metadata();

        // event and object selection
        bool pass_event_level_selection();
        void fill_nominal_objects();
        bool pass_object_level_selection();

        // summaries
        std::string counter_summary();
        std::string timer_summary();

        ////////////////////////////////////////////
        // storing object variables
        ////////////////////////////////////////////
        void fill_nt_variables();

        // Event
        void fill_event_variables();
        std::vector<float> get_mc_weights(const xAOD::EventInfo* eventinfo);

        // Electron
        void fill_electron_variables();
        void store_electron(const xAOD::Electron& in, int ele_idx);

        // Muon
        void fill_muon_variables();
        void store_muon(const xAOD::Muon& in, const xAOD::MuonContainer& muons);

        // Jet
        void fill_jet_variables();
        void store_jet(const xAOD::Jet& in);

        // Tau
        void fill_tau_variables();
        void store_tau(const xAOD::TauJet& in);

        // MET
        void fill_met_variables(SusyNtSys sys = NtSys::NOM);

        void clear_event();


    private :

        bool m_flags_checked;

        TFile* m_file_outtree;
        TTree* m_outtree; 

        TH1D* h_passTrigLevel; // histo storing event level fired triggers

        Susy::SusyNtObject m_susyNt; // SusyNt interface

        // timer
        TStopwatch m_timer;

        // SUSY final state
        int m_susyFinalState;

        // cutflow counters
        const std::vector<std::string> cutflow_labels();
        std::vector< size_t > m_cutstageCounters;

        // object counters
        uint    n_pre_ele;
        uint    n_pre_muo;
        uint    n_pre_tau;
        uint    n_pre_jet;
        uint    n_base_ele;
        uint    n_base_muo;
        uint    n_base_tau;
        uint    n_base_jet;
        uint    n_sig_ele;
        uint    n_sig_muo;
        uint    n_sig_tau;
        uint    n_sig_jet;


}; // clas SusyNtMaker
} // namespace Susy


#endif
