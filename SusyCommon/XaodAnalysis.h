#ifndef SusyCommon_XaodAnalysis_h
#define SusyCommon_XaodAnalysis_h


//Infrastructure
#ifdef ROOTCORE
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#endif

//ASG
#include "AsgTools/ToolHandle.h"

//xAOD
#include "xAODEventInfo/EventInfo.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/xAODTruthHelpers.h"
#include "xAODCore/ShallowCopy.h"
#include "TrigConfHLTData/HLTTriggerElement.h"
#include "TrigConfxAOD/xAODConfigTool.h"
#include "xAODTrigger/TrigNavigation.h"
#include "xAODTrigEgamma/TrigElectron.h"
#include "xAODTrigEgamma/TrigElectronContainer.h"

//CP systematics
#include "PATInterfaces/SystematicVariation.h"
#include "PATInterfaces/SystematicRegistry.h"
#include "PATInterfaces/SystematicCode.h"

//Tools
#include "ElectronEfficiencyCorrection/AsgElectronEfficiencyCorrectionTool.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "ElectronPhotonSelectorTools/AsgPhotonIsEMSelector.h"
class AsgElectronChargeIDSelectorTool;
#include "PileupReweighting/PileupReweightingTool.h"
#include "AsgTools/ToolHandle.h"
#include "MuonEfficiencyCorrections/MuonEfficiencyScaleFactors.h"
#include "MuonSelectorTools/MuonSelectionTool.h"
#include "IsolationSelection/IsolationSelectionTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TauAnalysisTools/TauSelectionTool.h"
#include "TauAnalysisTools/TauEfficiencyCorrectionsTool.h"
#include "TauAnalysisTools/TauTruthMatchingTool.h"
#include "TauAnalysisTools/TauTruthTrackMatchingTool.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

//ROOT
#include "TSelector.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTreeFormula.h"
#include "TBits.h"
class TDirectory;


//SUSY
#include "SUSYTools/SUSYObjDef_xAOD.h"

//SusyNtuple
#include "SusyNtuple/ElectronId.h"
#include "SusyNtuple/MuonId.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/SusyNt.h"
#include "SusyNtuple/SusyNtSys.h"
#include "SusyNtuple/TriggerTools.h"

//SusyCommon
#include "SusyCommon/SusyObjId.h"
#include "SusyCommon/MCType.h"

using namespace Susy;
using namespace NtSys;

// fw declarations
namespace TrigConf {
    class xAODConfigTool;
}
namespace Trig {
    class TrigDecisionTool;
    class FeatureContainer;
}


namespace Susy {

    const double MeV2GeV=1.0e-3;

    class XaodAnalysis : public TSelector
    {
        public :

            // TSelector Overrides
            virtual void Init(TTree* tree);
            virtual Bool_t Notify() { return kTRUE; }
            virtual void SlaveBegin(TTree* tree);
            virtual Bool_t Process(Long64_t entry);
            virtual void SlaveTerminate(){};
            virtual void Terminate();
            virtual Int_t version() const { return 3; } // need >=2 for ROOT

            ///////////////////////////////////////////////////////////////////
            // XaodAnalysis
            ///////////////////////////////////////////////////////////////////

            XaodAnalysis();
            virtual ~XaodAnalysis();

            virtual void set_debug(int dbg_level);
            virtual bool dbg() { return m_dbg; }

            virtual bool set_chain(TChain* chain);
            virtual TChain* chain() { return m_input_chain; }

            virtual bool set_input_container(std::string name);
            std::string input_container() { return m_input_container_name; }
            virtual void set_output_container(std::string name);
            std::string output_container() { return m_output_container_name; }

            virtual bool data_or_mc_from_name(const TString &s);

            virtual bool set_mctype(MCType type);
            MCType mc_type() { return m_mc_type; }
            virtual bool mc() { return m_isMC; }
            virtual bool data15() { return m_is_data15; }
            virtual bool data16() { return m_is_data16; }

            virtual void set_af2(bool is_af2) { m_is_af2 = is_af2; }
            virtual bool af2() { return m_is_af2; }

            virtual void set_write(bool write) { m_write_ntuple = write; }
            virtual bool fill_nt() { return m_write_ntuple; }

            virtual void run_systematics(bool run_sys) { m_sys = run_sys; }
            virtual bool sys() { return m_sys; }

            // method to collect the sumw information from CutBookKeepers
            void get_sumw(TTree* tree);
            bool collect_cutbooks(xAOD::TEvent& event, int file_idx);

            // check what time of sample this is
            bool is_derivation_from_metadata(TTree* tree);
            TDirectory* get_directory_from_chain(TTree* tree);

            DataStream stream_from_input_container(const TString &s, bool isMC);

            // tool initialization
            void initialize_local_tools();
            std::string default_grl_file();
            bool initialize_GRL_tool();
            void initialize_electron_tools();
            void initialize_photon_tools();
            void initialize_muon_tools();
            void initialize_tau_tools();
            void initialize_isolation_tools();

            // initialize SUSYTools
            void initialize_SUSYTools();

            // systematics
            void get_systematic_list();


        private :
            int m_dbg; // verbosity level
            bool m_isMC;
            bool m_is_af2;
            bool m_write_ntuple; // produce the output susyNt file
            bool m_sys; // run systematics

            TChain* m_input_chain; // input TChain of DAOD

            std::string m_input_container_name;
            std::string m_output_container_name;

            MCType m_mc_type;

            // data stream
            bool m_is_derivation;
            DataStream m_stream;
            bool m_is_data15;
            bool m_is_data16;

            // interfaces with the xAOD EDM
            xAOD::TEvent m_event;
            xAOD::TStore m_store;

            // path to ROOTCORE data/ dir
            std::string m_data_dir;

            // sumw counters
            uint64_t m_nEventsProcessed;
            double m_sumOfWeights;
            double m_sumOfWeightsSquared;


            //////////////////////////////////////////////
            // local ASG tools
            //////////////////////////////////////////////

            // GRL
            std::string m_grl_filename;
            GoodRunsListSelectionTool* m_grl_tool;

            // electron ID
            asg::AnaToolHandle<IAsgElectronLikelihoodTool> m_elecSelLikelihoodVeryLoose;
            asg::AnaToolHandle<IAsgElectronLikelihoodTool> m_elecSelLikelihoodLoose;
            asg::AnaToolHandle<IAsgElectronLikelihoodTool> m_elecSelLikelihoodLooseBLayer;
            asg::AnaToolHandle<IAsgElectronLikelihoodTool> m_elecSelLikelihoodMedium;
            asg::AnaToolHandle<IAsgElectronLikelihoodTool> m_elecSelLikelihoodTight;

            // photon ID
            asg::AnaToolHandle<IAsgPhotonIsEMSelector> m_photonSelLoose;
            asg::AnaToolHandle<IAsgPhotonIsEMSelector> m_photonSelTight;

            // muon ID
            asg::AnaToolHandle<CP::IMuonSelectionTool> m_muonSelectionToolVeryLoose;
            asg::AnaToolHandle<CP::IMuonSelectionTool> m_muonSelectionToolLoose;
            asg::AnaToolHandle<CP::IMuonSelectionTool> m_muonSelectionToolMedium;
            asg::AnaToolHandle<CP::IMuonSelectionTool> m_muonSelectionToolTight;

            // tau ID
            asg::AnaToolHandle<TauAnalysisTools::ITauSelectionTool> m_tauSelToolLoose;
            asg::AnaToolHandle<TauAnalysisTools::ITauSelectionTool> m_tauSelToolMedium;
            asg::AnaToolHandle<TauAnalysisTools::ITauSelectionTool> m_tauSelToolTight;
            TauAnalysisTools::TauTruthMatchingTool *m_tauTruthMatchingTool;

            // lepton/photon isolation
            asg::AnaToolHandle<CP::IIsolationSelectionTool> m_isoToolGradientLooseTight;
            asg::AnaToolHandle<CP::IIsolationSelectionTool> m_isoToolGradientTightCalo;
            asg::AnaToolHandle<CP::IIsolationSelectionTool> m_isoToolLooseTrackOnlyLoose;
            asg::AnaToolHandle<CP::IIsolationSelectionTool> m_isoToolLoose;
            asg::AnaToolHandle<CP::IIsolationSelectionTool> m_isoToolTight;

            //////////////////////////////////////////////
            // SUSYTools instances
            //////////////////////////////////////////////
            bool m_run_oneST;
            ST::SUSYObjDef_xAOD* m_susyObj[SusyObjId::Invalid];
            SusyObjId m_eleIDDefault;

            //////////////////////////////////////////////
            // CP systematics
            //////////////////////////////////////////////
            std::vector<CP::SystematicSet> sysList;
            std::vector<ST::SystInfo> systInfoList;

    }; // class XaodAnalysis

} // namespace Susy




//_________________[don't go below here]_________________[...or else]______//
#endif
