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
#include "ElectronPhotonSelectorTools/AsgElectronChargeIDSelectorTool.h"
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

//std/stl
#include <tuple>

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
            XaodAnalysis();
            virtual ~XaodAnalysis();

            // TSelector Overrides
            virtual void Init(TTree* tree);
            virtual Bool_t Notify() { return kTRUE; }
            virtual void SlaveBegin(TTree* tree);
            virtual Bool_t Process(Long64_t entry);
            virtual void SlaveTerminate(){};
            virtual void Terminate();
            // TSelector Version defines TSelector flow: https://root.cern.ch/developing-tselector#Version
            virtual Int_t Version() const { return 2; } // need >=2 for ROOT

            ///////////////////////////////////////////////////////////////////
            // XaodAnalysis
            ///////////////////////////////////////////////////////////////////


            virtual void set_debug(int dbg_level);
            virtual int dbg() { return m_dbg; }

            virtual bool set_chain(TChain* chain);
            virtual TChain* chain() { return m_input_chain; }

            virtual bool set_input_container(std::string name);
            std::string input_container() { return m_input_container_name; }
            virtual void set_output_container(std::string name);
            std::string output_container() { return m_output_container_name; }

            virtual void set_production_command(std::string cmd) { m_production_command = cmd; }
            virtual std::string production_command() { return m_production_command; }

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

            virtual void set_production_tag(std::string tag) { m_nt_tag = tag; }
            virtual std::string production_tag() { return m_nt_tag; }

            virtual void set_output_name(std::string outputname) { m_output_filename = outputname; }
            virtual std::string output_name() { return m_output_filename; }

            virtual void run_systematics(bool run_sys) { m_sys = run_sys; }
            virtual bool sys() { return m_sys; }

            virtual void set_nlep_filter(int nlep);
            virtual int nlep_for_filter() { return m_nlep_filter; }

            virtual void set_trig_filter(bool doit);
            virtual bool do_trig_filter() { return m_filter_trig; }

            virtual bool do_event_filter() { return m_filter; }

            virtual void set_cont_taus(bool save_them) { m_saveContTaus = save_them; }
            virtual bool cont_taus() { return m_saveContTaus; }

            virtual void set_store_truth(bool store_em) { m_store_truth = store_em; }
            virtual bool store_truth() { return m_store_truth; }

            // method to collect the sumw information from CutBookKeepers
            void get_sumw(TTree* tree);
            bool collect_cutbooks(xAOD::TEvent& event, int file_idx);

            // check what time of sample this is
            bool is_derivation_from_metadata(TTree* tree);
            TDirectory* get_directory_from_chain(TTree* tree);
            TString get_derivation_type(xAOD::TEvent& event);

            DataStream stream_from_input_container(const TString &s, bool isMC);

            // tool initialization
            void initialize_local_tools();
            std::string default_grl_file();
            bool initialize_GRL_tool();
            void initialize_electron_tools();
            void initialize_chargeflip_tagger();
            void initialize_photon_tools();
            void initialize_muon_tools();
            void initialize_tau_tools();
            void initialize_isolation_tools();

            // initialize SUSYTools
            void initialize_SUSYTools();

            // systematics
            void get_systematic_list();

            // event cleaning flags
            void fill_event_cleaning_flags();
            bool passGRL(const xAOD::EventInfo* ei);
            bool passTTCVeto(const xAOD::EventInfo* ei);
            bool passLarErr(const xAOD::EventInfo* ei);
            bool passTileErr(const xAOD::EventInfo* ei);
            bool passSCTErr(const xAOD::EventInfo* ei);
            bool passGoodVtx();

            void fill_object_cleaning_flags();
            bool passBadJet(ST::SystInfo sysInfo, SusyNtSys sys);
            bool passBadMuon(ST::SystInfo sysInfo, SusyNtSys sys);
            bool passCosmic(ST::SystInfo sysInfo, SusyNtSys sys);


            ///////////////////////////////////////////////////////////////////
            // Accessing the xAOD containers/objects
            ///////////////////////////////////////////////////////////////////

            // clear containers
            void clear_containers(bool delete_nominal = true);
            void clear_output_objects(bool delete_nominal = true);
            void delete_shallow_copies(bool delete_nominal = true);

            virtual std::vector<std::string> xaodTriggers();


            // fill our xAOD object containers
            virtual void retrieve_xaod_collections();

            virtual const xAOD::EventInfo* xaodEventInfo();
            
            virtual const xAOD::VertexContainer* xaodVertices();

            virtual xAOD::ElectronContainer* xaodElectrons(ST::SystInfo sysInfo, SusyNtSys sys = NtSys::NOM);

            virtual xAOD::MuonContainer* xaodMuons(ST::SystInfo sysInfo, SusyNtSys sys = NtSys::NOM);

            virtual xAOD::JetContainer* xaodJets(ST::SystInfo sysInfo, SusyNtSys sys = NtSys::NOM);

            virtual xAOD::TauJetContainer* xaodTaus(ST::SystInfo sysInfo, SusyNtSys sys = NtSys::NOM);

            virtual xAOD::PhotonContainer* xaodPhotons(ST::SystInfo sysInfo, SusyNtSys sys = NtSys::NOM);

            virtual xAOD::MissingETContainer* xaodMET() { return m_metContainer; }

            virtual xAOD::MissingETContainer* xaodTrackMET() { return m_trackMetContainer; }

            virtual const xAOD::TruthParticleContainer* xaodTruthParticles();

            ///////////////////////////////////////////////////////////////////
            // Filling xAOD objects
            ///////////////////////////////////////////////////////////////////
            void fill_objects(SusyNtSys sys, ST::SystInfo sysInfo);
            void fill_baseline_objects(SusyNtSys sys, ST::SystInfo sysInfo);
            void fill_signal_objects(SusyNtSys sys, ST::SystInfo sysInfo);

            ///////////////////////////////////////////////////////////////////
            // Output to SusyNtObject
            ///////////////////////////////////////////////////////////////////
            void sample_event_triggers();

            bool dilepton_trigger_matched(const xAOD::IParticle* part1, const xAOD::IParticle* part2, std::string chain = "");

            // electrons
            bool eleIsOfType(const xAOD::Electron& in, ElectronId id);
            TBits matchElectronTriggers(const xAOD::Electron& in);

            // muons
            TBits matchMuonTriggers(const xAOD::Muon& in);
            std::map<std::string, std::vector<unsigned int>> getDiMuTrigMap(const xAOD::Muon &in, const xAOD::MuonContainer &muons);
            bool muIsOfType(const xAOD::Muon &in, MuonId id);

        protected :
            int m_dbg; // verbosity level
            std::string m_production_command;
            bool m_isMC;
            bool m_is_af2;
            bool m_write_ntuple; // produce the output susyNt file
            std::string m_nt_tag;
            std::string m_output_filename; // name of output susyNt file
            bool m_sys; // run systematics

            TString m_derivation;

            int m_nlep_filter;
            bool m_filter_trig;
            bool m_filter;
            bool m_saveContTaus;
            bool m_store_truth;

            TChain* m_input_chain; // input TChain of DAOD

            std::string m_input_container_name;
            std::string m_output_container_name;

            MCType m_mc_type;

            // data stream
            bool m_is_derivation;
            DataStream m_stream;
            bool m_is_data15;
            bool m_is_data16;

            std::vector<std::string> m_triggerNames;


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

            // charge ID
            asg::AnaToolHandle<IAsgElectronLikelihoodTool> m_electronChargeIDTool;

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

            //////////////////////////////////////////////
            // xAOD containers/objects
            //////////////////////////////////////////////

            // interfaces with the xAOD EDM
            xAOD::TEvent m_event;
            xAOD::TStore m_store;

            static const xAOD::EventInfo* retrieveEventInfo(xAOD::TEvent &e, bool dbg);
            const xAOD::EventInfo* m_xaodEventInfo;

            // primary vertices
            static const xAOD::VertexContainer* retrieveVertices(xAOD::TEvent &e, bool dbg);
            const xAOD::VertexContainer* m_xaodVertices;

            // electrons
            xAOD::ElectronContainer* m_xaodElectrons;
            xAOD::ShallowAuxContainer* m_xaodElectronsAux;

            // muons
            xAOD::MuonContainer* m_xaodMuons;
            xAOD::ShallowAuxContainer* m_xaodMuonsAux;

            // jets
            xAOD::JetContainer* m_xaodJets;
            xAOD::ShallowAuxContainer* m_xaodJetsAux;

            // taus
            xAOD::TauJetContainer* m_xaodTaus;
            xAOD::ShallowAuxContainer* m_xaodTausAux;

            // photons
            xAOD::PhotonContainer* m_xaodPhotons;
            xAOD::ShallowAuxContainer* m_xaodPhotonsAux;

            // MET
            virtual void retrieveXaodMet(ST::SystInfo sysInfo, SusyNtSys sys = NtSys::NOM);
            xAOD::MissingETContainer* m_metContainer;
            xAOD::MissingETAuxContainer* m_metAuxContainer;

            // track MET
            virtual void retrieveXaodTrackMet(ST::SystInfo sysInfo, SusyNtSys sys = NtSys::NOM);
            xAOD::MissingETContainer* m_trackMetContainer;
            xAOD::MissingETAuxContainer* m_trackMetAuxContainer;

            // truth particles
            const xAOD::TruthParticleContainer* retrieveTruthParticles(xAOD::TEvent& e, bool dbg);
            const xAOD::TruthEventContainer* m_xaodTruthEvent;
            const xAOD::TruthParticleContainer* m_xaodTruthParticles;
            xAOD::TruthParticleAuxContainer* m_xaodTruthParticlesAux;

            const xAOD::TruthParticleContainer* m_xaodTruthTauParticles;
            xAOD::TruthParticleAuxContainer* m_xaodTruthTauParticlesAux;

            // containers at nominal scale (needed for systematic variation calculations)

            // electrons nom
            xAOD::ElectronContainer* m_xaodElectrons_nom;
            xAOD::ShallowAuxContainer* m_xaodElectronsAux_nom;

            // muons nom
            xAOD::MuonContainer* m_xaodMuons_nom;
            xAOD::ShallowAuxContainer* m_xaodMuonsAux_nom;

            // jets nom
            xAOD::JetContainer* m_xaodJets_nom;
            xAOD::ShallowAuxContainer* m_xaodJetsAux_nom;

            // taus nom
            xAOD::TauJetContainer* m_xaodTaus_nom;
            xAOD::ShallowAuxContainer* m_xaodTausAux_nom;

            // photons nom
            xAOD::PhotonContainer* m_xaodPhotons_nom;
            xAOD::ShallowAuxContainer* m_xaodPhotonsAux_nom;
        
            //////////////////////////////////////////////
            // output "objects"
            //////////////////////////////////////////////

            uint32_t m_cutFlags; // event cleaning cuts

            std::vector<int> m_preElectrons;
            std::vector<int> m_preMuons;
            std::vector<int> m_preLeptons;
            std::vector<int> m_preJets;
            std::vector<int> m_preTaus;
            std::vector<int> m_contTaus;
            std::vector<int> m_prePhotons;

            std::vector<int> m_baseElectrons;
            std::vector<int> m_baseMuons;
            std::vector<int> m_baseLeptons;
            std::vector<int> m_baseJets;
            std::vector<int> m_baseTaus;
            std::vector<int> m_basePhotons;

            std::vector<int> m_sigElectrons;
            std::vector<int> m_sigMuons;
            std::vector<int> m_sigLeptons;
            std::vector<int> m_sigJets;
            std::vector<int> m_sigTaus;
            std::vector<int> m_sigPhotons;

            // nominal
            std::vector<int> m_preElectrons_nom;
            std::vector<int> m_preMuons_nom;
            std::vector<int> m_preLeptons_nom;
            std::vector<int> m_preJets_nom;
            std::vector<int> m_preTaus_nom;
            std::vector<int> m_contTaus_nom;
            std::vector<int> m_prePhotons_nom;
            
            //////////////////////////////////////////////
            // trigger bits
            //////////////////////////////////////////////
            TBits m_evtTrigBits;
            static const size_t m_nTriggerBits=64;

    }; // class XaodAnalysis

} // namespace Susy




//_________________[don't go below here]_________________[...or else]______//
#endif
