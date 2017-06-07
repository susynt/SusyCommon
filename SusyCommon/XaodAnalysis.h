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
            virtual Int_t version() const { return 2; } // need >=2 for ROOT

            ///////////////////////////////////////////////////////////////////
            // XaodAnalysis
            ///////////////////////////////////////////////////////////////////

            XaodAnalysis();
            virtual ~XaodAnalysis(){};

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


        private :
            int m_dbg; // verbosity level
            bool m_isMC;

            TChain* m_input_chain; // input TChain of DAOD

            std::string m_input_container_name;
            std::string m_output_container_name;

            MCType m_mc_type;



    }; // class XaodAnalysis

} // namespace Susy




//_________________[don't go below here]_________________[...or else]______//
#endif
