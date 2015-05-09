//  -*- c++ -*-
#ifndef SUSY_SODACCESSOR_H
#define SUSY_SODACCESSOR_H

#include "SUSYTools/SUSYObjDef_xAOD.h"

namespace Susy
{
/**
   \brief A thin wrapper around SUSYObjDef_xAOD providing access to its protected datamembers

   Sometimes we need to access some of the protected datamembers of
   SUSYObjDef_xAOD.  This class can be used as a drop-in replacement
   for SUSYObjDef_xAOD and implements const getter for all of the
   attributes.

   See test_SodAccessor.cxx for an example showing how to use it.

   davide.gerbaudo@gmail.com
   May 2015
 */
class SodAccessor : public ST::SUSYObjDef_xAOD {
public:
    SodAccessor(const std::string& name):
        SUSYObjDef_xAOD(name) {}

    const int dataSource() const { return m_dataSource; }
    const bool doJetArea() const { return m_doJetArea; }
    const bool doJetGSC() const { return m_doJetGSC; }
    const bool is8TeV() const { return m_is8TeV; }
    const bool isDerived() const { return m_isDerived; }
    const int jesNPset() const { return m_jesNPset; }
		
    const std::string eleTerm    () const { return m_eleTerm    ; }
    const std::string gammaTerm  () const { return m_gammaTerm  ; }
    const std::string tauTerm    () const { return m_tauTerm    ; }
    const std::string jetTerm    () const { return m_jetTerm    ; }
    const std::string muonTerm   () const { return m_muonTerm   ; }
    const std::string inputMap   () const { return m_inputMap   ; }
    const std::string inputMETCont   () const { return m_inputMETCont   ; }
    const std::string outMETTerm   () const { return m_outMETTerm   ; }
    const std::string eleId() const { return m_eleId; }
    const std::string eleIdBaseline() const { return m_eleIdBaseline; }
    const std::string phId() const { return m_phId; }
    const std::string muId() const { return m_muId; }
    const std::string tauId() const { return m_tauId; }
		
    const bool doEle() const { return m_doEle; }
    const bool doMuon() const { return m_doMuon; }
    const bool doTau() const { return m_doTau; }
    const bool doGamma() const { return m_doGamma; }
    const bool debug() const { return m_debug; }

/*

  DG 2015-05-08
  For the accessors below in principle we need to include the
  ToolHandle header as well as the ones for each tool.
  Do we really need them?

    const ToolHandle<IPseudoJetGetter> lcget() const { return m_lcget; }
    const ToolHandle<IEventShapeTool> edtool() const { return m_edtool; }
    //
    const ToolHandle<IJetCalibrationTool> jetCalibTool() const { return m_jetCalibTool; }
    const ToolHandle<IJERTool> jerTool() const { return m_jerTool; }
    const ToolHandle<IJERSmearingTool> jerSmearingTool() const { return m_jerSmearingTool; }
    const ToolHandle<ICPJetUncertaintiesTool> jetUncertaintiesTool() const { return m_jetUncertaintiesTool; }
    const ToolHandle<IJetSelector> jetCleaningTool() const { return m_jetCleaningTool; }
    //
    const ToolHandle<CP::IMuonEfficiencyScaleFactors> muonEfficiencySFTool() const { return m_muonEfficiencySFTool; }
    const ToolHandle<CP::IMuonCalibrationAndSmearingTool> muonCalibrationAndSmearingTool() const { return m_muonCalibrationAndSmearingTool; }
    const ToolHandle<CP::IMuonSelectionTool> muonSelectionTool() const { return m_muonSelectionTool; }
    //
    const ToolHandle<IAsgElectronEfficiencyCorrectionTool> elecEfficiencySFTool_reco() const { return m_elecEfficiencySFTool_reco; }
    const ToolHandle<IAsgElectronEfficiencyCorrectionTool> elecEfficiencySFTool_id() const { return m_elecEfficiencySFTool_id; }
    const ToolHandle<IAsgElectronEfficiencyCorrectionTool> elecEfficiencySFTool_trig() const { return m_elecEfficiencySFTool_trig; }
    const ToolHandle<CP::IEgammaCalibrationAndSmearingTool> egammaCalibTool() const { return m_egammaCalibTool; }
    const ToolHandle<IAsgElectronLikelihoodTool> elecSelLikelihood() const { return m_elecSelLikelihood; }
    const ToolHandle<IAsgElectronIsEMSelector>   elecSelIsEM() const { return m_elecSelIsEM; }
    const ToolHandle<IAsgElectronLikelihoodTool> elecSelLikelihoodBaseline() const { return m_elecSelLikelihoodBaseline; }
    const ToolHandle<IAsgElectronIsEMSelector>   elecSelIsEMBaseline() const { return m_elecSelIsEMBaseline; }
    const ToolHandle<IAsgPhotonIsEMSelector>     photonSelIsEM() const { return m_photonSelIsEM; }
    const ToolHandle<IAsgPhotonEfficiencyCorrectionTool> photonEfficiencySFTool() const { return m_photonEfficiencySFTool; }
    //
    const ToolHandle<TauAnalysisTools::ITauSelectionTool> tauSelTool() const { return m_tauSelTool; }
    const ToolHandle<TauAnalysisTools::ITauSmearingTool> tauSmearingTool() const { return m_tauSmearingTool; }
    const ToolHandle<TauAnalysisTools::ITauEfficiencyCorrectionsTool> tauEffTool() const { return m_tauEffTool; }
    //
    const ToolHandle<IBTaggingEfficiencyTool> btagTool() const { return m_btagTool; }
    //
    const ToolHandle<IMETRebuilder> metRebuilder() const { return m_metRebuilder; }
    const ToolHandle<IMETSystematicsTool> metSystTool() const { return m_metSystTool; }

    const EventShapeCopier escopier() const { return m_escopier; }
*/

}; // SodAccessor
} // Susy

#endif
