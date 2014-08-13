#include "SusyCommon/get_object_functions.h"

#include "D3PDReader/ElectronD3PDObject.h"
#include "D3PDReader/EventInfoD3PDObject.h"
#include "D3PDReader/EventShapeD3PDObject.h"
#include "D3PDReader/MuonD3PDObject.h"
#include "D3PDReader/PhotonD3PDObject.h"
#include "D3PDReader/JetD3PDObject.h"
#include "D3PDReader/MissingETCompositionD3PDObject.h"
#include "D3PDReader/PrimaryVertexD3PDObject.h"
#include "D3PDReader/TauD3PDObject.h"

#include "egammaAnalysisUtils/egammaTriggerMatching.h"
#include "egammaAnalysisUtils/EnergyRescaler.h"
#include "ElectronEfficiencyCorrection/TElectronEfficiencyCorrectionTool.h"

#include "SUSYTools/JetID.hpp"

#include "TLorentzVector.h"

#include <iostream>

//#include "CalibrationDataInterface/CalibrationDataInterfaceROOT.h"
//#include "SUSYTools/FakeMetEstimator.h"
//#include "SUSYTools/BTagCalib.h"



#include <iostream>

using namespace std;




using namespace std;

//----------------------------------------------------------
vector<int> get_taus_baseline(D3PDReader::TauD3PDObject* taus, SUSYObjDef& susyObj,
                              float ptCut, float etaCut,
                              SUSYTau::IDLevel jetBDTLevel, SUSYTau::IDLevel eleBDTLevel,
                              SUSYTau::IDLevel muonLevel, SystErr::Syste sys, bool isD3PD1512)
{
    vector<int> taus_base;
    TLorentzVector tmpLV; // Used to calculate E in p1328 d3pds
    for(int iTau=0; iTau < taus->n(); iTau++){
        const D3PDReader::TauD3PDObjectElement* tau = & (*taus)[iTau];
        tmpLV.SetPtEtaPhiM(tau->pt(), tau->eta(), tau->phi(), tau->m());
        if(susyObj.FillTau(iTau, tau->pt(), tau->eta(), tau->leadTrack_eta(), tau->phi(), tau->Et(),
                           tau->charge(), tau->numTrack(),
                           tau->JetBDTSigLoose(), tau->JetBDTSigMedium(), tau->JetBDTSigTight(),
                           tau->EleBDTLoose(), tau->EleBDTMedium(), tau->EleBDTTight(),
                           tau->BDTEleScore(), tau->muonVeto(),
                           jetBDTLevel, eleBDTLevel, muonLevel,
                           ptCut, etaCut, sys, isD3PD1512))
        {
            taus_base.push_back(iTau);
        }
    }

    return taus_base;
}
//----------------------------------------------------------
vector<int> get_taus_signal(D3PDReader::TauD3PDObject* taus, vector<int>& taus_base,
                            SUSYObjDef& susyObj, float ptCut, float etaCut,
                            SUSYTau::IDLevel jetBDTLevel, SUSYTau::IDLevel eleBDTLevel,
                            SUSYTau::IDLevel muonLevel, SystErr::Syste sys, bool isD3PD1512)
{
    vector<int> taus_signal;
    TLorentzVector tmpLV; // Used to calculate E in p1328 d3pds
    for(unsigned int i=0; i < taus_base.size(); i++){
        int iTau = taus_base[i];
        const D3PDReader::TauD3PDObjectElement* tau = & (*taus)[iTau];
        tmpLV.SetPtEtaPhiM(tau->pt(), tau->eta(), tau->phi(), tau->m());

        // Easiest way to do this is to call SUSYObjDef::FillTau again.
        // The TLV will be calculated again, but oh well.
        // TODO: change to SUSYObjDef.IsSignalTau
        if(susyObj.FillTau(iTau, tau->pt(), tau->eta(), tau->leadTrack_eta(), tau->phi(), tau->Et(),
                           tau->charge(), tau->numTrack(),
                           tau->JetBDTSigLoose(), tau->JetBDTSigMedium(), tau->JetBDTSigTight(),
                           tau->EleBDTLoose(), tau->EleBDTMedium(), tau->EleBDTTight(),
                           tau->BDTEleScore(), tau->muonVeto(),
                           jetBDTLevel, eleBDTLevel, muonLevel,
                           ptCut, etaCut, sys))
        {
            taus_signal.push_back(iTau);
        }
    }

    return taus_signal;
}
//----------------------------------------------------------
vector<int> get_electrons_met(D3PDReader::MissingETCompositionD3PDObject *elMetEgamma10NoTau, SUSYObjDef &susyobj)
{
    vector<int> electrons_met;
    for (int iel=0; iel<elMetEgamma10NoTau->n(); iel++) {
        //if(electrons->MET_Egamma10NoTau_wet()->at(iel).at(0) != 0){
        // DG 2014-06-01: ntupcommon migration, I think that now we don't need to check this anymore...to be tested
        if(true){
            electrons_met.push_back(iel);
        }
    }
    return electrons_met;
}
//----------------------------------------------------------
vector<int> get_electrons_baseline(D3PDReader::ElectronD3PDObject *electrons,
                                   D3PDReader::MissingETCompositionD3PDObject *metEgamma10NoTau,
                                   bool kIsData, int run_number,
                                   SUSYObjDef &susyobj, float etcut, float etacut, SystErr::Syste el_syst)
{
    vector<int> electrons_base;
    for (int iel=0; iel<electrons->n(); iel++) {
        const D3PDReader::ElectronD3PDObjectElement* electron = & (*electrons)[iel];
        // electron energy resolution syst is done in susytools
        float wet = metEgamma10NoTau->wet()->at(iel).at(0);
        if (susyobj.FillElectron(iel,
                                 electron->eta(), electron->phi(),
                                 electron->cl_eta(), electron->cl_phi(), electron->cl_E(),
                                 electron->tracketa(), electron->trackphi(),
                                 electron->author(), electron->mediumPP(),
                                 electron->OQ(), electron->nPixHits(), electron->nSCTHits(),
                                 wet, etcut, etacut, el_syst))
        {
            electrons_base.push_back(iel);
        }
    }
    return electrons_base;
}
//----------------------------------------------------------
// Get isolation corrections
float elPtConeCorr(int elIdx, D3PDReader::ElectronD3PDObject *electrons,  vector<int> electrons_base,
                   D3PDReader::MuonD3PDObject *muons, vector<int> muons_base, SUSYObjDef &susyobj,
                   bool removeLeps)
{

    TLorentzVector lv = susyobj.GetElecTLV(elIdx);
    double ptcone = electrons->ptcone30()->at(elIdx);

    if (removeLeps){
	//Loop over base electrons
	for(uint iEl=0; iEl<electrons_base.size(); iEl++){
	    if(elIdx == electrons_base[iEl]) continue;
	    TLorentzVector lv1 = susyobj.GetElecTLV(electrons_base[iEl]);
	    float dR = lv.DeltaR(lv1);
	    if(dR < 0.3) ptcone -= electrons->trackpt()->at(electrons_base[iEl]);
	}
	// Loop over base muons
	for(uint iMu=0; iMu<muons_base.size(); iMu++){
	    TLorentzVector lv1 = susyobj.GetMuonTLV(muons_base[iMu]);
	    float dR = lv.DeltaR(lv1);
	    if(dR < 0.3){
		double idTrackPt = muons->id_qoverp_exPV()->at(muons_base[iMu])!=0.?
		    fabs(sin(muons->id_theta_exPV()->at(muons_base[iMu]))/muons->id_qoverp_exPV()->at(muons_base[iMu])) : 0.;
		ptcone -= idTrackPt;
	    }
	}
    }
    return ptcone;
}
//----------------------------------------------------------
float elEtTopoConeCorr(int elIdx, D3PDReader::ElectronD3PDObject *electrons,
                       vector<int> electrons_base,  SUSYObjDef &susyobj, int nGoodVtx,
                       bool isData, bool removeLeps) {
    float mev2gev=1000.0;
    TLorentzVector lv = susyobj.GetElecTLV(elIdx);
    float etconeSlope = isData? 0.02015*mev2gev : 0.01794*mev2gev;
    float etcone = (electrons->topoEtcone30_corrected()->at(elIdx) - etconeSlope*nGoodVtx);
    if (removeLeps){
        for(uint iEl=0; iEl<electrons_base.size(); iEl++){
            if(elIdx == electrons_base[iEl]) continue;
            TLorentzVector lv1 = susyobj.GetElecTLV(electrons_base[iEl]);
            float dR = lv.DeltaR(lv1);
            if( dR < 0.28  ) {
                etcone -= electrons->cl_E()->at(electrons_base[iEl]) / cosh( electrons->cl_eta()->at(electrons_base[iEl]));
            }
        }
    }
    return etcone;
}
//----------------------------------------------------------
vector<int> get_electrons_signal(D3PDReader::ElectronD3PDObject *electrons, vector<int> electrons_base,
                                 D3PDReader::MuonD3PDObject *muons, vector<int> muons_base,
                                 unsigned int nGoodVtx, bool isData, SUSYObjDef &susyobj,
                                 float etcut, float ptconeRelCut, float etconeRelCut, float d0Sig_cutval,
                                 float z0SinTheta_cutval, bool removeLeps)
{
    vector<int> electrons_signal;
    for (unsigned int iel=0; iel<electrons_base.size(); iel++) {
        int eleIdx = electrons_base[iel];

        // Proxy classes lead to cleaner code
        const D3PDReader::ElectronD3PDObjectElement* electron = & (*electrons)[eleIdx];
        const TLorentzVector* lv = & susyobj.GetElecTLV(eleIdx);

        // Impact parameter - we might move to trackIPEstimate_*_unbiasedpvunbiased variables in the future
        double d0         = electron->trackd0pv();
        double errD0      = electron->tracksigd0pv();
        double z0         = electron->trackz0pv();
        double d0Sig      = fabs(d0/errD0);
        double z0SinTheta = fabs(z0*sin(lv->Theta()));

        int isTightPP = 0;

        // Use D3PD isEM
        isTightPP = electron->tightPP();

        // New isolation criteria, based on latest version of
        // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/SUSYDirectGauginos#Lepton_isolation_and_impact_para

	float ptcone30Corr = elPtConeCorr(eleIdx, electrons, electrons_base, muons, muons_base , susyobj, removeLeps);
	float ptcone30Rel = ptcone30Corr / lv->Pt();

	float etcone30Corr = elEtTopoConeCorr(eleIdx, electrons, electrons_base, susyobj, nGoodVtx, isData, removeLeps);
        float etcone30RelFullCorr = etcone30Corr / lv->Pt();

        // Need to apply the isolation cuts separately from SUSYTools,
        // but we can trick it into applying the correct d0 and z0 cuts
        //if(susyobj.IsSignalElectron(eleIdx, isTightPP, electron->ptcone20(), eled0sig, elez0sig, etcut, 0.1, d0_cutval, z0_cutval))
        if(susyobj.IsSignalElectron(eleIdx, isTightPP, 0, d0Sig, z0SinTheta, etcut, 0, d0Sig_cutval, z0SinTheta_cutval)) {

            //if(ptcone30Rel < 0.16 && etcone30RelFullCorr < 0.18)
            if(ptcone30Rel < ptconeRelCut && etcone30RelFullCorr < etconeRelCut){
                electrons_signal.push_back(eleIdx);
            }

        }
    }
    return electrons_signal;
}
//----------------------------------------------------------
vector<int> get_muons_baseline(D3PDReader::MuonD3PDObject *muons, bool kIsData, SUSYObjDef &susyobj, float ptcut,
                               float etacut, SystErr::Syste mu_syst)
{
    //cout << "get_muons_baseline" << endl;
    vector<int> muon_base;
    for(int imuon=0; imuon<muons->n(); imuon++) {
        // Proxy classes lead to cleaner code
        const D3PDReader::MuonD3PDObjectElement* muon = & (*muons)[imuon];
        if( susyobj.FillMuon(imuon,
                             muon->pt(), muon->eta(), muon->phi(),
                             muon->me_qoverp_exPV(), muon->id_qoverp_exPV(),
                             muon->me_theta_exPV(), muon->id_theta_exPV(), muon->id_theta(),
                             muon->charge(),
                             muon->isCombinedMuon(), muon->isSegmentTaggedMuon(),
                             muon->loose(), muon->nPixHits(),
                             muon->nPixelDeadSensors(), muon->nPixHoles(), muon->nSCTHits(),
                             muon->nSCTDeadSensors(), muon->nSCTHoles(), muon->nTRTHits(),
                             muon->nTRTOutliers(),
                             ptcut, etacut, mu_syst) )
            muon_base.push_back(imuon);
    }

    return muon_base;
}
//----------------------------------------------------------
// ptcone isolation correction - subtract overlapping leptons
float muPtConeCorr(int muIdx, D3PDReader::MuonD3PDObject *muons, vector<int> muons_base,
                   D3PDReader::ElectronD3PDObject *electrons,  vector<int> electrons_base, SUSYObjDef &susyobj,
                   int nGoodVtx, bool isData, bool removeLeps)
{
    float mev2gev=1000.0;
    TLorentzVector lv = susyobj.GetMuonTLV(muIdx);
    double ptconeSlope = isData? 0.01098*mev2gev : 0.00627*mev2gev;
    double ptcone = muons->ptcone30()->at(muIdx) - ptconeSlope*nGoodVtx;
    if (removeLeps){
        for(uint iEl=0; iEl<electrons_base.size(); iEl++){
            TLorentzVector lv1 = susyobj.GetElecTLV(electrons_base[iEl]);
            float dR = lv.DeltaR(lv1);
            if(dR < 0.3) ptcone -= electrons->trackpt()->at(electrons_base[iEl]);
        }
        for(uint iMu=0; iMu<muons_base.size(); iMu++){
            if(muIdx == muons_base[iMu]) continue;
            TLorentzVector lv1 = susyobj.GetMuonTLV(muons_base[iMu]);
            float dR = lv.DeltaR(lv1);
            if(dR < 0.3){
                double idTrackPt = muons->id_qoverp_exPV()->at(muons_base[iMu])!=0.?
                    fabs(sin(muons->id_theta_exPV()->at(muons_base[iMu]))/ muons->id_qoverp_exPV()->at(muons_base[iMu])) : 0.;
                ptcone -= idTrackPt;
            }
        }
    }
    return ptcone;
}
//----------------------------------------------------------
vector<int> get_muons_signal(D3PDReader::MuonD3PDObject *muons, vector<int> muons_base,
			     D3PDReader::ElectronD3PDObject *electrons,  vector<int> electrons_base,
			     unsigned int nGoodVtx, bool isData, SUSYObjDef &susyobj,
                             float ptcut, float ptconeRelCut, float d0Sig_cutval, float z0SinTheta_cutval,
			     bool removeLeps)
{
    vector<int> muon_signal;
    for(unsigned int imuon=0; imuon<muons_base.size(); imuon++) {
        int muIdx = muons_base[imuon];

        // Proxy classes lead to cleaner code
        const D3PDReader::MuonD3PDObjectElement* muon = & (*muons)[muIdx];
        const TLorentzVector* lv = & susyobj.GetMuonTLV(muIdx);

        // Impact parameter - we might move to trackIPEstimate_*_unbiasedpvunbiased variables in the future
        double d0         = muon->d0_exPV();
        double errD0      = sqrt(muon->cov_d0_exPV());
        double z0         = muon->z0_exPV();
        double d0Sig      = fabs(d0/errD0);
        double z0SinTheta = fabs(z0*sin(lv->Theta()));

        // New isolation criteria, based on latest version of
        // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/SUSYDirectGauginos#Lepton_isolation_and_impact_para
	float ptcone30Corr = muPtConeCorr(muIdx, muons, muons_base, electrons, electrons_base, susyobj, nGoodVtx, isData, removeLeps);
	double ptcone30RelCorr = ptcone30Corr / lv->Pt();


        // At the moment, there is only one muon iso cut, so we can get SUSYTools to apply it for us.
        // However, we need to apply the d0 and z0 cuts separately
        if(susyobj.IsSignalMuon(muIdx, ptcone30RelCorr, ptcut, ptconeRelCut)) {
            //if(kIs2L || muod0sig < d0_cutval)
            if(d0Sig < d0Sig_cutval && z0SinTheta < z0SinTheta_cutval){
                muon_signal.push_back(muIdx);
            }
        }
    }

    return muon_signal;
}
//----------------------------------------------------------
vector<int> get_jet_baseline(D3PDReader::JetD3PDObject *jets, D3PDReader::PrimaryVertexD3PDObject *vertex,
                             D3PDReader::EventInfoD3PDObject *eventInfo, D3PDReader::EventShapeD3PDObject *evtShape,
                             bool kIsData, SUSYObjDef &susyobj,
                             float etcut, float etacut, SystErr::Syste whichsyste, bool nosmear, vector<int> &goodjets)
{
    vector<int> jet_index;
    int local_truth_label = 0;
    for(int iJet=0; iJet<jets->n(); iJet++) {
        const D3PDReader::JetD3PDObjectElement* jet = & (*jets)[iJet];
        susyobj.FillJet(iJet, jet->pt(), jet->eta(), jet->phi(), jet->E(),
                        jet->constscale_eta(), jet->constscale_phi(), jet->constscale_E(), jet->constscale_m(),
                        jet->ActiveAreaPx(), jet->ActiveAreaPy(), jet->ActiveAreaPz(), jet->ActiveAreaE(),
                        evtShape->rhoKt4LC(),
                        eventInfo->averageIntPerXing(),
                        vertex->nTracks());

        if(!kIsData) local_truth_label = jet->flavor_truth_label();
        susyobj.ApplyJetSystematics(iJet, jet->constscale_eta(),
                                    local_truth_label,
                                    eventInfo->averageIntPerXing(),
                                    vertex->nTracks(), whichsyste);

        if(susyobj.IsGoodJet(iJet, jet->constscale_eta(),
                             jet->emfrac(), jet->hecf(), jet->LArQuality(), jet->HECQuality(), jet->AverageLArQF(),
                             jet->Timing(), jet->sumPtTrk_pv0_500MeV(), jet->fracSamplingMax(), jet->SamplingMax(), jet->NegativeE(),
                             eventInfo->RunNumber(), etcut, etacut, JetID::VeryLooseBad))
            goodjets.push_back(iJet);

        float jetPt = susyobj.GetJetTLV(iJet).Pt();
        if(jetPt>etcut && fabs(jet->eta())<etacut && jet->E()>0) {
            jet_index.push_back(iJet);
        }
    }
    return jet_index;
}
//----------------------------------------------------------
vector<int> get_jet_signal(D3PDReader::JetD3PDObject *jets, SUSYObjDef &susyobj, vector<int> jets_base, float ptcut, float etacut, float jvfcut) {
    vector<int> jets_signal;
    for(unsigned int ijet=0; ijet<jets_base.size(); ijet++) {
        int jetID = jets_base.at(ijet);
        float jetPt = susyobj.GetJetTLV(jetID).Pt();
        bool isSignalJet = (jetPt>ptcut &&
                            fabs(jets->eta()->at(jetID))<etacut &&
                            jets->jvtxf()->at(jetID)>jvfcut);

        if(isSignalJet) {
            jets_signal.push_back(jetID);
        }
    }
    return jets_signal;
}
//----------------------------------------------------------
bool check_jet_tileHotSpot(D3PDReader::JetD3PDObject* jet, vector<int>& baseline_jets, SUSYObjDef& susyObj,
                           bool isData, int run_number)
{
    if(isData)
    {
        for(unsigned int i=0; i < baseline_jets.size(); i++){
            int idx = baseline_jets.at(i);
            const TLorentzVector* jetLV = & susyObj.GetJetTLV(idx);

            if(susyObj.isHotTile(run_number, jet->fracSamplingMax()->at(idx),
                                 jet->SamplingMax()->at(idx), jetLV->Eta(), jetLV->Phi())) return true;
        }
    }
    return false;
}
//----------------------------------------------------------
vector<int> get_photons_baseline(D3PDReader::PhotonD3PDObject *photons, SUSYObjDef &susyObj,
				 float ptcut, float etacut, SystErr::Syste whichsyste, int Quality)
{
    vector<int> photons_base;
    unsigned int ph_Quality; // DG 2014-05-30 why is this not initialized? are we ok with '0'?
    int CutEnabButFailTight=0;
    for(int iph=0; iph<photons->n(); iph++){
        const D3PDReader::PhotonD3PDObjectElement *ph = & (*photons)[iph];
        if( susyObj.FillPhoton(iph,
                               ph->eta(), ph->phi(),
                               ph->cl_eta(), ph->cl_phi(), ph->cl_E(), ph->etas2(),
                               ph->isEM(), ph->OQ(), ph->isConv(),
                               ph->reta(), ph->rphi(), ph->Ethad1(), ph->Ethad(),
                               ph->E277(), ph->weta2(), ph->f1(),
                               ph->emaxs1(), ph->Emax2(), ph->Emins1(),
                               ph->fside(), ph->wstot(), ph->ws3(),
                               ph_Quality, CutEnabButFailTight, ptcut, etacut,
                               whichsyste, Quality)) {
            photons_base.push_back(iph);
        }
    } // end loop over photons
    return photons_base;
}
//----------------------------------------------------------
vector<int> get_photons_signal(D3PDReader::PhotonD3PDObject *photons, vector<int> photons_base,
			       SUSYObjDef &susyObj, int nPV, float ptcut,
                               float isocut, unsigned int isoType)
{
    vector<int> photons_signal;

    for(unsigned int iph=0; iph<photons_base.size(); iph++){

        int idx = photons_base.at(iph);
        if( susyObj.IsSignalPhotonCompIso(idx,
				          photons->ED_median()->at(idx),
				          photons->Etcone40()->at(idx),
				          photons->cl_E()->at(idx),
				          photons->cl_eta()->at(idx),
				          photons->etas2()->at(idx),
				          photons->etap()->at(idx),
				          photons->isConv()->at(idx),
				          photons->Etcone40_corrected()->at(idx),
				          nPV,
				          ptcut,
				          isocut,
				          isoType) )
        {
	    // Save Photon
	    photons_signal.push_back(idx);
        }

    } // end loop over base photons

    return photons_signal;
}
//----------------------------------------------------------
bool PassesPtNCut(vector<int> signal_muons, vector<int> signal_electrons, SUSYObjDef &susyObj, int n, float cutval){

    vector<float> signalLeptonPt;
    for(unsigned int i=0;i<signal_electrons.size();i++){
        signalLeptonPt.push_back(susyObj.GetElecTLV(signal_electrons[i]).Pt());
    }
    for(unsigned int i=0;i<signal_muons.size();i++){
        signalLeptonPt.push_back(susyObj.GetMuonTLV(signal_muons[i]).Pt());
    }

    // Sort the pt's
    sort(signalLeptonPt.begin(),signalLeptonPt.end(),greater<float>());

    if((int)signalLeptonPt.size()<n) return false;
    if(signalLeptonPt[n-1]<cutval) return false;

    return true;
}
//----------------------------------------------------------
bool IsBadLooseMinus(vector<int> jets_baseline, D3PDReader::JetD3PDObject *jets, SUSYObjDef &susyobj) {

    for(unsigned int ijet=0; ijet<jets_baseline.size(); ijet++) {
        int jetID = jets_baseline.at(ijet);
        float jetPt = susyobj.GetJetTLV(jetID).Pt();
        bool isBadJet = JetID::isBadJet(JetID::VeryLooseBad,
                                        jets->emfrac()->at(jetID),
                                        jets->hecf()->at(jetID),
                                        jets->LArQuality()->at(jetID),
                                        jets->HECQuality()->at(jetID),
                                        jets->Timing()->at(jetID),
                                        jets->sumPtTrk_pv0_500MeV()->at(jetID)/1000., //GeV
                                        jets->emscale_eta()->at(jetID),
                                        jetPt/1000., //(jets->pt()->at(jetID))/1000., //jetPt/1000., //(jets->pt()->at(jetID))/1000., // GeV
                                        jets->fracSamplingMax()->at(jetID),
                                        jets->NegativeE()->at(jetID), // MeV
                                        jets->AverageLArQF()->at(jetID));
        if(isBadJet) return true;
    }
    return false;
}
//----------------------------------------------------------
bool IsBadJetEvent(D3PDReader::JetD3PDObject *jets, vector<int> jets_baseline, float jet_base_ptcut, SUSYObjDef &susyObj) {
    // This is poorly written.  It loops over jets and then calls a function which loops over the jets again...
    // Technically it still works but that's mostly just luck.  The logic is misguided.
    for(unsigned int iJet=0; iJet<jets_baseline.size(); iJet++) {
        int jetID = jets_baseline.at(iJet);
        bool isbadlooseminus = false;
        isbadlooseminus = IsBadLooseMinus(jets_baseline,jets,susyObj);
        float jetpT = susyObj.GetJetTLV(jetID).Pt(); //jets->pt()->at(jetID)
        if(isbadlooseminus && (jetpT > jet_base_ptcut)) {
            return true;
        }
    }
    return false;
}
//----------------------------------------------------------
bool IsBadMuonEvent(SUSYObjDef &susyObj, D3PDReader::MuonD3PDObject *muon, vector<int> baseline_muons, float qoverpcut) {
    for(unsigned int i=0; i<baseline_muons.size(); i++) {
        if(susyObj.IsBadMuon(muon->qoverp_exPV()->at(baseline_muons[i]), muon->cov_qoverp_exPV()->at(baseline_muons[i]), qoverpcut)) {
            return true;
        }
    }
    return false;
}
//----------------------------------------------------------
bool IsCosmic(SUSYObjDef &susyObj, D3PDReader::MuonD3PDObject *muon, vector<int> baseline_muons, float z0cut, float d0cut) {
    for(unsigned int i=0; i<baseline_muons.size(); i++) {
        const D3PDReader::MuonD3PDObjectElement* muonElement = & (*muon)[baseline_muons[i]];
        //if(susyObj.IsCosmicMuon(muon->z0_exPV()->at(baseline_muons[i]),muon->d0_exPV()->at(baseline_muons[i]), z0cut, d0cut))
        if(susyObj.IsCosmicMuon(muonElement->z0_exPV(), muonElement->d0_exPV(), z0cut, d0cut)) {
            return true;
        }
    }
    return false;
}
//----------------------------------------------------------
bool PrimaryVertexCut(SUSYObjDef &susyObj, D3PDReader::PrimaryVertexD3PDObject *vertex) {
    if(!susyObj.IsGoodVertex(vertex->nTracks())) {
        //cout<<"Bad Vertex"<<endl;
        return false;
    }
    return true;
}
//----------------------------------------------------------
