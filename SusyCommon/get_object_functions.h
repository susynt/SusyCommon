// Dear emacs, this is -*- c++ -*-
/**

Functions to select objects from D3PDs

These functions were originally in MultiLep -> ElectronTools, MuonTools, TauTools

author: Steve Farrell <sfarrell@cern.ch>

Migrated from MultiLep to SusyCommon, May 2014
davide.gerbaudo@gmail.com

\todo : enclose them in a namespace
\todo : get rid of GeV preprocessor const

*/

#include "SUSYTools/SUSYObjDef.h"

#include "TLorentzVector.h"
#include "TVector2.h"

#include <vector>

#ifndef get_object_functions_h
#define get_object_functions_h


namespace D3PDReader
{
class TauD3PDObject;
class ElectronD3PDObject;
class EventInfoD3PDObject;
class EFTriggerD3PDObject;
class MuonD3PDObject;
class JetD3PDObject;
class PhotonD3PDObject;
class PrimaryVertexD3PDObject;
class TauD3PDObject;
}


#define GeV 1000.

// Select baseline taus
vector<int> get_taus_baseline(D3PDReader::TauD3PDObject* taus, SUSYObjDef& susyObj,
                              float ptCut=20.*GeV, float etaCut=2.47,
                              SUSYTau::IDLevel jetBDTLevel=SUSYTau::TauLoose,
                              SUSYTau::IDLevel eleBDTLevel=SUSYTau::TauLoose,
                              SUSYTau::IDLevel muonLevel=SUSYTau::TauLoose,
                              SystErr::Syste sys=SystErr::NONE,
                              bool isD3PD1512=false);

// select signal taus
vector<int> get_taus_signal(D3PDReader::TauD3PDObject* taus, vector<int>& taus_base,
                            SUSYObjDef& susyObj, float ptCut=20.*GeV, float etaCut=2.47,
                            SUSYTau::IDLevel jetBDTLevel=SUSYTau::TauMedium,
                            SUSYTau::IDLevel eleBDTLevel=SUSYTau::TauMedium,
                            SUSYTau::IDLevel muonLevel=SUSYTau::TauMedium,
                            SystErr::Syste sys=SystErr::NONE,
                            bool isD3PD1512=false);


// met electrons
vector<int> get_electrons_met(D3PDReader::ElectronD3PDObject *electrons, SUSYObjDef &susyobj);

// baseline electrons
vector<int> get_electrons_baseline(D3PDReader::ElectronD3PDObject* electrons, bool kIsData, int run_number, SUSYObjDef& susyobj, float etcut,
				   float etacut, SystErr::Syste el_syst);

// signal electrons
vector<int> get_electrons_signal(D3PDReader::ElectronD3PDObject *electrons, vector<int> electrons_base,
                                 D3PDReader::MuonD3PDObject *muons, vector<int> muons_base, unsigned int nGoodVtx,
                                 bool isData, SUSYObjDef &susyobj, float etcut, float ptconeRelCut=0.16, float etconeRelCut=0.18,
                                 float d0Sig_cutval=5., float z0SinTheta_cutval=0.4, bool removeLeps = true);


// signal electrons in trigger plateau
vector<int> get_electrons_signal_triggerplateau(D3PDReader::ElectronD3PDObject *electrons, vector<int> electrons_signal, SUSYObjDef &susyobj,
                                                float ptcut_triggerplateau);

// electron trigger matching function
bool check_electron_triggermatch(D3PDReader::ElectronD3PDObject *electrons, int eItr, D3PDReader::EFTriggerD3PDObject *efmuontrigger, vector<int>* trigger,
                                 float DeltaRMin);

// function to get information about muon electron (fired or not and the correct plateauPtThreshold)
bool getElectronTriggerInformation(int RunNumber, float &plateauPtThreshold1L, float &plateauPtThreshold2L,
                                   bool &triggered1L, bool &triggered2L, bool &triggered2L_2,
                                   vector<int> **trig_EF_el_EF, vector<int> **trig_EF_el_EF_2e, D3PDReader::EFTriggerD3PDObject *efeltrigger,
                                   bool debug, bool kIs2L);

bool getElectronMuonTriggerInformation(int RunNumber, float &plateauPtThreshold1E1M_e, float &plateauPtThreshold1E1M_mu,
                                       bool &triggered1E1M, bool &triggered1E1M_2,
                                       D3PDReader::EFTriggerD3PDObject *eftrigger, vector<int> **trig_EF_el_EF_emu, vector<int> **trig_EF_mu_EF_emu,
                                       bool debug);



// all muons
vector<int> get_muons_all(D3PDReader::MuonD3PDObject *muons);

// baseline muons
vector<int> get_muons_baseline(D3PDReader::MuonD3PDObject *muons, bool kIsData, SUSYObjDef &susyobj,
                               float etcut, float etacut, SystErr::Syste mu_syst);

// signal muons
vector<int> get_muons_signal(D3PDReader::MuonD3PDObject *muons, vector<int> muons_base,
                             D3PDReader::ElectronD3PDObject *electrons,  vector<int> electrons_base,
                             unsigned int nGoodVtx, bool isData, SUSYObjDef &susyobj, float ptcut,
                             float ptconeRelCut=0.12, float d0Sig_cutval=3., float z0SinTheta_cutval=1.,
                             bool removeLeps=true);


vector<int> get_jet_baseline(D3PDReader::JetD3PDObject *jet, D3PDReader::PrimaryVertexD3PDObject *vertex, D3PDReader::EventInfoD3PDObject *eventinfo, 
                             bool kIsData, SUSYObjDef &susyobj, float etcut, float etacut, SystErr::Syste whichsyste, bool nosmear, vector<int> &goodjets);

vector<int> get_jet_signal(D3PDReader::JetD3PDObject *jet, SUSYObjDef &susyobj, vector<int> jet_base, float ptcut, float etacut, float jvfcut);

// Temporary, until there is a better SUSYTools function for this
bool check_jet_tileHotSpot(D3PDReader::JetD3PDObject* jet, vector<int>& baseline_jets, SUSYObjDef& susyObj, bool isData, int run_number);

// Baseline Photons
vector<int> get_photons_baseline(D3PDReader::PhotonD3PDObject *photons, SUSYObjDef &susyObj,
                                 float ptcut, float etacut, SystErr::Syste whichsyste, int Quality);

// Signal Photons
vector<int> get_photons_signal(D3PDReader::PhotonD3PDObject *photons, vector<int> photons_base, 
                               SUSYObjDef &susyObj, int nPV, float ptcut = 130000., 
                               float iso = 4000., unsigned int isoType = 1);

//----------------------------------------------------------
// A few additional utility functions
//----------------------------------------------------------

// define function to calculate MET 2D vector; Now you can control met container via metFlavor
TVector2 GetMetVector(SUSYObjDef &susyobj, D3PDReader::JetD3PDObject *jet, 
                      D3PDReader::MuonD3PDObject *muon, D3PDReader::ElectronD3PDObject *electron, 
                      D3PDReader::METD3PDObject *met, D3PDReader::EventInfoD3PDObject *eventInfo, 
                      vector<int> baseline_muons, vector<int> baseline_electrons, 
                      vector<int> all_electrons, SystErr::Syste whichsyste, 
                      SUSYMet::met_definition metFlavor = SUSYMet::Default, 
                      bool doMuonElossCorrection = false, bool doEgammaJetFix = false);

bool PassesPtNCut(vector<int> signal_muons, vector<int> signal_electrons, SUSYObjDef &susyobj, int n, float cutval);
bool IsBadJetEvent(D3PDReader::JetD3PDObject* jet, std::vector< int > baseline_jet, float jet_base_ptcut, SUSYObjDef &susyobj);
bool IsCosmic(SUSYObjDef &susyobj, D3PDReader::MuonD3PDObject *muon, vector<int> baseline_muons, float z0cut, float d0cut);
bool IsBadMuonEvent(SUSYObjDef &susyObj, D3PDReader::MuonD3PDObject *muon, vector<int> baseline_muons, float qoverpcut);
bool PrimaryVertexCut(SUSYObjDef &susyobj, D3PDReader::PrimaryVertexD3PDObject *vertex);
// 
template<class D3PDObject1, class D3PDObject2> 
vector<int> overlap_removal(SUSYObjDef &susyobj, D3PDObject1 *o1, vector<int> indices_1, 
                            D3PDObject2 *o2, vector<int> indices_2, float dr, 
                            bool sameType, bool removeSoft) 
{
    vector<int> survivors;
    for(unsigned int i=0; i<indices_1.size(); i++) {
        bool is_overlap = false;
        for(unsigned int j=0; j<indices_2.size(); j++) {
            TLorentzVector object1;
            TLorentzVector object2;
            // Fill the TLV with the corrected values. Different functions for each type of particle.
            if(dynamic_cast<D3PDReader::ElectronD3PDObject *>( o1 ))  object1 = susyobj.GetElecTLV(indices_1[i]);
            else if(dynamic_cast<D3PDReader::MuonD3PDObject *>( o1 )) object1 = susyobj.GetMuonTLV(indices_1[i]);
            else if(dynamic_cast<D3PDReader::JetD3PDObject *>( o1 ))  object1 = susyobj.GetJetTLV(indices_1[i]);
            else object1.SetPtEtaPhiM(o1->pt()->at(indices_1[i]), o1->eta()->at(indices_1[i]), 
                                      o1->phi()->at(indices_1[i]), o1->m()->at(indices_1[i]));
            if(dynamic_cast<D3PDReader::ElectronD3PDObject *>( o2 ))  object2 = susyobj.GetElecTLV(indices_2[j]);
            else if(dynamic_cast<D3PDReader::MuonD3PDObject *>( o2 )) object2 = susyobj.GetMuonTLV(indices_2[j]);
            else if(dynamic_cast<D3PDReader::JetD3PDObject *>( o2 ))  object2 = susyobj.GetJetTLV(indices_2[j]);
            else object2.SetPtEtaPhiM(o2->pt()->at(indices_2[j]), o2->eta()->at(indices_2[j]), 
                                      o2->phi()->at(indices_2[j]), o2->m()->at(indices_2[j]));

            // Don't remove an object against itself (for electron-electron overlap) 
            if(sameType && i==j) continue;
            if (removeSoft) {
                if ( (object1.E()/cosh(object1.Eta())) > (object2.E()/cosh(object2.Eta()))) {
                    continue;    // remove lowest Et
                }
                if(object1.DeltaR(object2) <= dr) {
                    is_overlap = true;
                }
            } else {

                if(object1.DeltaR(object2) <= dr) {
                    is_overlap = true;
                }
            }
        }
        if(is_overlap) {
            continue;
        }
        survivors.push_back(indices_1[i]);
    }
    return survivors;
}


// A template to Remove the SFOS pair below a certain cut
// can be called with muons and electrons
template<class particle> vector<int> RemoveSFOSPair(SUSYObjDef &susyobj, particle *particles, vector<int> baseline_particle, float cut) {
    vector<int> baseline_withoutSFOS;
    vector<unsigned int> to_remove;
    TLorentzVector particle_1;
    TLorentzVector particle_2;
    for(unsigned int i = 0; i < baseline_particle.size(); i++) {
        for(unsigned int j = i+1; j < baseline_particle.size(); j++) {
            if( particles->charge()->at(baseline_particle[i]) * particles->charge()->at(baseline_particle[j]) < 0 ) {
                // 4-vectors.. Should be smeared and corrected
                if (dynamic_cast<D3PDReader::ElectronD3PDObject *>( particles )) {
                    particle_1 = susyobj.GetElecTLV(baseline_particle.at(i));
                    particle_2 = susyobj.GetElecTLV(baseline_particle.at(j));
                }
                else if (dynamic_cast<D3PDReader::MuonD3PDObject *>( particles )) {
                    particle_1 = susyobj.GetMuonTLV(baseline_particle.at(i));
                    particle_2 = susyobj.GetMuonTLV(baseline_particle.at(j));
                }
                // For taus, or other default particles, use the d3pd variables to calculate Mll
                else{
                    particle_1.SetPtEtaPhiM(particles->pt()->at(baseline_particle[i]), particles->eta()->at(baseline_particle[i]), 
                                            particles->phi()->at(baseline_particle[i]), particles->m()->at(baseline_particle[i]));
                    particle_2.SetPtEtaPhiM(particles->pt()->at(baseline_particle[j]), particles->eta()->at(baseline_particle[j]), 
                                            particles->phi()->at(baseline_particle[j]), particles->m()->at(baseline_particle[j]));
                }
                if(sqrt((particle_1+particle_2)*(particle_1+particle_2))<cut) {
                    // flag both indices of the pair for removal
                    to_remove.push_back(i);
                    to_remove.push_back(j);
                }
            }
        }
    }
    // write the ones that are not to remove in the vector without the pairs
    for(unsigned int i=0; i<baseline_particle.size(); i++) {
        bool isInPair = false;
        for(unsigned int j=0; j<to_remove.size(); j++) {
            if(i==to_remove[j]) {
                isInPair=true;
            }
        }
        if(!isInPair) {
            baseline_withoutSFOS.push_back(baseline_particle[i]);
        }
    }
    return baseline_withoutSFOS;
}

#undef GeV

#endif
