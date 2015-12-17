// version 1.1, 21/09/2015
#include "SusyCommon/ss3l_chargeflip.h"
#include "SusyNtuple/ss3l_chargeflip_types.h"

#include "xAODEgamma/Electron.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h" // DG-2015-12-01 not used but needed to compile

#include <vector>

using std::vector;

//----------------------------------------------------------
double deltaR(double eta1,double phi1,double eta2,double phi2)
{
    double deltaEta = eta2-eta1;
    double deltaPhi = fabs(phi2-phi1);
    while(deltaPhi>TMath::Pi()) deltaPhi-=2*TMath::Pi();
    return sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);
}
//----------------------------------------------------------
void fillElectronChargeFlip( xAOD::ElectronContainer* reco_electrons,
                             const xAOD::TruthParticleContainer* container_particles,
                             int mcChannelNumber)
{
    //El_chargeFlip.clear();

    // First recover truth electrons
    //  const xAOD::EventInfo* eventInfo = nullptr;
    //  xEvent.retrieve( eventInfo, "EventInfo").ignore();
    int status_cut = -1;
    // switch(eventInfo->mcChannelNumber())
    switch(mcChannelNumber)
	{
	case 361100: case 361101: case 361102: // PowhegPythia W+ -> lnu
	case 361103: case 361104: case 361105: // PowhegPythia W- -> lnu
	case 361106: case 361107: case 361108: // PowhegPythia Z -> ll
	    status_cut = 62;
	    break;
	default:
	    status_cut = 2;
	};

    // const xAOD::TruthParticleContainer *container_particles = nullptr;
    // xEvent.retrieve(container_particles,"TruthParticles").ignore();
    vector<const xAOD::TruthParticle*> prompt_electrons;
    for(auto parent : *container_particles)
	{
	    int pid = parent->absPdgId();
	    if(((pid<23 || pid>37) && (pid<1000001 || pid>2999999)) || parent->p4().M()<12000.) continue;
	    if(status_cut!=-1 && parent->status()!=status_cut) continue;
	    if(!parent->hasDecayVtx()) continue;
	    auto vtx = parent->decayVtx();
	    const int ndec = vtx->nOutgoingParticles();
	    if(!ndec) continue;
	    for(int j=0;j<ndec;++j)
		{
		    auto child = vtx->outgoingParticle(j);
		    int cid = child->absPdgId();
		    if(cid==pid) break;
		    while(cid==15) // catch tau-mediated prompt electrons
			{
			    if(!child->hasDecayVtx()) break;
			    auto tauvtx = child->decayVtx();
			    const int ntaudec = tauvtx->nOutgoingParticles();
			    for(int k=0;k<ntaudec;++k)
				{
				    auto tauchild = tauvtx->outgoingParticle(k);
				    cid = tauchild->absPdgId();
				    if(cid==11 || cid==15)
					{
					    child = tauchild;
					    break;
					}
				}
			}
		    if(cid==11)
			{
			    prompt_electrons.push_back(child);
			}
		}
	}

    // Then try to match those to reconstructed electrons
    //  El_chargeFlip.resize(reco_electrons.size(),CHARGE_UNKNOWN);
    //  for(int j=0;j<reco_electrons.size();++j)
    //  {    //   auto reco_electron = electrons[j];
    for(auto reco_electron : *reco_electrons)
	{
	    int nMatched = 0, sumMatchedQ = 0;
	    for(auto prompt_electron : prompt_electrons)
		{
		    double dr = deltaR(reco_electron->eta(),reco_electron->phi(),prompt_electron->eta(),prompt_electron->phi());
		    if(dr<0.1)
			{
			    nMatched++;
			    sumMatchedQ += prompt_electron->charge();
			}
		    else if(dr<0.2 &&
			    //   e->auxdata< int >("truthType")==4 &&
			    //   e->auxdata< int >("truthOrigin")==5)
			    reco_electron->auxdata< int >("truthType")==4 &&
			    reco_electron->auxdata< int >("truthOrigin")==5)
			{
			    //  int originbkg = e->auxdata< int >("bkgTruthOrigin");
			    int originbkg = reco_electron->auxdata< int >("bkgTruthOrigin");
			    if(originbkg==10 || (originbkg>=12 && originbkg<=15) || originbkg==22)
				{
				    nMatched++;
				    sumMatchedQ += prompt_electron->charge();
				}
			}
		}
	    if(nMatched==1)
		{
		    if(reco_electron->charge()==sumMatchedQ)
			// El_chargeFlip[j] = CHARGE_CORRECT;
			reco_electron->auxdata< int >("chargeFlip") = CHARGE_CORRECT;
		    else
			// El_chargeFlip[j] = CHARGE_FLIPPED;
			reco_electron->auxdata< int >("chargeFlip") = CHARGE_FLIPPED;
		}
	    else if(nMatched>1)
		{
		    int sumR = reco_electron->charge()*nMatched;
		    if(sumR==sumMatchedQ)
			//El_chargeFlip[j] = CHARGE_MAYBE_CORRECT;
			reco_electron->auxdata< int >("chargeFlip") = CHARGE_MAYBE_CORRECT;
		    else if(sumR==-sumMatchedQ)
			// El_chargeFlip[j] = CHARGE_MAYBE_FLIPPED;
			reco_electron->auxdata< int >("chargeFlip") = CHARGE_MAYBE_FLIPPED;
		    else
			// El_chargeFlip[j] = CHARGE_AMBIGUOUS;
			reco_electron->auxdata< int >("chargeFlip") = CHARGE_AMBIGUOUS;
		}
	    else //El_chargeFlip[j] = CHARGE_UNKNOWN;
		reco_electron->auxdata< int >("chargeFlip") = CHARGE_UNKNOWN;


	    if(reco_electron->auxdata< int >("chargeFlip") == CHARGE_UNKNOWN) // El_chargeFlip[j]==CHARGE_UNKNOWN)
		{
		    int type = reco_electron->auxdata< int >("truthType");
		    //int origin = e->auxdata< int >("truthOrigin");
		    //int originbkg = e->auxdata< int >("bkgTruthOrigin");
		    int mpid = reco_electron->auxdata< int >("bkgMotherPdgId");
		    if(type==4 && abs(mpid)==11) // origin==5 && originbkg==40
			{
			    if(reco_electron->charge()==(mpid==-11))
				//El_chargeFlip[j] = CHARGE_MCTRUTHCLASSIFIER_CORRECT;
				reco_electron->auxdata< int >("chargeFlip") = CHARGE_MCTRUTHCLASSIFIER_CORRECT;
			    else
				// El_chargeFlip[j] = CHARGE_MCTRUTHCLASSIFIER_FLIPPED;
				reco_electron->auxdata< int >("chargeFlip") = CHARGE_MCTRUTHCLASSIFIER_FLIPPED;
			}
		}
	}
}
//----------------------------------------------------------
