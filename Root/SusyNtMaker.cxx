#include "SusyCommon/SusyNtMaker.h"

//SusyCommon
#include "SusyCommon/SusyObjId.h"
#include "SusyCommon/ss3l_chargeflip.h"

//SusyNtuple
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/TriggerTools.h"

//xAOD
//#include "EventPrimitives/EventPrimitivesHelper.h"
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
    m_file_outtree(0),
    m_outtree(0),
    m_susyNt(0)
{

}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::Init(TTree *tree)
{
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::SlaveBegin(TTree* tree)
{
    cout << "SusyNtMaker::SlaveBegin" << endl;
    XaodAnalysis::SlaveBegin(tree);
    XaodAnalysis::Init(tree);

    return;
}
//////////////////////////////////////////////////////////////////////////////
Bool_t SusyNtMaker::Process(Long64_t entry)
{
    cout << "SusyNtMaker::Process" << endl;

    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void SusyNtMaker::Terminate()
{
    cout << "SusyNtMaker::Terminate" << endl;

    return;
}
