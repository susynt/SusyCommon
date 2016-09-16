
// SusyCommon
#include "SusyCommon/XaodAnalysis.h"
//#include "SusyCommon/get_object_functions.h"

// #include "egammaAnalysisUtils/egammaTriggerMatching.h"
// #include "D3PDReader/JetD3PDObject.h"
#include "xAODBase/IParticleHelpers.h" // setOriginalObjectLink
#include "xAODEgamma/EgammaxAODHelpers.h"

#include "SusyNtuple/RecoTruthClassification.h"

#include "TChainElement.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TSystem.h"

#include <limits>
#include <algorithm> // copy_if, transform
#include <iterator> // back_inserter
#include <numeric> // accumulate

#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"




using namespace std;

using Susy::XaodAnalysis;


#undef CHECK
#define CHECK( ARG )                                                \
    do {                                                            \
        const bool result = ARG;                                    \
        if( ! result ) {                                            \
            ::Error( "XaodAnalysis", "Failed to execute: \"%s\"",   \
                     #ARG );                                        \
            exit(-1);                                               \
        }                                                           \
    } while( false )


//----------------------------------------------------------
XaodAnalysis::XaodAnalysis() :
    m_input_chain(0),
    //m_sample(""),
    m_triggerSet("run2"),
    m_stream(Stream_Unknown),
    m_isDerivation(false), // dantrim event shape
    m_nEventsProcessed(0),
    m_sumOfWeights(0),
    m_sumOfWeightsSquared(0),
    m_isAF2(false),
    m_is8TeV(true),
    m_selectPhotons(false),
    m_selectTaus(false),
    m_selectTruth(false),
    m_doMetMuCorr(false),
    m_doMetFix(false),
    m_lumi(LUMI_A_A4),
    m_mcRun(0),
    m_mcLB(0),
    m_sys(false),
//        m_eleMediumSFTool(0),
//        m_pileup(0),
//        m_pileup_up(0),
//        m_pileup_dn(0),
//        m_susyXsec(0),
//        m_hforTool(),
    m_grl(NULL),
    m_tree(NULL),
    m_entry(0),
    m_dbg(0),
    m_isMC(false),
    m_isData15(false),
    m_isData16(false),
    m_isMC15b(false),
    m_isMC15c(false),
    m_flagsAreConsistent(false),
    m_flagsHaveBeenChecked(false),
    m_event(xAOD::TEvent::kClassAccess),
    //m_event(xAOD::TEvent::kAthenaAccess), ///> dantrim April 28 2016 -- to get CutBookkeepers needs kAthenaAccess or kClassAcess
    m_store(),
    m_eleIDDefault(eleTightLLH),
	m_electronEfficiencySFTool(0),
    m_elecSelLikelihoodVeryLoose(0),
    m_elecSelLikelihoodLoose(0),
    m_elecSelLikelihoodLooseBLayer(0),
    m_elecSelLikelihoodMedium(0),
    m_elecSelLikelihoodTight(0),
    m_photonSelLoose(0),
    m_photonSelTight(0),
	m_pileupReweightingTool(0),
    m_muonSelectionToolVeryLoose(0),
    m_muonSelectionToolLoose(0),
    m_muonSelectionToolMedium(0),
    m_muonSelectionToolTight(0),
    m_isoToolGradientLooseTight(0),
    m_isoToolGradientTightCalo(0),
    m_isoToolLooseTrackOnlyLoose(0),
    m_isoToolLoose(0),
    m_isoToolTight(0),
	m_tauTruthMatchingTool(0),
	m_tauTruthTrackMatchingTool(0),
    m_tauSelToolLoose(0),
    m_tauSelToolMedium(0),
    m_tauSelToolTight(0),
    m_polreweight(nullptr),
    m_evtTrigBits(m_nTriggerBits)
{
    clearOutputObjects();
    clearContainerPointers();


}
//----------------------------------------------------------
void XaodAnalysis::Init(TTree *tree)
{
    bool verbose = m_dbg>0;
    xAOD::Init("Susy::XaodAnalysis").ignore();

    m_isMC = XaodAnalysis::isSimuFromSamplename(m_inputContainerName);

    if(!m_isMC && !(m_isData15 || m_isData16)) {
        cout << "XaodAnalysis::Init    ERROR Inconsistent options! This sample ("
             << m_inputContainerName << ") is flagged as NOT MC but neither the "
             << "data15 nor data16 flags are true!" << endl;
        exit(1);
    }

    if(!m_isMC && m_isMC15b) { 
        m_isMC15b = false;
    }
    if(!m_isMC && m_isMC15c) {
        m_isMC15c = false;
    }

    if(m_isMC){
        // get the inital (pre-skimmed) counters
        TObjArray* chainFiles = m_input_chain->GetListOfFiles();
        TIter next(chainFiles);
        TChainElement *chFile=0;
            while (( chFile=(TChainElement*)next() )) {
                cout << "XaodAnalysis::Init    CutBookkeeper info for: " << chFile->GetTitle() << endl;
                TFile* f = TFile::Open(chFile->GetTitle());
                m_event.readFrom(f);
                m_event.getEntry(0);
                if(!XaodAnalysis::getCutBookkeeperInfo(m_event)) exit(1);
                f->Close();
                f->Delete();
            } // while
        cout << "XaodAnalysis::Init    CutBookkeeper info totals: " << endl;
        cout << "    > m_nEventsProcessed   : " << m_nEventsProcessed << endl;
        cout << "    > m_sumOfWeights       : " << m_sumOfWeights << endl;
        cout << "    > m_sumOfWeightsSquared: " << m_sumOfWeightsSquared << endl;
    } // if MC

    m_event.readFrom(tree);
    m_isDerivation = XaodAnalysis::isDerivationFromMetaData(tree, verbose);
    m_stream = XaodAnalysis::streamFromSamplename(m_inputContainerName, m_isMC);

    // get the directory for the data/
    char *tmparea=getenv("ROOTCOREBIN");
    char* TestArea = getenv("TestArea");

    if (tmparea != NULL) {
        m_data_dir = tmparea;
        m_data_dir = m_data_dir + "/data/";
    }
    else if (TestArea != NULL ) {/// Athena
        tmparea = TestArea;
        m_data_dir = tmparea;
        m_data_dir = m_data_dir + "/";
    } else {
        cout << " RootCore area not set up " << endl
             <<"Exiting... "<<endl << endl;
        exit(1);
    }

    // GRL tool (if data)
    if(!m_isMC) {
        if(!initGrlTool()) {
            cout << "XaodAnalysis::Init ERROR    Unable to initialize GRL tool. Exiting." << endl;
            exit(1);
        }
    }

    cout << "----------------------------------------------------------" << endl;
    cout << "XaodAnalysis::Init    Treating sample as " << (m_isMC ? ( ( m_isMC15c ? "mc15c" : (m_isMC15b ? "mc15b" : "mc15a") )) : ( m_isData16 ? "data16_13TeV" : "data15_13TeV" ) ) << endl;
    cout << "----------------------------------------------------------" << endl;

    // initialize SUSYTools
    initSusyTools();
    // initialize other analysis tools
    initLocalTools();
    // grab the systematics list for our tools
    if(m_isMC && m_sys) getSystematicList();
    else{
        ST::SystInfo infodef;
        infodef.affectsKinematics = false;
        infodef.affectsWeights = false;
        infodef.affectsType = ST::Unknown;
        systInfoList.push_back(infodef);
    }
}
//----------------------------------------------------------
XaodAnalysis::~XaodAnalysis()
{
    cout<<"~XaodAnalysis : todo"<<endl;
}
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::SlaveBegin(TTree *tree)
{

    if(m_dbg) cout << "XaodAnalysis::SlaveBegin" << endl;

    if(m_isMC){
    // do nothing
    }
}


/*--------------------------------------------------------------------------------*/
// Main process loop function - This is just an example for testing
/*--------------------------------------------------------------------------------*/
Bool_t XaodAnalysis::Process(Long64_t entry)
{
    static Long64_t chainEntry = -1;
    chainEntry++;
    m_event.getEntry(entry);
    retrieveCollections();
    if(m_dbg || chainEntry%10000==0)
        {
            const xAOD::EventInfo* eventinfo = xaodEventInfo();
            cout<<"run "<<eventinfo->eventNumber()<<" event "<<eventinfo->runNumber()<<endl;
        }
    // Object selection
    // SusyNtSys sys = NtSys_NOM;
    // selectObjects(sys);
    // buildMet();
    clearOutputObjects();
    deleteShallowCopies();
    return kTRUE;
}

/*--------------------------------------------------------------------------------*/
// The Terminate() function is the last function to be called during
// a query. It always runs on the client, it can be used to present
// the results graphically or save the results to file.
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::Terminate()
{
    if(m_dbg) cout << "XaodAnalysis::Terminate" << endl;

    if(m_isMC){
        // delete m_susyXsec;
        // delete m_pileup;
        // delete m_pileup_up;
        // delete m_pileup_dn;
    }

    delete m_electronEfficiencySFTool;
    delete m_elecSelLikelihoodVeryLoose;
    delete m_elecSelLikelihoodLoose;
    delete m_elecSelLikelihoodLooseBLayer;
    delete m_elecSelLikelihoodMedium;
    delete m_elecSelLikelihoodTight;
    delete m_photonSelLoose;
    delete m_photonSelTight;

    delete m_pileupReweightingTool;
    delete m_muonSelectionToolVeryLoose;
    delete m_muonSelectionToolLoose;
    delete m_muonSelectionToolMedium;
    delete m_muonSelectionToolTight;
    delete m_tauTruthMatchingTool;
    delete m_TauEffEleTool;
    delete m_tauTruthTrackMatchingTool;
    delete m_tauSelToolLoose;
    delete m_tauSelToolMedium;
    delete m_tauSelToolTight;

    delete m_isoToolGradientLooseTight;
    delete m_isoToolGradientTightCalo;
    delete m_isoToolLooseTrackOnlyLoose;
    delete m_isoToolLoose;
    delete m_isoToolTight;
    delete m_polreweight;

    // dantrim trig
  //  delete m_trigTool;
  //  delete m_configTool;
  //  //delete m_trigEgammaMatchTool;
    
  // for(int i=TightLH; i<=LooseLH; i++){
    for(int i : Susy::leptonIds()){
        delete m_susyObj[i];
    }
}
//----------------------------------------------------------
XaodAnalysis& XaodAnalysis::initSusyTools()
{
    for(int susyObjId : Susy::leptonIds()){
        string idName = SusyObjId2str(static_cast<SusyObjId>(susyObjId));
        bool isEle = isEleObj(static_cast<SusyObjId>(susyObjId)); 
        string name = "SUSYObjDef_xAOD_" + idName;
        m_susyObj[susyObjId] = new ST::SUSYObjDef_xAOD(name);
        cout << "------------------------------------------------------------" << endl;
        cout << "XaodAnalysis::initSusyTools    " << name <<endl;
        cout << "------------------------------------------------------------" << endl;

        string config_type = isEle ? "ele" : "muo";
        string config_name = "SusyCommon/SUSYTools_SusyNt_" + config_type + idName + ".conf";
        cout << "XaodAnalysis::initSusyTools    Using SUSYTools configuration file: " << config_name << endl; 
        //CHECK( m_susyObj[susyObjId]->setProperty("ConfigFile", "SusyCommon/SUSYTools_SusyNt.conf") );
        CHECK( m_susyObj[susyObjId]->setProperty("ConfigFile", config_name) );
        
        // set the verbosity level of SUSYTools
        m_susyObj[susyObjId]->msg().setLevel(m_dbg ? MSG::DEBUG : MSG::WARNING);
    
        // set "datasource" for determining run configuration in SUSYTools
        ST::ISUSYObjDef_xAODTool::DataSource datasource = !m_isMC ? ST::ISUSYObjDef_xAODTool::Data : ( m_isAF2 ? ST::ISUSYObjDef_xAODTool::AtlfastII : ST::ISUSYObjDef_xAODTool::FullSim);
    
        //ST::SettingDataSource datasource = !m_isMC ? ST::Data : (m_isAF2 ? ST::AtlfastII : ST::FullSim);
        CHECK( m_susyObj[susyObjId]->setProperty("DataSource",datasource) );

        ///////////////////////////////////////
        // set up the pileup reweighting tool
        // inside of ST
        ///////////////////////////////////////
        // prw config files
        std::vector<std::string> prwFiles;
        if(!m_isMC15b && !m_isMC15c) {
            prwFiles.push_back("dev/SUSYTools/merged_prw.root");
        }
        else if(m_isMC15b && !m_isMC15c) {
            prwFiles.push_back("dev/SUSYTools/merged_prw_mc15b.root");
        }
        else if(!m_isMC15b && m_isMC15c) {
            prwFiles.push_back("dev/SUSYTools/merged_prw_mc15c_July27_afterFix.root");
        }
        else {
            cout << "XaodAnalysis::initSusyTools    "
                 << "Inconsisent MC options when selecting PRW config! Exitting." << endl;
            exit(1);
        }
        m_susyObj[susyObjId]->setProperty("PRWConfigFiles", prwFiles);
        // data luminosity profile
        std::vector<std::string> lumicalcFiles;
        lumicalcFiles.push_back(m_data_dir+"SusyCommon/ilumicalc_histograms_None_297730-306451_OflLumi-13TeV-005.root"); // data16: updated Sept 15 2016 with GRL v82-pro20-12 (18.857 fb-1)
        lumicalcFiles.push_back(m_data_dir+"SusyCommon/ilumicalc_histograms_None_276262-284484.root"); // data15: updated with GRL v75
        m_susyObj[susyObjId]->setProperty("PRWLumiCalcFiles", lumicalcFiles); 

        cout << "XaodAnalysis::initSusyTools    Configuring SUSYTools' PRW tool: " << endl;
        cout << "XaodAnalysis::initSusyTools     + prw config files..." << endl;
        for(int i = 0; i < (int)prwFiles.size(); i++) {
            cout << "XaodAnalysis::initSusyTools        > " << prwFiles[i] << endl;
        }
        cout << "XaodAnalysis::initSusyTools     + lumi calc files..." << endl;
        for(int i = 0; i < (int)lumicalcFiles.size(); i++) {
            cout << "XaodAnalysis::initSusyTools        > " << lumicalcFiles[i] << endl;
        }


        if(m_susyObj[susyObjId]->initialize() != StatusCode::SUCCESS){
            cout << "XaodAnalysis: Cannot intialize SUSYObjDef_xAOD...Aborting" << endl;
            abort();
        }
        
        std::cout << " INITIALIZED SUSYTOOLS with properties " << std::endl;
        for(auto& x:m_susyObj[susyObjId]->getPropertyMgr()->getProperties()){
            if(x.second->typeName()=="string"){
                string foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                cout << " Property << " << x.first << ": " << foo << endl;
            }
            else if(x.second->typeName()=="int"){
                int foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                cout << " Property << " << x.first << ": " << foo << endl;
            }
            else if(x.second->typeName()=="float"){
                float foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                cout << " Property << " << x.first << ": " << foo << endl;
            }
            else if(x.second->typeName()=="double"){
                double foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                cout << " Property << " << x.first << ": " << foo << endl;
            }
            else if(x.second->typeName()=="bool"){
                bool foo;
                m_susyObj[susyObjId]->getPropertyMgr()->getProperty(x.first, foo);
                string value = foo ? "True" : "False";
                cout << " Property << " << x.first << ": " << value << endl;
            }
        }
    //    CHECK( m_susyObj[i]->SUSYToolsInit() );
    }

    return *this;
}

//----------------------------------------------------------
XaodAnalysis& XaodAnalysis::initLocalTools()
{


    //initPileupTool(); // this is now done in SUSYTools (as of ST-00-06-15)
    initElectronTools();
    initPhotonTools();
    initMuonTools();
    initTauTools();
    initIsoTools();
    initStopPolReweightTool();

 //   initTrigger(); // now using SUSYTools for all our trigger needs

    return *this;
}
//----------------------------------------------------------
void XaodAnalysis::initStopPolReweightTool()
{
  if(m_polreweight == nullptr) { 
     m_polreweight = new StopPolarization::PolarizationReweight; 
     m_polreweight->setUnitMeV(); // set MeV
     m_polreweight->setMassW(80399.); 
     m_polreweight->setWidthW(2085.);
     m_polreweight->setMassZ(91187.6);
     m_polreweight->setWidthZ(2495.2);
     m_polreweight->setMassTop(172500.);
     m_polreweight->setWidthTop(1333.13);
     m_polreweight->setMassWThreshold(0.);
     m_polreweight->setMassZThreshold(0.);
     m_polreweight->setMassTopThreshold(54000.);
     std::string generatorName = "MadGraphPythia8";
     m_polreweight->setGeneratorName(generatorName);
     m_polreweight->setDecayPythia(true);
     m_polreweight->setPhaseSpaceOnly(true);
  }
}
//----------------------------------------------------------
void XaodAnalysis::initPileupTool()
{
    // This function is obsolete now that SUSYTools implements its own tool

    m_pileupReweightingTool = new CP::PileupReweightingTool("PileupReweightingTool");

    std::vector<std::string> prwFiles;
    std::vector<std::string> lumicalcFiles;

    if(!m_isMC15b && !m_isMC15c) {
        prwFiles.push_back("dev/SUSYTools/merged_prw.root"); //group 25ns
    }
    else if(m_isMC15b && !m_isMC15c) {
        prwFiles.push_back("dev/SUSYTools/merged_prw_mc15b.root");
    }
    else if(!m_isMC15b && m_isMC15c) {
        prwFiles.push_back("dev/SUSYTools/merged_prw_mc15c_July27_afterFix.root");
    }
    else {
        cout << "XaodAnalysis::initPileupTool    "
             << "Inconsistent MC options when setting up PRW config! Exiting." << endl;
        exit(1);
    }
    lumicalcFiles.push_back(m_data_dir+"SusyCommon/ilumicalc_histograms_None_297730-306451_OflLumi-13TeV-005.root"); // data16: updated Sept 15 2016 with GRL v82-pro20-12 (18.857 fb-1)
    lumicalcFiles.push_back(m_data_dir+"SusyCommon/ilumicalc_histograms_None_276262-284484.root"); // data15: updated with GRL v79


    cout << "XaodAnalysis::initPileupTool    Configuring SusyCommon PRW tool: " << endl;
    cout << "XaodAnalysis::initPileupTool     + prw config files..." << endl;
    for(int i = 0; i < (int)prwFiles.size(); i++) {
        cout << "XaodAnalysis::initPileupTool        > " << prwFiles[i] << endl;
    }
    cout << "XaodAnalysis::initPileupTool     + lumi calc files..." << endl;
    for(int i = 0; i < (int)lumicalcFiles.size(); i++) {
        cout << "XaodAnalysis::initPileupTool        > " << lumicalcFiles[i] << endl;
    }

    CHECK(m_pileupReweightingTool->setProperty("ConfigFiles", prwFiles));
    CHECK(m_pileupReweightingTool->setProperty("LumiCalcFiles", lumicalcFiles));
    CHECK(m_pileupReweightingTool->setProperty("DefaultChannel", 410000));
    CHECK(m_pileupReweightingTool->setProperty("DataScaleFactor",     1./ 1.16));
    CHECK(m_pileupReweightingTool->setProperty("DataScaleFactorUP",   1.));
    CHECK(m_pileupReweightingTool->setProperty("DataScaleFactorDOWN", 1./ 1.23));

    CHECK( m_pileupReweightingTool->initialize() );

}
//----------------------------------------------------------
void XaodAnalysis::initElectronTools()
{
    // Initialize the electron likelihood ID tools for each of the
    // working points

    // dantrim June 27 2016 -- the EG people decide to do this for us no, imagine that :)
    //std::string confDir                 = "ElectronPhotonSelectorTools/offline/";
    //std::string tight_LLH_conf          = "mc15_20160113/ElectronLikelihoodTightOfflineConfig2015.conf";
    //std::string medium_LLH_conf         = "mc15_20150712/ElectronLikelihoodMediumOfflineConfig2015.conf";
    //std::string loose_LLHBLayer_conf    = "mc15_20150712/ElectronLikelihoodLooseOfflineConfig2015_CutBL.conf"; 
    //std::string loose_LLH_conf          = "mc15_20150712/ElectronLikelihoodLooseOfflineConfig2015.conf";
    //std::string very_loose_LLH_conf     = "mc15_20150712/ElectronLikelihoodVeryLooseOfflineConfig2015.conf";


    std::string wp_veryloose = "VeryLooseLHElectron";
    std::string wp_loose = "LooseLHElectron";
    std::string wp_loose_blayer = "LooseBLLHElectron";
    std::string wp_medium = "MediumLHElectron";
    std::string wp_tight = "TightLHElectron";

    // VeryLooseLLH
    m_elecSelLikelihoodVeryLoose = new AsgElectronLikelihoodTool("AsgElectronLikelihoodToolVeryLoose");
    CHECK( m_elecSelLikelihoodVeryLoose->setProperty("WorkingPoint", wp_veryloose) );
    CHECK( m_elecSelLikelihoodVeryLoose->initialize() );

    // LooseLLH
    m_elecSelLikelihoodLoose = new AsgElectronLikelihoodTool("AsgElectronLikelihoodToolLoose");
    CHECK( m_elecSelLikelihoodLoose->setProperty("WorkingPoint", wp_loose) );
    CHECK( m_elecSelLikelihoodLoose->initialize() );

    // LooseAndBLayerLLH
    m_elecSelLikelihoodLooseBLayer = new AsgElectronLikelihoodTool("AsgElectronLikelihoodToolLooseBLayer");
    CHECK( m_elecSelLikelihoodLooseBLayer->setProperty("WorkingPoint", wp_loose_blayer) ); 
    CHECK( m_elecSelLikelihoodLooseBLayer->initialize() );
    
    // MediumLLH
    m_elecSelLikelihoodMedium = new AsgElectronLikelihoodTool("AsgElectronLikelihoodToolMedium");
    CHECK( m_elecSelLikelihoodMedium->setProperty("WorkingPoint", wp_medium) );
    CHECK( m_elecSelLikelihoodMedium->initialize() );

    // TightLLH
    m_elecSelLikelihoodTight = new AsgElectronLikelihoodTool("AsgElectronLikelihoodToolTight");
    CHECK( m_elecSelLikelihoodTight->setProperty("WorkingPoint", wp_tight) );
    CHECK( m_elecSelLikelihoodTight->initialize() );

}
//----------------------------------------------------------
void XaodAnalysis::initPhotonTools()
{
    // Initialize photon selection tools
    //std::string config_location = "ElectronPhotonSelectorTools/offline/mc15_20150712/";
    //std::string config_loose = "PhotonIsEMLooseSelectorCutDefs.conf";
    //std::string config_tight = "PhotonIsEMTightSelectorCutDefs.conf";

    //  > loose
    m_photonSelLoose = new AsgPhotonIsEMSelector("SusyCommonPhoLoose");
    CHECK( m_photonSelLoose->setProperty("WorkingPoint", "LoosePhoton") );
    CHECK( m_photonSelLoose->initialize() );

    //  > tight
    m_photonSelTight = new AsgPhotonIsEMSelector("SusyCommonPhoTight");
    CHECK( m_photonSelTight->setProperty("WorkingPoint", "TightPhoton") );
    CHECK( m_photonSelTight->initialize() );

}
//----------------------------------------------------------
void XaodAnalysis::initMuonTools()
{

    // Initialize muon selection tools

    m_muonSelectionToolVeryLoose = new CP::MuonSelectionTool("MuonSelectionTool_VeryLoose");
    CHECK( m_muonSelectionToolVeryLoose->setProperty( "MaxEta", 2.7 ) );
    CHECK( m_muonSelectionToolVeryLoose->setProperty( "MuQuality", int(xAOD::Muon::VeryLoose) ));// Warning: includes bad muons!
    CHECK( m_muonSelectionToolVeryLoose->initialize() );

    m_muonSelectionToolLoose = new CP::MuonSelectionTool("MuonSelectionTool_Loose");
    CHECK( m_muonSelectionToolLoose->setProperty( "MaxEta", 2.7 ) );
    CHECK( m_muonSelectionToolLoose->setProperty( "MuQuality", int(xAOD::Muon::Loose) ));
    CHECK( m_muonSelectionToolLoose->initialize() );
    
    m_muonSelectionToolMedium = new CP::MuonSelectionTool("MuonSelectionTool_Medium");
    CHECK( m_muonSelectionToolMedium->setProperty( "MaxEta", 2.7 ) );
    CHECK( m_muonSelectionToolMedium->setProperty( "MuQuality", int(xAOD::Muon::Medium) ));
    CHECK( m_muonSelectionToolMedium->initialize() );

    m_muonSelectionToolTight = new CP::MuonSelectionTool("MuonSelectionTool_Tight");
    CHECK( m_muonSelectionToolTight->setProperty( "MaxEta", 2.7 ) );
    CHECK( m_muonSelectionToolTight->setProperty( "MuQuality", int(xAOD::Muon::Tight) ));
    CHECK( m_muonSelectionToolTight->initialize() );
}
//----------------------------------------------------------
void XaodAnalysis::initTauTools()
{
    // Let's select some taus
    // > loose
    m_tauSelToolLoose = new TauAnalysisTools::TauSelectionTool("SusyCommonTauSelectionTool_Loose");
    CHECK( m_tauSelToolLoose->setProperty("ConfigPath", m_data_dir + "SUSYTools/tau_selection_loose.conf") );
    CHECK( m_tauSelToolLoose->initialize() );
    // > medium
    m_tauSelToolMedium = new TauAnalysisTools::TauSelectionTool("SusyCommonTauSelectionTool_Medium");
    CHECK( m_tauSelToolMedium->setProperty("ConfigPath", m_data_dir + "SUSYTools/tau_selection_medium.conf") );
    CHECK( m_tauSelToolMedium->initialize() );
    // > tight
    m_tauSelToolTight = new TauAnalysisTools::TauSelectionTool("SusyCommonTauSelectionTool_Tight");
    CHECK( m_tauSelToolTight->setProperty("ConfigPath", m_data_dir + "SUSYTools/tau_selection_tight.conf") );
    CHECK( m_tauSelToolTight->initialize() );

    // Truth Matching for all of the infos
    m_tauTruthMatchingTool = new TauAnalysisTools::TauTruthMatchingTool("SusyCommonTauTruthMatchingTool");
    CHECK( m_tauTruthMatchingTool->initialize() );
    m_tauTruthMatchingTool->msg().setLevel( MSG::INFO );

    //m_tauTruthMatchingTool = new TauAnalysisTools::TauTruthMatchingTool("TauTruthMatchingTool");  // memory leak check
    //m_tauTruthMatchingTool->msg().setLevel(m_dbg ? MSG::DEBUG : MSG::ERROR);                      // memory leak check
    //CHECK(m_tauTruthMatchingTool->initialize());                                                  // memory leak check
   // // dantrim Jul 2 :: following TauTruthTool setup found in SUSYTools
   // m_tauTruthMatchingTool = new TauAnalysisTools::TauTruthMatchingTool("TauTruthMatchingTool");
   // m_tauTruthMatchingTool->msg().setLevel(MSG::INFO);
   // CHECK( m_tauTruthMatchingTool->initialize() );   
 
    //AT: 05-07-15: consider adding ?
    //m_tauTruthMatchingTool->setProperty("SampleType", (int)TauAnalysisTools::SHERPA);
    //m_tauTruthMatchingTool->setProperty("SampleType", (int)TauAnalysisTools::PYTHIA);//default


    #warning TauEfficiencyTool should be initialized per TauSelTool WP
    #warning Tau SF should not be used currently
    m_TauEffEleTool = new TauAnalysisTools::TauEfficiencyCorrectionsTool("TauEfficiencyTool");
    m_TauEffEleTool->msg().setLevel(MSG::ERROR); // don't need warnings about histogram binning

    //CHECK(m_TauEffEleTool->setProperty("EfficiencyCorrectionType", (int) TauAnalysisTools::SFContJetID)); // this is mc12 rec, use default mc15
    CHECK(m_TauEffEleTool->initialize());
}
//----------------------------------------------------------
void XaodAnalysis::initIsoTools()
{
    // see: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/IsolationSelectionTool

    // Gradient Loose WP for leptons, FixedCutTight WP for photons
    m_isoToolGradientLooseTight = new CP::IsolationSelectionTool( "IsolationSelectionTool_GradientLoose" );
    CHECK( m_isoToolGradientLooseTight->setProperty( "ElectronWP",   "GradientLoose") );
    CHECK( m_isoToolGradientLooseTight->setProperty( "MuonWP",       "GradientLoose") );
    CHECK( m_isoToolGradientLooseTight->setProperty( "PhotonWP",     "FixedCutTight") );
    CHECK( m_isoToolGradientLooseTight->initialize() );
    
    // Gradient WP for leptons, FixedCutTightCaloOnly WP for photons
    m_isoToolGradientTightCalo = new CP::IsolationSelectionTool( "IsolationSelectionTool_Gradient" );
    CHECK( m_isoToolGradientTightCalo->setProperty( "ElectronWP",   "Gradient") );
    CHECK( m_isoToolGradientTightCalo->setProperty( "MuonWP",       "Gradient") );
    CHECK( m_isoToolGradientTightCalo->setProperty( "PhotonWP",     "FixedCutTightCaloOnly") );
    CHECK( m_isoToolGradientTightCalo->initialize() );

    // LooseTrackOnly WP, Loose WP for photons
    m_isoToolLooseTrackOnlyLoose = new CP::IsolationSelectionTool( "IsolationSelectionTool_LooseTrackOnly" );
    CHECK( m_isoToolLooseTrackOnlyLoose->setProperty( "ElectronWP",   "LooseTrackOnly") );
    CHECK( m_isoToolLooseTrackOnlyLoose->setProperty( "MuonWP",       "LooseTrackOnly") );
    CHECK( m_isoToolLooseTrackOnlyLoose->setProperty( "PhotonWP",     "FixedCutLoose") );
    CHECK( m_isoToolLooseTrackOnlyLoose->initialize() );
    
    // Loose WP for leptons, FixedCutTight WP for photons
    m_isoToolLoose = new CP::IsolationSelectionTool( "IsolationSelectionTool_Loose" );
    CHECK( m_isoToolLoose->setProperty( "ElectronWP",   "Loose") );
    CHECK( m_isoToolLoose->setProperty( "MuonWP",       "Loose") );
    CHECK( m_isoToolLoose->setProperty( "PhotonWP",     "FixedCutTight") );
    CHECK( m_isoToolLoose->initialize() );
    
    // Tight WP for leptons, FixedCutTight WP for photons
    m_isoToolTight = new CP::IsolationSelectionTool( "IsolationSelectionTool_Tight" );
    CHECK( m_isoToolTight->setProperty( "ElectronWP",   "Tight") );
    CHECK( m_isoToolTight->setProperty( "MuonWP",       "Tight") );
    CHECK( m_isoToolTight->setProperty( "PhotonWP",     "FixedCutTight") );
    CHECK( m_isoToolTight->initialize() );

}
//----------------------------------------------------------
void XaodAnalysis::initTrigger()
{
 //  Use SUSYTools implementation
 //   // dantrim trig
 //   m_configTool = new TrigConf::xAODConfigTool("xAODConfigTool");
 //   ToolHandle<TrigConf::ITrigConfigTool> configHandle(m_configTool);
 //   CHECK( configHandle->initialize() );
 //   
 //   m_trigTool = new Trig::TrigDecisionTool("TrigDecTool");
 //   m_trigTool->setProperty("ConfigTool", configHandle);
 //   m_trigTool->setProperty("TrigDecisionKey", "xTrigDecision");
 //   m_trigTool->setProperty("OutputLevel", MSG::ERROR).ignore(); // dantrim Mar 16 2015 -- tool outputs extraneous errors due to extraneous tool, ignore them
 //   CHECK( m_trigTool->initialize() );

 // Trigger matching is now done via SUSYTools

 //   // Tool for egamma matching
 //   m_trigEgammaMatchTool = new Trig::TrigEgammaMatchingTool("TrigEgammaMatchTool");
 //   CHECK( m_trigEgammaMatchTool->setProperty("TriggerTool", m_trigDec) );
 //   CHECK( m_trigEgammaMatchTool->initialize() );

}
/*--------------------------------------------------------------------------------*/
// Get the list of recommended systematics from CP
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::getSystematicList()
{
    if(m_dbg>=5) cout << "getSystematicList" << endl;
    //Get from SUSYTools the list of systematics and
    //what each systematics affects: weight/kin and object type
    systInfoList = m_susyObj[m_eleIDDefault]->getSystInfoList();
}

//----------------------------------------------------------
const xAOD::EventInfo* XaodAnalysis::retrieveEventInfo(xAOD::TEvent &e, bool dbg)
{
    const xAOD::EventInfo* evt = NULL;
    e.retrieve(evt, "EventInfo");
    if(dbg){
        if(evt) cout<<"XaodAnalysis::retrieveEventInfo: retrieved"<<endl;
        else    cout<<"XaodAnalysis::retrieveEventInfo: failed"<<endl;
    }
    return evt;
}
//----------------------------------------------------------
const xAOD::EventInfo* XaodAnalysis::xaodEventInfo()
{
    if(m_xaodEventInfo==NULL){
        m_xaodEventInfo = retrieveEventInfo(m_event, m_dbg);
    }
    return m_xaodEventInfo;
}
//----------------------------------------------------------
void XaodAnalysis::retrieveXaodTrackMet(ST::SystInfo sysInfo, SusyNtSys sys)
{
    if(m_dbg>=5) cout << "retrieveXaodTrackMet " << SusyNtSysNames[sys] << endl;
    
    m_trackMetContainer = new xAOD::MissingETContainer;
    m_trackMetAuxContainer = new xAOD::MissingETAuxContainer;
    m_trackMetContainer->setStore( m_trackMetAuxContainer );
    if(m_dbg>=5) cout << "Made trackMetContainer pointers " << m_trackMetContainer << " eleID " << m_eleIDDefault << endl;

    xAOD::ElectronContainer* electrons = xaodElectrons(sysInfo, sys);
    xAOD::MuonContainer*     muons     = xaodMuons(sysInfo, sys);
    xAOD::JetContainer*      jets      = xaodJets(sysInfo, sys);

    // GetTrackMET( met, jet, elec, muon)
    m_susyObj[m_eleIDDefault]->GetTrackMET(*m_trackMetContainer,
                                            jets,
                                            electrons,
                                            muons);

    if(m_dbg>=5) cout << "Rebult TrackMet with "
                      << " ele size " << electrons->size()
                      << " muo size " << muons->size()
                      << " jet size " << jets->size()
                      << endl;
}
    

//const xAOD::MissingETContainer* XaodAnalysis::xaodMET_Track()
//{
//    if (m_metTrackContainer == NULL){
//        if(m_dbg>=5) std::cout << "Retrieving track MET " << std::endl;
//        m_metTrackContainer = retrieveMET_Track(m_event, m_dbg);
//    }
//    return m_metTrackContainer;
//}
//----------------------------------------------------------
xAOD::MuonContainer* XaodAnalysis::xaodMuons(ST::SystInfo sysInfo, SusyNtSys sys)
{
    bool syst_affectsMuons     = ST::testAffectsObject(xAOD::Type::Muon, sysInfo.affectsType);
    if(sys!=NtSys::NOM && syst_affectsMuons){
        if(m_xaodMuons==NULL){
            CHECK( m_susyObj[m_eleIDDefault]->GetMuons(m_xaodMuons, m_xaodMuonsAux) );
        }
        if(m_dbg>=5) cout << "xaodMuo "<< m_xaodMuons->size() << endl;
        return m_xaodMuons;
    }
    else{
        if(m_xaodMuons_nom==NULL){
            CHECK( m_susyObj[m_eleIDDefault]->GetMuons(m_xaodMuons_nom, m_xaodMuonsAux_nom) );
            if(m_dbg>=5) cout << "xaodMuo_nom " << m_xaodMuons_nom->size() << endl;
        }
        return m_xaodMuons_nom;
    }
    return NULL;
}
//----------------------------------------------------------
xAOD::ElectronContainer* XaodAnalysis::xaodElectrons(ST::SystInfo sysInfo, SusyNtSys sys)
{
    bool syst_affectsElectrons = ST::testAffectsObject(xAOD::Type::Electron, sysInfo.affectsType);
    if(sys!=NtSys::NOM && syst_affectsElectrons){
        if(m_xaodElectrons==NULL){
            CHECK( m_susyObj[m_eleIDDefault]->GetElectrons(m_xaodElectrons, m_xaodElectronsAux) );
        }
        if(m_dbg>=5) cout << "xaodEle " << m_xaodElectrons->size() << endl;
        return m_xaodElectrons;
    }
    else{
        if(m_xaodElectrons_nom==NULL){
            CHECK( m_susyObj[m_eleIDDefault]->GetElectrons(m_xaodElectrons_nom, m_xaodElectronsAux_nom) );
            if(m_dbg>=5) cout << "xaodEle_nom " << m_xaodElectrons_nom->size() << endl;
        }
        return m_xaodElectrons_nom;
    }
    return NULL;
}
//----------------------------------------------------------
xAOD::TauJetContainer* XaodAnalysis::xaodTaus(ST::SystInfo sysInfo, SusyNtSys sys)
{
    bool syst_affectsTaus      = ST::testAffectsObject(xAOD::Type::Tau, sysInfo.affectsType);
    if(sys!=NtSys::NOM && syst_affectsTaus){
        if(m_xaodTaus==NULL){
            CHECK( m_susyObj[m_eleIDDefault]->GetTaus(m_xaodTaus, m_xaodTausAux) );
        }
        if(m_dbg>=5) cout << "xaodTaus " << m_xaodTaus->size() << endl;
        return m_xaodTaus;
    }
    else{
        if(m_xaodTaus_nom==NULL){
            CHECK( m_susyObj[m_eleIDDefault]->GetTaus(m_xaodTaus_nom, m_xaodTausAux_nom) );
            if(m_dbg>=5) cout << "xaodTaus_nom " << m_xaodTaus_nom->size() << endl;
        }
        return m_xaodTaus_nom;
    }
    return NULL;
}
//----------------------------------------------------------
xAOD::JetContainer* XaodAnalysis::xaodJets(ST::SystInfo sysInfo, SusyNtSys sys)
{
    bool syst_affectsJets      = ST::testAffectsObject(xAOD::Type::Jet, sysInfo.affectsType);
    if(sys!=NtSys::NOM && syst_affectsJets){
        if(m_xaodJets==NULL){
            CHECK( m_susyObj[m_eleIDDefault]->GetJets(m_xaodJets, m_xaodJetsAux) );
        }
        if(m_dbg>=5) cout << "xaodJets " << m_xaodJets->size() << endl;
        return m_xaodJets;
    }
    else{
        if(m_xaodJets_nom==NULL){
            CHECK( m_susyObj[m_eleIDDefault]->GetJets(m_xaodJets_nom, m_xaodJetsAux_nom) );
            if(m_dbg>=5) cout << "xaodJets_nom " << m_xaodJets_nom->size() << endl;
        }
        return m_xaodJets_nom;
    }
    return NULL;
}
//----------------------------------------------------------
xAOD::PhotonContainer* XaodAnalysis::xaodPhotons(ST::SystInfo sysInfo, SusyNtSys sys)
{
    bool syst_affectsPhotons   = ST::testAffectsObject(xAOD::Type::Photon, sysInfo.affectsType);
    if(sys!=NtSys::NOM && syst_affectsPhotons){
        if(m_xaodPhotons==NULL){
            CHECK( m_susyObj[m_eleIDDefault]->GetPhotons(m_xaodPhotons, m_xaodPhotonsAux) );
        }
        if(m_dbg>=5) cout << "xaodPho " << m_xaodPhotons->size()  << endl;
        return m_xaodPhotons;
    }
    else{
        if(m_xaodPhotons_nom==NULL){
            CHECK( m_susyObj[m_eleIDDefault]->GetPhotons(m_xaodPhotons_nom, m_xaodPhotonsAux_nom) );
            if(m_dbg>=5) cout << "xaodPho_nom " << m_xaodPhotons_nom->size() << endl;
        }

        return m_xaodPhotons_nom;
    }

    return NULL;
}
//----------------------------------------------------------
const xAOD::TruthEventContainer* XaodAnalysis::retrieveTruthEvent(xAOD::TEvent &e, bool dbg)
{
    const xAOD::TruthEventContainer* truth = NULL;
    e.retrieve(truth, "TruthEvents");
    if(dbg){
        if(truth) cout<<"XaodAnalysis::retrieveTruthEvents: retrieved "<<endl;
        else      cout<<"XaodAnalysis::retrieveTruthEvents: failed"<<endl;
    }
    return truth;
}
//----------------------------------------------------------
const xAOD::TruthEventContainer* XaodAnalysis::xaodTruthEvent()
{
    if(m_xaodTruthEvent==NULL && m_isMC){
        m_xaodTruthEvent = retrieveTruthEvent(m_event, m_dbg);
    }
    return m_xaodTruthEvent;
}
//----------------------------------------------------------
const xAOD::TruthParticleContainer* XaodAnalysis::retrieveTruthParticles(xAOD::TEvent &e, bool dbg)
{
    const xAOD::TruthParticleContainer* truthP = NULL;
    e.retrieve(truthP, "TruthParticles");
    if(dbg){
        if(truthP) cout<<"XaodAnalysis::retrieveTruthParticles: retrieved "<<truthP->size()<<endl;
        else       cout<<"XaodAnalysis::retrieveTruthParticles: failed"<<endl;
    }
    return truthP;
}
//----------------------------------------------------------
const xAOD::TruthParticleContainer* XaodAnalysis::xaodTruthParticles()
{
    if(m_xaodTruthParticles==NULL && m_isMC){
        m_xaodTruthParticles = retrieveTruthParticles(m_event, m_dbg);
    }
    return m_xaodTruthParticles;
}

//----------------------------------------------------------
//const xAOD::TruthParticleContainer* XaodAnalysis::xaodTruthTauParticles()  // test memory leak check
//{
//    return m_tauTruthMatchingTool->getTruthTauContainer();
//}
//----------------------------------------------------------
//const xAOD::TruthParticleAuxContainer* XaodAnalysis::xaodTruthTauParticlesAux()  // test memory leak check
//{
//    return m_tauTruthMatchingTool->getTruthTauAuxContainer();
//}

//----------------------------------------------------------
/// temporary patch, see SUSYToolsTester.cxx @ SUSYTools-00-05-00-14
bool muon_is_safe_for_met(const xAOD::Muon_v1 *mu)
{
    return (mu->muonType()==xAOD::Muon::Combined ||
            mu->muonType()==xAOD::Muon::SegmentTagged ||
            mu->muonType()==xAOD::Muon::MuonStandAlone);
}
/*--------------------------------------------------------------------------------*/
// Build MissingEt
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::retrieveXaodMet( ST::SystInfo sysInfo, SusyNtSys sys)
{
    if(m_dbg>=5) cout << "retrieveXaodMet " << SusyNtSysNames[sys] << endl;

    m_metContainer = new xAOD::MissingETContainer;
    m_metAuxContainer = new xAOD::MissingETAuxContainer;
    m_metContainer->setStore( m_metAuxContainer );
    if(m_dbg>=5) cout << "Made metContainer pointers " << m_metContainer << " eleID " << m_eleIDDefault << endl; 

    xAOD::ElectronContainer* electrons = xaodElectrons(sysInfo,sys);
    xAOD::MuonContainer*     muons     = xaodMuons(sysInfo,sys);
    xAOD::JetContainer*      jets      = xaodJets(sysInfo,sys);
    xAOD::PhotonContainer*   photons   = xaodPhotons(sysInfo,sys);

    // GetMET(met, jet, elec, muon, gamma, taujet, doTST = true, doJVT=true, invis = 0)
    m_susyObj[m_eleIDDefault]->GetMET(*m_metContainer,
                                      jets,
                                      electrons,
                                      muons,
                                      photons,
                                      0);
   
    if(m_dbg>=5) cout <<"Rebuilt MET with " 
                      << " ele size " << electrons->size()
                      << " jets size " << jets->size()
                      << " muons size " << muons->size()
                      << std::endl;

}
//----------------------------------------------------------
const xAOD::VertexContainer* XaodAnalysis::retrieveVertices(xAOD::TEvent &e, bool dbg)
{
    const xAOD::VertexContainer* vtx = NULL;
    e.retrieve(vtx, "PrimaryVertices");
    if(dbg){
        if(!vtx) cout<<"XaodAnalysis::retrieveVertices: failed"<<endl;
    }
    return vtx;
}
//----------------------------------------------------------
const xAOD::VertexContainer* XaodAnalysis::xaodVertices()
{
    if(m_xaodVertices==NULL){
        m_xaodVertices = retrieveVertices(m_event, m_dbg);
    }
    return m_xaodVertices;
}

//----------------------------------------------------------
void XaodAnalysis::selectBaselineObjects(SusyNtSys sys, ST::SystInfo sysInfo)
{
    if(m_dbg>=5) cout << "selectBaselineObjects with sys=" <<  SusyNtSysNames[sys] << endl;

    //////////////////////////////////
    // Grab the object containers
    //////////////////////////////////
    xAOD::ElectronContainer* electrons = xaodElectrons(sysInfo,sys);
    xAOD::MuonContainer*     muons     = xaodMuons(sysInfo,sys);
    xAOD::JetContainer*      jets      = xaodJets(sysInfo,sys);
    xAOD::TauJetContainer*   taus      = xaodTaus(sysInfo,sys);

    //////////////////////////////////
    // Electrons
    //////////////////////////////////
    int iEl = -1;
    for(const auto& el : *electrons) {
        iEl++;
        if(m_dbg>=5) cout<<"El "
                         <<" pt " << el->pt()
                         <<" eta " << el->eta()
                         <<" phi " << el->phi()
                         <<endl;
      //  if(!el->auxdata< bool >("baseline")) continue;
        //AT:12/16/14 TO UPDATE Base Obj should be after overlap removal
        if( (bool)el->auxdata< char >("baseline")==1 ) m_baseElectrons.push_back(iEl); //&&
        //    (bool)el->auxdata< char >("passOR")==1 ) m_baseElectrons.push_back(iEl); 

        if(m_dbg>=5) cout<<"\t El passing"
                         <<" baseline? "<< (bool)(el->auxdata< char >("baseline"))
                         <<" signal? "<<   (bool)(el->auxdata< char >("signal"))
                         <<endl;
        //AT 05-02-15: Minimum kinematic for electrons
        //dantrim May 5 2015 - thresholds a la ElectronEfficiencyCorrection
        const xAOD::CaloCluster* cluster = el->caloCluster();
        double et = cluster->e()/cosh(cluster->eta());
        if( et * MeV2GeV > 7 )
            m_preElectrons.push_back(iEl);
    }
    if(m_dbg) cout<<"preElectrons["<<m_preElectrons.size()<<"]"<<endl;

    //////////////////////////////////
    // Muons
    //////////////////////////////////
    int iMu = -1;
    for(const auto& mu : *muons){
        iMu++;
        m_preMuons.push_back(iMu);
        if(m_dbg>=5) cout<<"Mu passing"
                         <<" baseline? "<< bool(mu->auxdata< char >("baseline"))
                         <<" signal? "<<   bool(mu->auxdata< char >("signal"))
                         <<" pt " << mu->pt()
                         <<" eta " << mu->eta()
                         <<" phi " << mu->phi()
                         <<endl;
        if( (bool)mu->auxdata< char >("baseline")==1 ) m_baseMuons.push_back(iMu); //&& 
        //    (bool)mu->auxdata< char >("passOR")==1 ) m_baseMuons.push_back(iMu); 
        // if(signal) m_sigMuons.push_back(iMu);
    }
    if(m_dbg) cout<<"preMuons["<<m_preMuons.size()<<"]"<<endl;

    //////////////////////////////////
    // For updated jet selection (as of SUSY,2.3.15a) 
    // we need OR flags for Jets in order to check for b-tagging
    // and "bad" jets
    //
    // set OverlapRemoval flags
    //  --> "passOR"
    //  signature: (ele, muo, jets, useSigLep=false, useIsoLep=false, doBjetOR=false)
    //////////////////////////////////
    m_susyObj[m_eleIDDefault]->OverlapRemoval(electrons, muons, jets); 

    //////////////////////////////////
    // Jets
    //////////////////////////////////
    int iJet=-1;
    for(const auto& jet : *jets){
        iJet++;
        if(jet->pt()*MeV2GeV > 20) m_preJets.push_back(iJet);
        // defaults ST::00-06-15: IsBadJet(const xAOD::Jet&, jvtcut=0.64)
        // defaults ST::00-06-15: IsSignalJet(const xAOD::Jet&, ptcut=20000., etacut=2.8, jvtcut=0.64)
        //    !!--> IsBadJet must be called before IsSignalJet !!
        //    !!--> IsSignalJet must be called before IsBJet !!
        //    !!--> IsBadJet AND IsSignalJet must be called AFTER OverlapRemoval !!

        // defaults ST::00-06-13: mv2c20 > -0.4434 and dec_signal(jet)==True for IsBJet==True
        if(m_dbg>=5) cout<<"Jet passing"
                         <<" baseline? "<< bool(jet->auxdata< char >("baseline")==1)
                         <<" signal? "<<   bool(jet->auxdata< char >("signal")==1)
                         <<" pt " << jet->pt()
                         <<" eta " << jet->eta()
                         <<" phi " << jet->phi()
                         <<endl;
        if((bool)jet->auxdata< char >("baseline")==1 ) m_baseJets.push_back(iJet); //&&
         //  (bool)jet->auxdata< char >("passOR")==1 ) m_baseJets.push_back(iJet); 
    }
    if(m_dbg) cout<<"preJets["<<m_preJets.size()<<"]"<<endl;

    //////////////////////////////////
    // Taus
    //////////////////////////////////
    int iTau=-1;
    //xAOD::TauJetContainer* taus = xaodTaus(sysInfo,sys);
    for(const auto& tau : *taus){
        iTau++;
        if(m_dbg>=5) cout<<"Tau passing"
                         <<" baseline? "<< bool(tau->auxdata< char >("baseline")==1)
                         <<" signal? "<< bool(tau->auxdata< char >("signal")==1)
                         <<" pt " << tau->pt()
                         <<" eta " << tau->eta()
                         <<" phi " << tau->phi()
                         <<" q   " << tau->charge()
                         <<endl;
        //Container tau: tau->pt()>20*GeV && abs(tau->eta())<2.47 ???  //AT TO ADD
        //if(tau->pt() * MeV2GeV > 20 && fabs(tau->eta())<2.47) m_preTaus.push_back(iTau);        

        // MJF: Apply very loose pre-selection to container taus
        int nTracks = tau->nTracks();
        if (std::abs(tau->charge()) == 1 && (nTracks==1 || nTracks==3 || nTracks==5)) {
            m_contTaus.push_back(iTau);
        }
        if((bool)tau->auxdata< char >("baseline")==1){
            m_preTaus.push_back(iTau);
            m_baseTaus.push_back(iTau);
        }
    }
    //if(m_dbg) cout<<"m_preTaus["<<m_preTaus.size()<<"]"<<endl;

    //////////////////////////////////
    // If Nom systematics keep track
    // of the pre_object indices
    //////////////////////////////////
    if(sys==NtSys::NOM){
        m_preElectrons_nom = m_preElectrons;
        m_preMuons_nom     = m_preMuons;
        m_preJets_nom      = m_preJets;
        m_contTaus_nom     = m_contTaus;
        m_preTaus_nom      = m_preTaus;
        m_preLeptons_nom   = m_preLeptons;
    }

/**/
    // Container object selection
    //-DG-if(m_selectTaus) m_contTaus = get_taus_baseline(xaodTaus(), m_susyObj, 20.*GeV, 2.47,
    //-DG-                                                SUSYTau::TauNone, SUSYTau::TauNone, SUSYTau::TauNone,
    //-DG-                                                susySys, true);

    // Preselection
    //-DG-m_preElectrons = get_electrons_baseline(xaodElectrons(), &m_event.el_MET_Egamma10NoTau,
    //-DG-                                        !m_isMC, m_event.eventinfo.RunNumber(), m_susyObj,
    //-DG-                                        7.*GeV, 2.47, susySys);
    //-DG-m_preMuons = get_muons_baseline(xaodMuons(), !m_isMC, m_susyObj,
    //-DG-                                6.*GeV, 2.5, susySys);
    // Removing eta cut for baseline jets. This is for the bad jet veto.
    //-DG-m_preJets = get_jet_baseline(jets, &m_event.vxp, &m_event.eventinfo, &m_event.Eventshape, !m_isMC, m_susyObj,
    //-DG-                             20.*GeV, std::numeric_limits<float>::max(), susySys, false, goodJets);

    // Selection for met muons
    // Diff with preMuons is pt selection
    //-DG-m_metMuons = get_muons_baseline(xaodMuons(), !m_isMC, m_susyObj,
    //-DG-                                10.*GeV, 2.5, susySys);

    // Preselect taus
    //-DG-if(m_selectTaus) m_preTaus = get_taus_baseline(xaodTaus(), m_susyObj, 20.*GeV, 2.47,
    //-DG-                                               SUSYTau::TauLoose, SUSYTau::TauLoose, SUSYTau::TauLoose,
    //-DG-                                               susySys, true);
    //-DG-performOverlapRemoval();

    // combine leptons
    //-DG-m_preLeptons    = buildLeptonInfos(xaodElectrons(), m_preElectrons, xaodMuons(), m_preMuons, m_susyObj);
    //-DG-m_baseLeptons   = buildLeptonInfos(xaodElectrons(), m_baseElectrons, xaodMuons(), m_baseMuons, m_susyObj);
}

/*--------------------------------------------------------------------------------*/
// perform overlap
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::performOverlapRemoval()
{
//-DG-  // e-e overlap removal
//-DG-  m_baseElectrons = overlap_removal(m_susyObj, xaodElectrons(), m_preElectrons, xaodElectrons(), m_preElectrons,
//-DG-                                    0.05, true, true);
//-DG-  // jet-e overlap removal
//-DG-  m_baseJets      = overlap_removal(m_susyObj, jets, m_preJets, xaodElectrons(), m_baseElectrons,
//-DG-                                    0.2, false, false);
//-DG-
//-DG-  if(m_selectTaus) {
//-DG-    // tau-e overlap removal
//-DG-    m_baseTaus    = overlap_removal(m_susyObj, xaodTaus(), m_preTaus, xaodElectrons(), m_baseElectrons, 0.2, false, false);
//-DG-    // tau-mu overlap removal
//-DG-    m_baseTaus    = overlap_removal(m_susyObj, xaodTaus(), m_baseTaus, xaodMuons(), m_preMuons, 0.2, false, false);
//-DG-  }
//-DG-
//-DG-  // e-jet overlap removal
//-DG-  m_baseElectrons = overlap_removal(m_susyObj, xaodElectrons(), m_baseElectrons, jets, m_baseJets,
//-DG-                                    0.4, false, false);
//-DG-
//-DG-  // m-jet overlap removal
//-DG-  m_baseMuons     = overlap_removal(m_susyObj, xaodMuons(), m_preMuons, jets, m_baseJets, 0.4, false, false);
//-DG-
//-DG-  // e-m overlap removal
//-DG-  vector<int> copyElectrons = m_baseElectrons;
//-DG-  m_baseElectrons = overlap_removal(m_susyObj, xaodElectrons(), m_baseElectrons, xaodMuons(), m_baseMuons,
//-DG-                                    0.01, false, false);
//-DG-  m_baseMuons     = overlap_removal(m_susyObj, xaodMuons(), m_baseMuons, xaodElectrons(), copyElectrons, 0.01, false, false);
//-DG-
//-DG-  // m-m overlap removal
//-DG-  m_baseMuons     = overlap_removal(m_susyObj, xaodMuons(), m_baseMuons, xaodMuons(), m_baseMuons, 0.05, true, false);
//-DG-
//-DG-  // jet-tau overlap removal
//-DG-  m_baseJets      = overlap_removal(m_susyObj, jets, m_baseJets, xaodTaus(), m_baseTaus, 0.2, false, false);
//-DG-
//-DG-  // remove SFOS lepton pairs with Mll < 12 GeV
//-DG-  m_baseElectrons = RemoveSFOSPair(m_susyObj, xaodElectrons(), m_baseElectrons, 12.*GeV);
//-DG-  m_baseMuons     = RemoveSFOSPair(m_susyObj, xaodMuons(), m_baseMuons,     12.*GeV);
//-DG-  //m_baseTaus      = RemoveSFOSPair(m_susyObj, xaodTaus(), m_baseTaus,      12.*GeV);
}

/*--------------------------------------------------------------------------------*/
// Signal object selection - do baseline selection first!
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::selectSignalObjects(SusyNtSys sys, ST::SystInfo sysInfo)
{
    if(m_dbg>=5) cout << "selectSignalObjects with sys=" <<  SusyNtSysNames[sys] << endl;

    //////////////////////////////////
    // Grab the object containers
    //////////////////////////////////
    xAOD::ElectronContainer* electrons = xaodElectrons(sysInfo, sys);
    xAOD::MuonContainer*     muons     = xaodMuons(sysInfo, sys);
    xAOD::JetContainer*      jets      = xaodJets(sysInfo, sys);
    xAOD::TauJetContainer*   taus      = xaodTaus(sysInfo, sys);
    xAOD::PhotonContainer*   photons   = xaodPhotons(sysInfo, sys);

    //////////////////////////////////
    // Electrons
    //////////////////////////////////
    int iEl = 0;
    //xAOD::ElectronContainer* electrons = xaodElectrons(sysInfo,sys);
    for(const auto& el : *electrons) {
        if( (bool)el->auxdata< char >("signal")==1 &&
            (bool)el->auxdata< char >("passOR")==1 )   m_sigElectrons.push_back(iEl);
            iEl++;
    }
    if(m_dbg) cout<<"m_sigElectrons["<<m_sigElectrons.size()<<"]"<<endl;

    //////////////////////////////////
    // Muons
    //////////////////////////////////
    int iMu = 0;
    //xAOD::MuonContainer* muons = xaodMuons(sysInfo,sys);
    for(const auto& mu : *muons){
        if( (bool)mu->auxdata< char >("signal")==1 &&
            (bool)mu->auxdata< char >("passOR")==1 )  m_sigMuons.push_back(iMu);
            iMu++;
    }
    if(m_dbg) cout<<"m_sigMuons["<<m_sigMuons.size()<<"]"<<endl;

    //////////////////////////////////
    // Jets
    //////////////////////////////////
    int iJet=0;
    //xAOD::JetContainer* jets = xaodJets(sysInfo,sys);
    for(const auto& jet : *jets){
        if(jet->pt()*MeV2GeV >20.0 &&
            jet->auxdata< char >("passOR")==1 &&
           !((bool)jet->auxdata< char >("bad")==1) //AT: Added 12/13/14
           // DG tmp-2014-11-02 (!jet->isAvailable("bad") || !jet->auxdata< bool >("bad"))
            )
            m_sigJets.push_back(iJet);
        iJet++;
    }
    if(m_dbg) cout<<"m_sigJets["<<m_sigJets.size()<<"]"<<endl;

    //////////////////////////////////
    // Taus
    //////////////////////////////////
    int iTau=0;
    //xAOD::TauJetContainer* taus = xaodTaus(sysInfo,sys);
    for(const auto& tau : *taus){
        if(tau->pt() * MeV2GeV >20.0 &&
           (bool)tau->auxdata< char >("signal")==1)
            // tau->auxdata< int >("passOR") && // tau not involved in OR?
            m_sigTaus.push_back(iTau);
        }
        iTau++;
    if(m_dbg) cout<<"m_sigTaus["<<m_sigTaus.size()<<"]"<<endl;

    //////////////////////////////////
    // Photons
    //////////////////////////////////
    int iPh=0;
    //xAOD::PhotonContainer* photons = xaodPhotons(sysInfo,sys);
    if(photons) {
        for(const auto& ph : *photons){
            if( (bool)ph->auxdata< char >("passCleaning") &&
                (bool)ph->auxdata< char >("passAmbiguity") )
                m_sigPhotons.push_back(iPh);
            iPh++;
        }
    }
    if(m_dbg && photons) cout<<"m_sigPhotons["<<m_sigPhotons.size()<<"]"<<endl;

}

/*--------------------------------------------------------------------------------*/
// Signal photons
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::selectSignalPhotons()
{
//-DG-  if(m_dbg>=5) cout << "selectSignalPhotons" << endl;
//-DG-
//-DG-  int phoQual = 2;      // Quality::Tight
//-DG-  uint isoType = 1;     // Corresponds to PTED corrected isolation
//-DG-  float etcone40CorrCut = 3*GeV;
//-DG-
//-DG-  vector<int> base_photons = get_photons_baseline(&m_event.ph, m_susyObj,
//-DG-                                                  20.*GeV, 2.47, SystErr::NONE, phoQual);
//-DG-
//-DG-  // Latest and Greatest
//-DG-  int nPV = getNumGoodVtx();
//-DG-  m_sigPhotons = get_photons_signal(&m_event.ph, base_photons, m_susyObj, nPV,
//-DG-                                    20.*GeV, etcone40CorrCut, isoType);
}
/*--------------------------------------------------------------------------------*/
// Truth object selection
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::selectTruthObjects()
{
//-DG-  if(m_dbg>=5) cout << "selectTruthObjects" << endl;
//-DG-
//-DG-  // ==>> First the truth particles
//-DG-  // Done under SusyNtMaker::fillTruthParticleVars
//-DG-
//-DG-  // ==>> Second the truth jets
//-DG-  for(int index=0; index < m_event.AntiKt4Truth.n(); index++) {
//-DG-      // const xAOD::JetD3PDObjectElement &trueJet = m_event.AntiKt4Truth[index];
//-DG-      // if( trueJet.pt()/GeV > 15. && fabs(trueJet.eta()) < 4.5) m_truJets.push_back(index);
//-DG-#warning truth not implemented
//-DG-  }
//-DG-
//-DG-  // ==>> Third and last the truth met
//-DG-  m_truMet.SetPxPyPzE(m_event.MET_Truth.NonInt_etx(), m_event.MET_Truth.NonInt_ety(), 0, m_event.MET_Truth.NonInt_sumet());
}
//----------------------------------------------------------
void XaodAnalysis::clearOutputObjects(bool deleteNominal)
{
    m_preElectrons.clear();
    m_preMuons.clear();
    m_preJets.clear();
    m_contTaus.clear();
    m_preTaus.clear();
    m_preLeptons.clear();
    m_baseElectrons.clear();
    m_baseMuons.clear();
    m_baseTaus.clear();
    m_baseLeptons.clear();
    m_baseJets.clear();
    m_sigElectrons.clear();
    m_sigMuons.clear();
    m_sigLeptons.clear();
    m_sigJets.clear();
    m_sigTaus.clear();
    m_cutFlags = 0;

    m_sigPhotons.clear();
    m_truParticles.clear();
    m_truJets.clear();

    m_metMuons.clear();

    if(deleteNominal){
        m_preElectrons_nom.clear();
        m_preMuons_nom.clear();
        m_preJets_nom.clear();
        m_contTaus_nom.clear();
        m_preTaus_nom.clear();
        m_preLeptons_nom.clear();
    }


}

/*--------------------------------------------------------------------------------*/
// Count number of good vertices
/*--------------------------------------------------------------------------------*/
uint XaodAnalysis::getNumGoodVtx()
{
/*
    xAOD::VertexContainer::const_iterator pv_itr = m_xaodVertices->begin();
    //AT-2014-10-31: Run2 harmonisation - we don't need to cut on nTrack for run@
    //https://cds.cern.ch/record/1700874/files/ATL-COM-PHYS-2014-451.pdf

    uint nVtx = 0;
    for(auto it=m_xaodVertices->begin(), end=m_xaodVertices->end(); it!=end; ++it){
        const xAOD::Vertex &vtx = **it;
        if(vtx.nTrackParticles() >=5 ) nVtx++;
    }
    return nVtx;
*/
    //AT:2014-10-31: To be change to this for run2 : 2
    return  xaodVertices()->size();

}

/*--------------------------------------------------------------------------------*/
// Count number of good vertices
/*--------------------------------------------------------------------------------*/
const xAOD::Vertex* XaodAnalysis::getPV()
{
    xAOD::Vertex* vtx = NULL;
    const xAOD::VertexContainer* vertices = xaodVertices();
    for(auto it=vertices->begin(), end=vertices->end(); it!=end; ++it){
        if((**it).vertexType()==xAOD::VxType::PriVtx)  vtx = *it;
    }
    return vtx;
}

/*--------------------------------------------------------------------------------*/
// Match reco jet to a truth jet
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::matchTruthJet(int iJet)
{
    // // Loop over truth jets looking for a match
    // const TLorentzVector &jetLV = m_susyObj[m_eleIDDefault]->GetJetTLV(iJet);
    // for(int i=0; i<m_event.AntiKt4Truth.n(); i++){
    //   // const xAOD::JetD3PDObjectElement &trueJet = m_event.AntiKt4Truth[i];
    //   // TLorentzVector trueJetLV;
    //   // trueJetLV.SetPtEtaPhiE(trueJet.pt(), trueJet.eta(), trueJet.phi(), trueJet.E());
    //   // if(jetLV.DeltaR(trueJetLV) < 0.3) return true;
    //   }
    return false;
}

/*--------------------------------------------------------------------------------*/
// Return electron type
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::eleIsOfType(const xAOD::Electron &in, ElectronId id)
{
    if     (id==ElectronId::VeryLooseLLH  && m_elecSelLikelihoodVeryLoose->accept(in))  return true;
    else if(id==ElectronId::LooseLLH  && m_elecSelLikelihoodLoose->accept(in))  return true;
    else if(id==ElectronId::LooseLLHBLayer && m_elecSelLikelihoodLooseBLayer->accept(in)) return true;
    else if(id==ElectronId::MediumLLH && m_elecSelLikelihoodMedium->accept(in)) return true;
    else if(id==ElectronId::TightLLH  && m_elecSelLikelihoodTight->accept(in))  return true;

    return false;
}
/*--------------------------------------------------------------------------------*/
// Return muon type
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::muIsOfType(const xAOD::Muon &in, MuonId id)
{
    if     (id==MuonId::VeryLoose && m_muonSelectionToolVeryLoose->accept(in))  return true;
    else if(id==MuonId::Loose     && m_muonSelectionToolLoose    ->accept(in))  return true;
    else if(id==MuonId::Medium    && m_muonSelectionToolMedium   ->accept(in))  return true;
    else if(id==MuonId::Tight     && m_muonSelectionToolTight    ->accept(in))  return true;
    return false;
} 
/*--------------------------------------------------------------------------------*/
// Get triggers
/*--------------------------------------------------------------------------------*/
std::vector<std::string> XaodAnalysis::xaodTriggers()
{
    if(m_triggerNames.size()==0){
        m_triggerNames = getTrigNames(m_triggerSet);
        return m_triggerNames;
    }
    else { return m_triggerNames; }
}
/*--------------------------------------------------------------------------------*/
// Event trigger flags
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::fillEventTriggers()
{
    if(m_dbg>=5) cout << "fillEventTriggers" << endl;
    m_evtTrigBits.ResetAllBits();
    std::vector<std::string> trigs = XaodAnalysis::xaodTriggers();
    for (unsigned int iTrig = 0; iTrig < trigs.size(); iTrig++) {
        if(m_susyObj[m_eleIDDefault]->IsTrigPassed(trigs[iTrig])) m_evtTrigBits.SetBitNumber(iTrig, true);
        //if(m_trigTool->isPassed(trigs[iTrig]))  m_evtTrigBits.SetBitNumber(iTrig, true);
    }
}
/*--------------------------------------------------------------------------------*/
// Muon trigger matching
/*--------------------------------------------------------------------------------*/
TBits XaodAnalysis::matchMuonTriggers(const xAOD::Muon &in)
{
    if(m_dbg>=10) cout << "XaodAnalysis::matchMuonTriggers" << endl;
    TBits muoTrigBits(m_nTriggerBits);
    muoTrigBits.ResetAllBits();
    std::vector<std::string> trigs = XaodAnalysis::xaodTriggers();
    for(unsigned int iTrig=0; iTrig<trigs.size(); iTrig++){

        std::size_t elfound = trigs[iTrig].find("el");
        std::size_t mufound = trigs[iTrig].find("mu");
        std::size_t mufound2 = trigs[iTrig].find("MU");
        std::size_t metfound = trigs[iTrig].find("xe");

        bool mu_chain = ((mufound != std::string::npos) || (mufound2 != std::string::npos));
        bool el_chain = (elfound != std::string::npos);
        bool met_chain = (metfound != std::string::npos); 
        bool dilepton = mu_chain && el_chain;

        if( (mu_chain || dilepton) && !met_chain)
            if(m_susyObj[m_eleIDDefault]->IsTrigMatched(&in, trigs[iTrig])) muoTrigBits.SetBitNumber(iTrig, true);
            //if(m_susyObj[m_eleIDDefault]->IsTrigMatched(&in, trigs[iTrig])) { cout << "     > muon match to : " << trigs[iTrig] << endl;  muoTrigBits.SetBitNumber(iTrig, true); }
    }
    return muoTrigBits;
}
/*--------------------------------------------------------------------------------*/
// Electron trigger matching
/*--------------------------------------------------------------------------------*/
TBits XaodAnalysis::matchElectronTriggers(const xAOD::Electron &in)
{
    // DA: maybe split up trigger list in terms of electron/muon/etc triggers? but then bit numbers
    // may be out of sync w.r.t. the stored histogram... in any case, the non-passed triggers are always
    // false and on SusyNtuple side the user will provide specific string for the ele/muo trigger that
    // is checked against the trig histo
    if(m_dbg>=10) cout << "XaodAnalysis::matchElectronTriggers" << endl;
    TBits eleTrigBits(m_nTriggerBits);
    eleTrigBits.ResetAllBits();
    std::vector<std::string> trigs = XaodAnalysis::xaodTriggers();
    for(unsigned int iTrig = 0; iTrig < trigs.size(); iTrig++) {

        std::size_t elfound = trigs[iTrig].find("HLT_e");
        std::size_t elfound2 = trigs[iTrig].find("HLT_2e");
        std::size_t mufound = trigs[iTrig].find("mu");
        std::size_t mufound2 = trigs[iTrig].find("MU");
        std::size_t metfound = trigs[iTrig].find("xe");

        bool mu_chain = ((mufound != std::string::npos) || (mufound2 != std::string::npos));
        bool el_chain = ((elfound != std::string::npos) || (elfound2 != std::string::npos));
        bool met_chain = (metfound != std::string::npos); 
        bool dilepton = mu_chain && el_chain;

        bool ismatch = false;
        if( (el_chain || dilepton) && !met_chain)
            ismatch = m_susyObj[m_eleIDDefault]->IsTrigMatched(&in, trigs[iTrig]);
        if(ismatch) {
            //cout << "     > ele match to : " << trigs[iTrig] << endl;
            eleTrigBits.SetBitNumber(iTrig,true);
        }
    }
    return eleTrigBits;
}
/*--------------------------------------------------------------------------------*/
// Electron trigger matching
/*--------------------------------------------------------------------------------*/
//void XaodAnalysis::matchElectronTriggers()
//{
//    if(m_dbg>=5) cout << "matchElectronTriggers" << endl;
//    for ( uint i=0; i < m_preElectrons.size(); i++ ) {
//        int iEl = m_preElectrons[i];
//        xAOD::ElectronContainer* nom_electrons = xaodElectrons(systInfoList[0]);
//        const xAOD::Electron* electron = nom_electrons->at(iEl);
//        long long flags = 0;
//        TLorentzVector* ele_tlv;
//        ele_tlv->SetPtEtaPhiM(electron->pt(), electron->eta(), electron->phi(), 0);
//        if(matchElectronTrigger(ele_tlv, "EF_e7_medium1"))     flags |= TRIG_e7_medium1;
//        if(matchElectronTrigger(ele_tlv, "EF_e12Tvh_loose1"))  flags |= TRIG_e12Tvh_loose1;
//        if(matchElectronTrigger(ele_tlv, "EF_e12Tvh_medium1")) flags |= TRIG_e12Tvh_medium1;
//
//        m_eleTrigFlags[iEl] = flags;
//    }
//}

/*--------------------------------------------------------------------------------*/
//bool XaodAnalysis::matchElectronTrigger(const TLorentzVector* lv, string chain)
//{
//    // TODO : figure out why they have software releases that do not have correctly updated aligned dependencies....
//
//    // get chain group representing the requesting trigger 
//    bool ele_ismatch = false;
//    auto cg = m_trigTool->getChainGroup(chain);
//    auto fc = cg->features();
//    auto eleFeatureContainers = fc.containerFeature<xAOD::TrigElectronContainer>();
//    for(auto &econt : eleFeatureContainers) {
//        for ( auto trige : *econt.cptr() ) {
//            TLorentzVector* trig_tlv;
//            trig_tlv->SetPtEtaPhiM( trige->pt(), trige->eta(), trige->phi(), 0);
//            float dR = lv->DeltaR(*trig_tlv); 
//            if ( dR < 0.15 ) ele_ismatch = true;
//        } // trige
//    } // econt
//    return ele_ismatch;
//}
        
        




//-DG--  for(uint i=0; i<m_preElectrons.size(); i++){
//-DG--    int iEl = m_preElectrons[i];
//-DG--    const TLorentzVector &lv = m_susyObj[m_eleIDDefault]->GetElecTLV(iEl);
//-DG--    // trigger flags
//-DG--    long long flags = 0;
//-DG--    // 2012 triggers only
//-DG--    if(matchElectronTrigger(lv, m_event.trig_EF_el.EF_e7T_medium1()))               { flags |= TRIG_e7_medium1; }
//-DG--    if(matchElectronTrigger(lv, m_event.trig_EF_el.EF_e12Tvh_loose1()))             { flags |= TRIG_e12Tvh_loose1; }
//-DG--    if(matchElectronTrigger(lv, m_event.trig_EF_el.EF_e12Tvh_medium1()))            { flags |= TRIG_e12Tvh_medium1; }
//-DG--    if(matchElectronTrigger(lv, m_event.trig_EF_el.EF_e24vh_medium1()))             { flags |= TRIG_e24vh_medium1; }
//-DG--    if(matchElectronTrigger(lv, m_event.trig_EF_el.EF_e24vhi_medium1()))            { flags |= TRIG_e24vhi_medium1; }
//-DG--    if(matchElectronTrigger(lv, m_event.trig_EF_el.EF_2e12Tvh_loose1()))            { flags |= TRIG_2e12Tvh_loose1; }
//-DG--    if(matchElectronTrigger(lv, m_event.trig_EF_el.EF_e24vh_medium1_e7_medium1()))  { flags |= TRIG_e24vh_medium1_e7_medium1; }
//-DG--    if(matchElectronTrigger(lv, m_event.trig_EF_el.EF_e12Tvh_medium1_mu8()))        { flags |= TRIG_e12Tvh_medium1_mu8; }
//-DG--    if(matchElectronTrigger(lv, m_event.trig_EF_el.EF_e18vh_medium1()))             { flags |= TRIG_e18vh_medium1; }
//-DG--    if(matchElectronTrigger(lv, m_event.trig_EF_el.EF_e18vh_medium1_2e7T_medium1())){ flags |= TRIG_e18vh_medium1_2e7T_medium1; }
//-DG--    if(matchElectronTrigger(lv, m_event.trig_EF_el.EF_2e7T_medium1_mu6()))          { flags |= TRIG_2e7T_medium1_mu6; }
//-DG--    if(matchElectronTrigger(lv, m_event.trig_EF_el.EF_e7T_medium1_2mu6()))          { flags |= TRIG_e7T_medium1_2mu6; }
//-DG--    if(matchElectronTrigger(lv, m_event.trig_EF_el.EF_e24vh_medium1_EFxe35_tclcw())){ flags |= TRIG_e24vh_medium1_EFxe35_tclcw; }
//-DG--    m_eleTrigFlags[iEl] = flags;
//-DG--  }
//}
/*--------------------------------------------------------------------------------*/
//bool XaodAnalysis::matchElectronTrigger(const TLorentzVector &lv, vector<int>* trigBools)
//{
//    // matched trigger index - not used
//    //static int indexEF = -1;
//    // Use function defined in egammaAnalysisUtils/egammaTriggerMatching.h
//    // return PassedTriggerEF(lv.Eta(), lv.Phi(), trigBools, indexEF, m_event.trig_EF_el.n(),
//    //                        m_event.trig_EF_el.eta(), m_event.trig_EF_el.phi());
//    return false;
//}

/*--------------------------------------------------------------------------------*/
// Muon trigger matching
/*--------------------------------------------------------------------------------*/
//void XaodAnalysis::matchMuonTriggers()
//{
//-DG--  if(m_dbg>=5) cout << "matchMuonTriggers" << endl;
//-DG--  for(uint i=0; i<m_preMuons.size(); i++){
//-DG--    int iMu = m_preMuons[i];
//-DG--    const TLorentzVector &lv = m_susyObj[m_eleIDDefault]->GetMuonTLV(iMu);
//-DG--    long long flags = 0;
//-DG--    // 2012 triggers only
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu8()) )         { flags |= TRIG_mu8; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu13()) )        { flags |= TRIG_mu13; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu18_tight()) )  { flags |= TRIG_mu18_tight; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu24i_tight()) ) { flags |= TRIG_mu24i_tight; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_2mu13()) )       { flags |= TRIG_2mu13; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu18_tight_mu8_EFFS()) )   { flags |= TRIG_mu18_tight_mu8_EFFS; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu18_tight_e7_medium1()) ) { flags |= TRIG_mu18_tight_e7_medium1; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu15()) )                  { flags |= TRIG_mu15; }
//-DG--    if(!m_isMC && m_event.eventinfo.RunNumber()>=206248 &&
//-DG--       matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_2mu8_EFxe40wMu_tclcw()))   { flags |= TRIG_2mu8_EFxe40wMu_tclcw; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu6()) )                  { flags |= TRIG_mu6; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_2mu6()) )                 { flags |= TRIG_2mu6; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu18_tight_2mu4_EFFS()) ) { flags |= TRIG_mu18_tight_2mu4_EFFS; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu4T()) )                 { flags |= TRIG_mu4T; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu24()) )                 { flags |= TRIG_mu24; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu4T_j65_a4tchad_xe70_tclcw_veryloose()) ) { flags |= TRIG_mu4T_j65_a4tchad_xe70_tclcw_veryloose; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_2mu4T_xe60_tclcw()) ) { flags |= TRIG_2mu4T_xe60_tclcw; }
//-DG--    if(m_event.trig_EF_trigmuonef.EF_2mu8_EFxe40_tclcw.IsAvailable() &&
//-DG--       matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_2mu8_EFxe40_tclcw()) ) { flags |= TRIG_2mu8_EFxe40_tclcw; }
//-DG--    if( matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu24_j65_a4tchad_EFxe40_tclcw()) ) { flags |= TRIG_mu24_j65_a4tchad_EFxe40_tclcw; }
//-DG--    if(m_event.trig_EF_trigmuonef.EF_mu24_j65_a4tchad_EFxe40wMu_tclcw.IsAvailable() &&
//-DG--       matchMuonTrigger(lv, m_event.trig_EF_trigmuonef.EF_mu24_j65_a4tchad_EFxe40wMu_tclcw()) ) { flags |= TRIG_mu24_j65_a4tchad_EFxe40wMu_tclcw; }
//-DG--    m_muoTrigFlags[iMu] = flags;
//-DG--  }
//}
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::matchMuonTrigger(const TLorentzVector &lv, vector<int>* passTrig)
{
//-DG--  // loop over muon trigger features
//-DG--  for(int iTrig=0; iTrig < m_event.trig_EF_trigmuonef.n(); iTrig++){
//-DG--
//-DG--    // Check to see if this feature passed chain we want
//-DG--    if(passTrig->at(iTrig)){
//-DG--
//-DG--      // Loop over muon EF tracks
//-DG--      TLorentzVector lvTrig;
//-DG--      for(int iTrk=0; iTrk < m_event.trig_EF_trigmuonef.track_n()->at(iTrig); iTrk++){
//-DG--
//-DG--        lvTrig.SetPtEtaPhiM( m_event.trig_EF_trigmuonef.track_CB_pt()->at(iTrig).at(iTrk),
//-DG--                             m_event.trig_EF_trigmuonef.track_CB_eta()->at(iTrig).at(iTrk),
//-DG--                             m_event.trig_EF_trigmuonef.track_CB_phi()->at(iTrig).at(iTrk),
//-DG--                             0 );       // only eta and phi used to compute dR anyway
//-DG--        // Require combined offline track...?
//-DG--        if(!m_event.trig_EF_trigmuonef.track_CB_hasCB()->at(iTrig).at(iTrk)) continue;
//-DG--        float dR = lv.DeltaR(lvTrig);
//-DG--        if(dR < 0.15){
//-DG--          return true;
//-DG--        }
//-DG--
//-DG--      } // loop over EF tracks
//-DG--    } // trigger object passes chain?
//-DG--  } // loop over trigger objects
//-DG--
//-DG--  // matching failed
    return false;
}

/*--------------------------------------------------------------------------------*/
// Tau trigger matching
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::matchTauTriggers()
{
//-DG--  if(m_dbg>=5) cout << "matchTauTriggers" << endl;
//-DG--  for(uint i=0; i<m_preTaus.size(); i++){
//-DG--
//-DG--    int iTau = m_preTaus[i];
//-DG--    const TLorentzVector &lv = m_susyObj[m_eleIDDefault]->GetTauTLV(iTau);
//-DG--
//-DG--    // trigger flags
//-DG--    long long flags = 0;
//-DG--
//-DG--    // tau20_medium1
//-DG--    if( matchTauTrigger(lv, m_event.trig_EF_tau.EF_tau20_medium1()) ){
//-DG--      flags |= TRIG_tau20_medium1;
//-DG--    }
//-DG--    // tau20Ti_medium1
//-DG--    if( matchTauTrigger(lv, m_event.trig_EF_tau.EF_tau20Ti_medium1()) ){
//-DG--      flags |= TRIG_tau20Ti_medium1;
//-DG--    }
//-DG--    // tau29Ti_medium1
//-DG--    if( matchTauTrigger(lv, m_event.trig_EF_tau.EF_tau29Ti_medium1()) ){
//-DG--      flags |= TRIG_tau29Ti_medium1;
//-DG--    }
//-DG--    // tau29Ti_medium1_tau20Ti_medium1
//-DG--    if( matchTauTrigger(lv, m_event.trig_EF_tau.EF_tau29Ti_medium1_tau20Ti_medium1()) ){
//-DG--      flags |= TRIG_tau29Ti_medium1_tau20Ti_medium1;
//-DG--    }
//-DG--    // tau20Ti_medium1_e18vh_medium1
//-DG--    if( matchTauTrigger(lv, m_event.trig_EF_tau.EF_tau20Ti_medium1_e18vh_medium1()) ){
//-DG--      flags |= TRIG_tau20Ti_medium1_e18vh_medium1;
//-DG--    }
//-DG--    // tau20_medium1_mu15
//-DG--    if( matchTauTrigger(lv, m_event.trig_EF_tau.EF_tau20_medium1_mu15()) ){
//-DG--      flags |= TRIG_tau20_medium1_mu15;
//-DG--    }
//-DG--
//-DG--    // assign the trigger flags for this tau
//-DG--    m_tauTrigFlags[iTau] = flags;
//-DG--  }
}
/*--------------------------------------------------------------------------------*/
bool XaodAnalysis::matchTauTrigger(const TLorentzVector &lv, vector<int>* passTrig)
{
//-DG--  // loop over tau trigger features
//-DG--  for(int iTrig=0; iTrig < m_event.trig_EF_tau.n(); iTrig++){
//-DG--    // Check to see if this feature passed chain we want
//-DG--    if(passTrig->at(iTrig)){
//-DG--      // Now, try to match offline tau to this online tau
//-DG--      static TLorentzVector trigLV;
//-DG--      trigLV.SetPtEtaPhiM(m_event.trig_EF_tau.pt()->at(iTrig), m_event.trig_EF_tau.eta()->at(iTrig),
//-DG--                          m_event.trig_EF_tau.phi()->at(iTrig), m_event.trig_EF_tau.m()->at(iTrig));
//-DG--      float dR = lv.DeltaR(trigLV);
//-DG--      if(dR < 0.15) return true;
//-DG--    }
//-DG--  }
//-DG--  // matching failed
    return false;
}

//----------------------------------------------------------
int XaodAnalysis::truthElectronCharge(const xAOD::Electron &in)
{
    int type   = xAOD::TruthHelpers::getParticleTruthType(in);
    int origin = xAOD::TruthHelpers::getParticleTruthOrigin(in);

    if(m_dbg>15) std::cout << "check Charge flip ele " << in.pt()*MeV2GeV 
                           << " type " << type << " origin " << origin << " nTrk " << in.nTrackParticles() << endl;

    if(isPromptElectron(type,origin)){
        const xAOD::TruthParticle* truthEle = xAOD::TruthHelpers::getTruthParticle(in);
        if(m_dbg>15) std::cout << "Truth Prompt ele " << truthEle->pdgId() << " " << in.charge()  << endl;
        if (truthEle->pdgId()==11) return -1;
        if (truthEle->pdgId()==-11) return 1;
    }
    else{        
        if(type==4 && origin==5 && in.nTrackParticles()>1){//Not sure if Type/Origin check really needed
            for(uint itrk=0; itrk< in.nTrackParticles(); itrk++ ){
                const xAOD::TrackParticle* trackParticle = in.trackParticle(itrk);
                static SG::AuxElement::Accessor<int> acc_truthType("truthType");
                static SG::AuxElement::Accessor<int> acc_truthOrigin("truthOrigin");
                if(m_dbg>15) std::cout << "\tNon-prompt ele " << trackParticle << endl;
                if(trackParticle){//Just in case track got lost in slimming.
                    int extraType=0;
                    int extraOrigin=0;
                    if(acc_truthType.isAvailable(*trackParticle))   extraType   = acc_truthType(*trackParticle);
                    if(acc_truthOrigin.isAvailable(*trackParticle)) extraOrigin = acc_truthOrigin(*trackParticle);
                    if(m_dbg>15) std::cout << "\t\t pt " << trackParticle->pt()*MeV2GeV
                                           << " type " << extraType << " & origin " << extraOrigin << endl;
                    if(isPromptElectron(extraType,extraOrigin)) {
                        const xAOD::TruthParticle* truthEle = xAOD::TruthHelpers::getTruthParticle(*trackParticle);
                        if(m_dbg>15)
                            std::cout << " \t\t Found charged flipped ?" << truthEle->pdgId() << " " << in.charge() << endl;
                        if (truthEle->pdgId()==11) return -1;
                        if (truthEle->pdgId()==-11) return 1;
                    }
                }
                else {
                    if(m_dbg>15) std::cout << "\t Don't have track " << endl;
                }
            }
        }
    }
    if(m_dbg>15) std::cout << "Cannot determined charge " << std::endl;
    return 0;
}

//----------------------------------------------------------
bool XaodAnalysis::isChargeFlip(int recoCharge, int truthCharge)
{
    return ( (recoCharge*truthCharge) > 0);
}
//----------------------------------------------------------
int XaodAnalysis::classifyTau(const xAOD::TauJet &in)
{
    int type=-1; 
    int origin=-1;
    
    for(uint itrk=0; itrk< in.nTracks(); itrk++){
            const xAOD::TrackParticle* trackParticle = in.track( itrk );
            static SG::AuxElement::Accessor<int> acc_truthType("truthType");
            static SG::AuxElement::Accessor<int> acc_truthOrigin("truthOrigin");
            if(acc_truthType.isAvailable(*trackParticle))   
                type   = acc_truthType(*trackParticle);
            if(acc_truthOrigin.isAvailable(*trackParticle)) 
                origin = acc_truthOrigin(*trackParticle);
            if(m_dbg>15) std::cout<<"Tau object. itrk= "<<itrk<<" trackParticle= "<<trackParticle
                                  <<" Type= "<<type<<" Origin= "<<origin<<std::endl; 
            // electron faking a tau // (do we need muon faking tau? never saw that)
            if (isPromptElectron(type,origin)) return ELEC; 
            // electron conversion faking a tau // maybe just origin==5 ?
            else if (type==4&&(origin==5||origin==6||origin==7)) return CONV; 
    }
            
// REL20
/*
   int JetPartonID = (in.jet())->auxdata< int >("PartonTruthLabelID"); // ghost association
   int JetConeID   = (in.jet())->auxdata< int >("ConeTruthLabelID"); // cone association
    
   if (JetConeID ==15) return PROMPT; // 
   if (JetConeID ==4 || JetConeID ==5 ) return HF;
   if (JetPartonID ==1 || JetPartonID ==2 ||  JetPartonID ==3) return QJET; // Light Quark jet
   if (JetPartonID ==21) return GJET; // Gluon Jet
   return UK;
*/        
//REL19 -> can not separate Q and G, put everything in UK
   int JetConeID = (in.jet())->auxdata< int >("TruthLabelID");
   if (JetConeID ==15) return PROMPT; 
   if (JetConeID ==4 || JetConeID ==5 ) return HF;
   return UK;
}

//----------------------------------------------------------
bool XaodAnalysis::passGRL(const xAOD::EventInfo* eventinfo)
{
    return (m_isMC ||
            m_grl->passRunLB(eventinfo->runNumber(), eventinfo->lumiBlock()));
}
//----------------------------------------------------------
bool XaodAnalysis::passTTCVeto(const xAOD::EventInfo* eventinfo)
{
    bool eventPassesTTC = eventinfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ? false : true;
    return eventPassesTTC;
    //   return (m_event.eventinfo.coreFlags() & 0x40000) == 0;
}
//----------------------------------------------------------
bool XaodAnalysis::passTileErr(const xAOD::EventInfo* eventinfo)
{
    // dantrim - add check of MC? >> seems to catch this, MC set to true
    bool eventPassesTileTrip = eventinfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ? false : true;
//	bool eventPassesTileTrip = (m_isMC || true); // SUSYToolsTester: move to xAOD tool
    return eventPassesTileTrip;
}
bool XaodAnalysis::passSCTErr(const xAOD::EventInfo* eventinfo)
{
    bool passSCTerr = eventinfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error ? false : true;
    return passSCTerr;
}
//----------------------------------------------------------
bool XaodAnalysis::passLarErr(const xAOD::EventInfo* eventinfo)
{
    // dantrim - add check of MC? >> seems to catch this, MC set to true
    bool eventPassesLarErr = eventinfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ? false : true;
    return eventPassesLarErr;
//    return m_isMC || (m_event.eventinfo.larError()!=2);
}
/*--------------------------------------------------------------------------------*/
// Check event level cleaning cuts like GRL, LarError, etc.
/*--------------------------------------------------------------------------------*/
void XaodAnalysis::assignEventCleaningFlags()
{
    const xAOD::EventInfo* eventinfo = xaodEventInfo();
    if(passGRL(eventinfo))            m_cutFlags |= ECut_GRL;
    if(passTTCVeto(eventinfo))        m_cutFlags |= ECut_TTC;
    if(passLarErr(eventinfo))         m_cutFlags |= ECut_LarErr; 
    if(passTileErr(eventinfo))        m_cutFlags |= ECut_TileErr;
    if(passSCTErr(eventinfo))         m_cutFlags |= ECut_SCTErr;
    if(passGoodVtx())                 m_cutFlags |= ECut_GoodVtx;
    if(passTileTrip())                m_cutFlags |= ECut_TileTrip;
}
//----------------------------------------------------------
void XaodAnalysis::assignObjectCleaningFlags(ST::SystInfo sysInfo, SusyNtSys sys)
{
    //AT check if m_cutFalgs save at each systematics ?
    if(passTileHotSpot())             m_cutFlags |= ECut_HotSpot;
    if(passBadJet(sysInfo, sys))      m_cutFlags |= ECut_BadJet;
    if(passBadMuon(sysInfo,sys))      m_cutFlags |= ECut_BadMuon;
    if(passCosmic(sysInfo,sys))       m_cutFlags |= ECut_Cosmic;
    if(passLarHoleVeto())             m_cutFlags |= ECut_SmartVeto;
}
//----------------------------------------------------------
bool XaodAnalysis::passLarHoleVeto()
{
    // LAr veto is not used anymore
    return true;
}
//----------------------------------------------------------
bool XaodAnalysis::passTileHotSpot()
{
    return false;
//-DG--  xAOD::JetContainer *jets =  xaodJets();
//-DG--  return !check_jet_tileHotSpot(jets, m_preJets, m_susyObj, !m_isMC, m_event.eventinfo.RunNumber());
}
//----------------------------------------------------------
bool XaodAnalysis::passBadJet(ST::SystInfo sysInfo, SusyNtSys sys)
{
    xAOD::JetContainer* jets = xaodJets(sysInfo, sys);
    bool pass_jetCleaning = true;
    for(auto &i : m_preJets) { 
        if((bool)jets->at(i)->auxdata< char >("bad")==1) pass_jetCleaning = false;
    } 
    return pass_jetCleaning;
//  return !IsBadJetEvent(jets, m_baseJets, 20.*GeV, m_susyObj);
}
//----------------------------------------------------------
bool XaodAnalysis::passGoodVtx()
{
    for(auto it=xaodVertices()->begin(), end=xaodVertices()->end(); it!=end; ++it){
        const xAOD::Vertex *vtx = *it;
        if(vtx->vertexType() == xAOD::VxType::PriVtx) return true;
    }
    return false;
}
//----------------------------------------------------------
bool XaodAnalysis::passTileTrip()
{
    return false;
//  return !m_susyObj[m_eleIDDefault]->IsTileTrip(m_event.eventinfo.RunNumber(), m_event.eventinfo.lbn(), m_event.eventinfo.EventNumber());
}
//----------------------------------------------------------
bool XaodAnalysis::passBadMuon(ST::SystInfo sysInfo, SusyNtSys sys)
{
    xAOD::MuonContainer* muons = xaodMuons(sysInfo, sys);
    bool pass_bad_muon = true;
    for(auto &i : m_preMuons) {
        if(muons->at(i)->auxdata<char>("baseline")==1 && muons->at(i)->auxdata<char>("bad")==1) { 
            pass_bad_muon = false;
        }
    }
   // for(auto &i : m_baseMuons) {
   //     if(m_susyObj[m_eleIDDefault]->IsBadMuon(*muons->at(i))) pass_bad_muon = false;
   // }
    return pass_bad_muon;
    //  return !IsBadMuonEvent(m_susyObj, xaodMuons(), m_preMuons, 0.2);
}
//----------------------------------------------------------
bool XaodAnalysis::passCosmic(ST::SystInfo sysInfo, SusyNtSys sys)
{
    xAOD::MuonContainer* muons = xaodMuons(sysInfo, sys);
    bool pass_cosmic = true;
    for(auto &i : m_baseMuons) {
        if((bool)muons->at(i)->auxdata< char >("passOR")==1 && (bool)muons->at(i)->auxdata< char >("cosmic")==1) pass_cosmic = false;
    }
    return pass_cosmic;


//    for(auto it=muons->begin(), end=muons->end(); it!=end; ++it){
//        const xAOD::Muon &mu = **it;
//        if(!mu.auxdata< bool >("baseline")) continue; // dantrim -- removing this, cutflow is in line with Maria : Feb 14 2015
//        if(!mu.auxdata< bool >("passOR")) continue; 
//        if(m_susyObj[m_eleIDDefault]->IsCosmicMuon(mu)) return false;
//
//     //   return(mu.auxdata<bool>("baseline") && mu.auxdata<bool>("passOR") && !m_susyObj[m_eleIDDefault]->IsCosmicMuon(mu));
//    }
//     return true;
    //  return !IsCosmic(m_susyObj, xaodMuons(), m_baseMuons, 1., 0.2);
}
//----------------------------------------------------------
float XaodAnalysis::getEventWeight(float lumi)
{
    if(!m_isMC) return 1;
    else
        return 1.0;
//  return m_event.eventinfo.mc_event_weight() * getXsecWeight() * getPileupWeight() * lumi / m_sumw;
}
//----------------------------------------------------------
float XaodAnalysis::getXsecWeight()
{
    return 1.0;
    // // Use user cross section if it has been set
    // if(m_xsec > 0) return m_xsec;

    // // Use SUSY cross section file
    // int id = m_event.eventinfo.mc_channel_number();
    // if(m_xsecMap.find(id) == m_xsecMap.end()) {
    //   m_xsecMap[id] = m_susyXsec->process(id);
    // }
    // return m_xsecMap[id].xsect() * m_xsecMap[id].kfactor() * m_xsecMap[id].efficiency();
}
//----------------------------------------------------------
double XaodAnalysis::getPileupWeight(const xAOD::EventInfo* eventinfo)
{
    double pile_up_w = 1.0;
    if(eventinfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ) {
        pile_up_w = m_pileupReweightingTool->getCombinedWeight(*eventinfo);
    //    pile_up_w = m_pileupReweightingTool->GetCombinedWeight(eventinfo->runNumber(),
    //                    eventinfo->mcChannelNumber(), eventinfo->averageInteractionsPerCrossing());
    }
    return pile_up_w;

    //m_pileupReweightingTool->execute();
    //return xaodEventInfo()->auxdata< double >( "PileupWeight" );
}
//----------------------------------------------------------
float XaodAnalysis::getPileupWeightUp()
{
    return 1.0;
    // return m_pileup_up->GetCombinedWeight(m_event.eventinfo.RunNumber(), m_event.eventinfo.mc_channel_number(), m_event.eventinfo.averageIntPerXing());
}
//----------------------------------------------------------
float XaodAnalysis::getPileupWeightDown()
{
    return 1.0;
    // return m_pileup_dn->GetCombinedWeight(m_event.eventinfo.RunNumber(), m_event.eventinfo.mc_channel_number(), m_event.eventinfo.averageIntPerXing());
}
//----------------------------------------------------------
float XaodAnalysis::getPDFWeight8TeV()
{
    return 1.0;
    // #ifdef USEPDFTOOL
    // float scale = m_event.mcevt.pdf_scale()->at(0);
    // float x1 = m_event.mcevt.pdf_x1()->at(0);
    // float x2 = m_event.mcevt.pdf_x2()->at(0);
    // int id1 = m_event.mcevt.pdf_id1()->at(0);
    // int id2 = m_event.mcevt.pdf_id2()->at(0);

    // // MultLeip function... Not working?
    // //return scaleBeamEnergy(*m_pdfTool, 21000, m_event.mcevt.pdf_scale()->at(0), m_event.mcevt.pdf_x1()->at(0),
    //                        //m_event.mcevt.pdf_x2()->at(0), m_event.mcevt.pdf_id1()->at(0), m_event.mcevt.pdf_id2()->at(0));
    // // Simple scaling
    // //return m_pdfTool->event_weight( pow(scale,2), x1, x2, id1, id2, 21000 );

    // // For scaling to/from arbitrary beam energy
    // m_pdfTool->setEventInfo( scale*scale, x1, x2, id1, id2 );
    // //return m_pdfTool->scale((3.5+4.)/3.5);
    // // possible typo correction?
    // return m_pdfTool->scale(4./3.5);

    // #else
    // return 1;
    // #endif
}
//----------------------------------------------------------
float XaodAnalysis::getLepSF(const vector<LeptonInfo>& leptons)
{
    // TODO: incorporate systematics
    float lepSF = 1;
    if(m_isMC){
        // // Loop over leptons
        // for(uint iLep=0; iLep<leptons.size(); iLep++){
        //   const LeptonInfo &lep = leptons[iLep];
        //   // Electrons
        //   if(lep.isElectron()){
        //       const xAOD::ElectronD3PDObjectElement* el = lep.getElectronElement();
        //     lepSF *= m_susyObj[m_eleIDDefault]->GetSignalElecSF(el->cl_eta(), lep.lv()->Pt(), true, true, false);
        //   }
        //   // Muons
        //   else{
        //     lepSF *= m_susyObj[m_eleIDDefault]->GetSignalMuonSF(lep.idx());
        //   }
        // }
    }
    return lepSF;
}
//----------------------------------------------------------
float XaodAnalysis::getBTagSF(const vector<int>& jets)
{
    return 1;
}
//----------------------------------------------------------
void XaodAnalysis::calcRandomRunLB()
{
//-DG--  if(m_pileup){
//-DG--    m_mcRun = m_pileup->GetRandomRunNumber(m_event.eventinfo.RunNumber());
//-DG--    m_mcLB = m_pileup->GetRandomLumiBlockNumber(m_mcRun);
//-DG--  }
}
//----------------------------------------------------------
int XaodAnalysis::getHFORDecision()
{
    //AT-2014-10-30 Not yet implemented for DC14
    return 1;
    // return m_hforTool.getDecision(m_event.eventinfo.mc_channel_number(),
    //                               m_event.mc.n(),
    //                               m_event.mc.pt(),
    //                               m_event.mc.eta(),
    //                               m_event.mc.phi(),
    //                               m_event.mc.m(),
    //                               m_event.mc.pdgId(),
    //                               m_event.mc.status(),
    //                               m_event.mc.vx_barcode(),
    //                               m_event.mc.parent_index(),
    //                               m_event.mc.child_index(),
    //                               HforToolD3PD::ALL); //HforToolD3PD::DEFAULT
}
//----------------------------------------------------------
void XaodAnalysis::dumpEvent()
{
    const xAOD::EventInfo* eventinfo = xaodEventInfo();
    cout<<(*eventinfo)<<endl;
}
//----------------------------------------------------------
void XaodAnalysis::dumpBaselineObjects()
{
    uint nEle = m_baseElectrons.size();
    uint nMu  = m_baseMuons.size();
    //uint nTau = m_baseTaus.size();
    uint nJet = m_baseJets.size();
    ST::SystInfo sysInfo = systInfoList[0]; // nominal
    cout.precision(2);
    if(nEle){
        xAOD::ElectronContainer* electrons = XaodAnalysis::xaodElectrons(sysInfo);
        cout << "Baseline electrons" << endl;
        for(auto& iEl : m_baseElectrons){
            const xAOD::Electron* el = electrons->at(iEl);
            cout << "  El : " << fixed
                 << " q " << setw(2) << (int) el->charge()
                 << " pt " << setw(6) << el->p4().Pt()/1000. //lv.Pt()/GeV
                 << " eta " << setw(5) << el->p4().Eta() //lv.Eta()
                 << " phi " << setw(5) << el->p4().Phi(); //lv.Phi();
            if(m_isMC) cout << " type " << setw(2) << xAOD::TruthHelpers::getParticleTruthType(*el)
                         << " origin " << setw(2)  << xAOD::TruthHelpers::getParticleTruthOrigin(*el);
            cout << endl;
        } // El
    } // nEle
    if(nMu){
        xAOD::MuonContainer* muons = XaodAnalysis::xaodMuons(sysInfo);
        cout << "Baseline muons" << endl;
        for(auto& iMu : m_baseMuons){
            const xAOD::Muon* mu = muons->at(iMu);
            cout << "  Mu : " << fixed
                 << " q " << setw(2) << (int) mu->charge() //(int) muo.charge()
                 << " pt " << setw(6) << mu->p4().Pt()/1000. //lv.Pt()/GeV
                 << " eta " << setw(5) << mu->p4().Eta() //lv.Eta()
                 << " phi " << setw(5) << mu->p4().Phi(); //lv.Phi();
            if(m_isMC) {
                const xAOD::TrackParticle* trackParticle = mu->primaryTrackParticle();
                if(trackParticle) {
                    static SG::AuxElement::Accessor<int> acc_truthType("truthType");
                    static SG::AuxElement::Accessor<int> acc_truthOrigin("truthOrigin");
                    int truthType = -999999, truthOrigin = -999999;
                    if(acc_truthType.isAvailable(*trackParticle)) truthType = acc_truthType(*trackParticle);
                    if(acc_truthOrigin.isAvailable(*trackParticle)) truthOrigin = acc_truthOrigin(*trackParticle);
                    cout << " type " << setw(2) << truthType << " origin " << setw(2) << truthOrigin;
                    cout << endl;
                } // if trackparticle
            } // if is MC
        } // iMu
    } // nMu
    if(nJet){
        xAOD::JetContainer* jets = XaodAnalysis::xaodJets(sysInfo);
        cout << "Baseline jets" << endl;
        for(auto& iJ : m_baseJets) {
            xAOD::Jet* jet = jets->at(iJ);
            cout << "  Jet : " << fixed
            << " pt " << setw(6) << jet->p4().Pt()/1000. //lv.Pt()/GeV
            << " eta " << setw(5) << jet->p4().Eta() //lv.Eta()
            << " phi " << setw(5) << jet->p4().Phi() //lv.Phi()
            << " mv1 " << (jet->btagging())->MV1_discriminant(); //jet.flavor_weight_MV1();
            cout << endl;
        } // iJ
    } // nJet
cout.precision(6);
cout.unsetf(ios_base::fixed);
}
//----------------------------------------------------------
void XaodAnalysis::dumpSignalObjects()
{
    uint nEle = m_sigElectrons.size();
    uint nMu  = m_sigMuons.size();
    // taus
    //uint nTau = m_sigTaus.size();
    uint nJet = m_sigJets.size();
    ST::SystInfo sysInfo = systInfoList[0]; // nominal

    cout.precision(2);
    if(nEle){
        xAOD::ElectronContainer* electrons = XaodAnalysis::xaodElectrons(sysInfo);
        cout << "Signal electrons" << endl;
        for(auto& iEl : m_sigElectrons) {
            const xAOD::Electron* el = electrons->at(iEl); 
            cout << "  El : " << fixed
                 << " q " << setw(2) << (int) el->charge()
                 << " pt " << setw(6) << el->p4().Pt()/1000. //lv.Pt()/GeV
                 << " eta " << setw(5) << el->p4().Eta() //lv.Eta()
                 << " phi " << setw(5) << el->p4().Phi(); //lv.Phi();
                if(m_isMC) cout << " type " << setw(2)   << xAOD::TruthHelpers::getParticleTruthType(*el)
                                << " origin " << setw(2) << xAOD::TruthHelpers::getParticleTruthOrigin(*el);
            cout << endl;
        } // iEl
    }
    if(nMu){
        xAOD::MuonContainer* muons = XaodAnalysis::xaodMuons(sysInfo);
        cout << "Signal muons" << endl;
        for(auto& iMu : m_sigMuons){
            const xAOD::Muon* mu = muons->at(iMu);
            cout << "  Mu : " << fixed
                 << " q " << setw(2) << (int) mu->charge() //(int) muo.charge()
                 << " pt " << setw(6) << mu->p4().Pt()/1000. //lv.Pt()/GeV
                 << " eta " << setw(5) << mu->p4().Eta() //lv.Eta()
                 << " phi " << setw(5) << mu->p4().Phi(); //lv.Phi();
            if(m_isMC) {
                const xAOD::TrackParticle* trackParticle = mu->primaryTrackParticle();
                if(trackParticle) {
                    static SG::AuxElement::Accessor<int> acc_truthType("truthType");
                    static SG::AuxElement::Accessor<int> acc_truthOrigin("truthOrigin");
                    int truthType = -999999, truthOrigin = -999999;
                    if(acc_truthType.isAvailable(*trackParticle)) truthType = acc_truthType(*trackParticle);
                    if(acc_truthOrigin.isAvailable(*trackParticle)) truthOrigin = acc_truthOrigin(*trackParticle);
                    cout << " type " << setw(2) << truthType << " origin " << setw(2) << truthOrigin;
                    cout << endl;
                } // if trackparticle
            } // if is MC
        } // iMu
    } // nMu
    if(nJet){
        xAOD::JetContainer* jets = XaodAnalysis::xaodJets(sysInfo);
        cout << "Signal jets" << endl;
        for(auto& iJ : m_sigJets) {
            xAOD::Jet* jet = jets->at(iJ);
            cout << "  Jet : " << fixed
            << " pt " << setw(6) << jet->p4().Pt()/1000. //lv.Pt()/GeV
            << " eta " << setw(5) << jet->p4().Eta() //lv.Eta()
            << " phi " << setw(5) << jet->p4().Phi() //lv.Phi()
            << " mv1 " << (jet->btagging())->MV1_discriminant(); //jet.flavor_weight_MV1();
            cout << endl;
        } // iJ
    } // nJet
cout.precision(6);
cout.unsetf(ios_base::fixed);
}
//----------------------------------------------------------
bool XaodAnalysis::runningOptionsAreValid()
{
    bool valid=true;
    bool isSimulation = xaodEventInfo()->eventType( xAOD::EventInfo::IS_SIMULATION );
    bool isData = !isSimulation;
    if(m_isMC != isSimulation) {
        valid=false;
        if(m_dbg)
            cout<<"XaodAnalysis::runningOptionsAreValid invalid isMC:"
                <<" (m_isMC: "<<m_isMC<<" != isSimulation: "<<isSimulation<<")"
                <<endl;
    }
    if(isData) { // check stream
        const std::vector< xAOD::EventInfo::StreamTag > &streams= xaodEventInfo()->streamTags();
        vector<string> streamnames(streams.size());
        std::transform(streams.begin(), streams.end(), streamnames.begin(),
                       [](const xAOD::EventInfo::StreamTag &s) {
                           cout << "AT:  stream " << s.name()<< endl;
                           return s.name();
                       });
        bool isPhysicsMain = (find(streamnames.begin(), streamnames.end(), "Main") != streamnames.end());
        bool consistentStream = (isPhysicsMain ? m_stream==Stream_PhysicsMain : false); 

        if(!consistentStream) {
            valid=false;
            if(m_dbg)
                cout << "XaodAnalysis::runningOptionsAreValid: inconsistent stream"
                     << " m_stream: "
                     << (m_stream==Stream_PhysicsMain   ? "Stream_PhysicsMain" : "unknown")
                     << " eventinfo: "
                     << accumulate(streamnames.begin(), streamnames.end(), std::string(),
                            [](const std::string& a, const std::string& b) -> std::string {
                                return a + (a.length() > 0 ? "," : "") + b;
                                })
                     << endl;

        } // !consistentStream
    } // isData
    if(m_dbg)
        cout<<"XaodAnalysis::runningOptionsAreValid(): "<<(valid?"true":"false")<<endl;
    return valid;
}
//----------------------------------------------------------
std::string XaodAnalysis::defaultGrlFile()
{
    string grl_file = "";
    if(m_isData15) {
        grl_file = "$ROOTCOREBIN/data/SusyCommon/data15_13TeV.periodAllYear_DetStatus-v79-repro20-02_DQDefects-00-02-02_PHYS_StandardGRL_All_Good_25ns.xml";
    }
    else if(m_isData16) {
        grl_file = "$ROOTCOREBIN/data/SusyCommon/data16_13TeV.periodAllYear_DetStatus-v82-pro20-12_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml";
    }
    else {
        cout << "XaodAnalysis::defaultGrlFile    ERROR Inconsistent data flags. Neither \"m_isData15\" nor \"m_isData15\" flags are set!" << endl;
        exit(1);
    }
    return grl_file;
}
//----------------------------------------------------------
bool XaodAnalysis::initGrlTool()
{
    bool success = true;
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    std::vector<std::string> grl_files;
    grl_files.push_back(XaodAnalysis::defaultGrlFile());
    cout << "XaodAnalysis::initGrlTool    Using GRL: " << grl_files[0] << endl;
    m_grl->setProperty("GoodRunsListVec", grl_files);
    m_grl->setProperty("PassThrough", false);
    if(!m_grl->initialize().isSuccess()) success = false;
    return success;
}
//----------------------------------------------------------
DataStream XaodAnalysis::streamFromSamplename(const TString &sample, bool isMC)
{
    bool isData(!isMC);
//    TString sample(s.c_str());
    DataStream stream = Stream_Unknown;
    if     (isMC) stream = Stream_MC;
    else if(sample.Contains("main", TString::kIgnoreCase) && isData) stream = Stream_PhysicsMain;
    else
        cout<<"XaodAnalysis::streamFromSamplename('"<<sample<<"',isData="<<(isData?"true":"false")<<")"
            <<" : cannot determine the stream, returning "<<streamName(stream)<<endl;
    return stream;
}
//----------------------------------------------------------
bool XaodAnalysis::isDataFromSamplename(const TString &sample)
{
    bool is_data = false;
    bool is_data16 = false;
    is_data = sample.Contains("data", TString::kIgnoreCase);
    is_data16 = sample.Contains("data16", TString::kIgnoreCase);

    if(is_data && !is_data16) m_isData15 = true;
    else if(is_data && is_data16) m_isData16 = true;

    return is_data;
}
//----------------------------------------------------------
bool XaodAnalysis::isSimuFromSamplename(const TString &s)
{
    bool isMCsample = !XaodAnalysis::isDataFromSamplename(s);
    stringstream sx;
    sx << "XaodAnalysis::isSimuFromSamplename    "
       << "Sample (" << s <<") treated as simulation? " << (isMCsample ? "YES" : "NO") << endl;
    if(!isMCsample) {
        sx << "XaodAnalysis::isSimuFromSamplename    "
           << " > " << (m_isData16 ? " data16_13Tev " : " data15_13TeV ") << endl;
    }
    return isMCsample;
}
//----------------------------------------------------------
bool XaodAnalysis::isDerivationFromMetaData(TTree* intree, bool verbose)
{
    // Following implementation in SUSYToolsTester
    bool is_derivation = false;
    TTree* metadata = nullptr;

    if(TDirectory* treeDir = getDirectoryFromTreeOrChain(intree, verbose)){
        if(treeDir->Get("MetaData")){
            metadata = dynamic_cast<TTree*>(treeDir->Get("MetaData"));
        }
    }
    if(metadata){
        metadata->LoadTree(0);
        is_derivation = !metadata->GetBranch("StreamAOD");
    } else {
       cout << "XaodAnalysis::isDerivationFromMetaData    cannot get MetaData tree" << endl;
    }

    if(is_derivation) cout << "Treating input as a derivation" << endl;

    return is_derivation;
/*
    if(TDirectory* treeDir = getDirectoryFromTreeOrChain(intree, verbose)){
        if(TObject *obj = treeDir->Get("MetaData")){
            metadata = static_cast<TTree*>(obj);
        }
    }
    if(metadata){
        TTreeFormula streamAOD("StreamAOD", "StreamAOD", metadata);
        // DG 2015-05-17 don't understand the logic here, should clarify...
        // why using a TTreeFormula that will cause warning msgs? can
        // we just check whether the branch is there?
        is_derived = (streamAOD.GetNcodes() < 1);
        if(verbose) cout<<"This file is "<<(is_derived ? "" : "NOT ")<<"a derivation"<<endl;
    } else {
        cout<<"XaodAnalysis::isDerivationFromMetaData: cannot get metadata tree"<<endl;
    }
    return is_derived;
*/
}
//----------------------------------------------------------
bool XaodAnalysis::getCutBookkeeperInfo(xAOD::TEvent& event)
{
    bool ok = true;
    const xAOD::CutBookkeeperContainer* completeCBC = 0;
    ok = event.retrieveMetaInput(completeCBC, "CutBookkeepers");
    if(!ok) {
        cout << "XaodAnalysis::getCutBookkeeperInfo    "
             << "Failed to get CBK metadata! Exitting." << endl;
        exit(1);
    }

    const xAOD::CutBookkeeper* allEventsCBK = 0;
    int maxcycle = -1;
    for(auto cbk : *completeCBC) {
        if(cbk->name() == "AllExecutedEvents" && cbk->inputStream()=="StreamAOD" && cbk->cycle() > maxcycle) {
            maxcycle = cbk->cycle();
            allEventsCBK = cbk;
        } // if
    } // cbk
    uint64_t nevents_ = -1;
    double sumw_ = -1.0;
    double sumw2_ = -1.0;
    if(allEventsCBK) {
        nevents_ = allEventsCBK->nAcceptedEvents();
        sumw_ = allEventsCBK->sumOfEventWeights();
        sumw2_ = allEventsCBK->sumOfEventWeightsSquared();
    } // all
    else {
        cout << "XaodAnalysis::getCutBookkeeperInfo    "
             << "\"AllExecutedEvents\" branch not found in CBK metadata! Exitting." << endl;
        exit(1);
    }

    m_nEventsProcessed += nevents_;
    m_sumOfWeights += sumw_;
    m_sumOfWeightsSquared += sumw2_;
    cout << "nAcceptedEvents: " << nevents_ << "  sumOfEventWeights: " << sumw_ << "  sumOfEventWeightsSquared: " << sumw2_ << endl;

    return ok;
}
//----------------------------------------------------------
TDirectory* XaodAnalysis::getDirectoryFromTreeOrChain(TTree* tree, bool verbose)
{
    TDirectory* dir = nullptr;
    if(tree){
        dir = tree->GetDirectory(); // probably a file
        if(dir){
            if(verbose) cout<<"got the directory from the tree : "<<dir->GetName()<<endl;
        } else {
            if(verbose) cout<<"trying to get the directory from a chain"<<endl;
            if(TChain *c = dynamic_cast<TChain*>(tree)){
                if(TChainElement* ce = static_cast<TChainElement*>(c->GetListOfFiles()->First())){
                    TFile *firstFile = TFile::Open(ce->GetTitle()); // not owned (TChain will close it?), see TChain::GetListOfFiles
                    dir = static_cast<TDirectory*>(firstFile);
                }
            }
        }
    }
    if(verbose)
        cout<<"getDirectoryFromTreeOrChain: got "<<(dir ? dir->GetName() : "NULL")<<endl;
    return dir;
}
//----------------------------------------------------------
void XaodAnalysis::selectObjects(SusyNtSys sys, ST::SystInfo sysInfo)
{
    selectBaselineObjects(sys, sysInfo);
    selectSignalObjects(sys,sysInfo);
//--DG-- todo     if(m_selectTruth) selectTruthObjects();
}
//----------------------------------------------------------
XaodAnalysis& XaodAnalysis::deleteShallowCopies(bool deleteNominal)
{
    if(m_dbg>5) cout << "deleteShallowCopies " << deleteNominal << endl;

    if(m_xaodMuons        ) delete m_xaodMuons;
    if(m_xaodMuonsAux     ) delete m_xaodMuonsAux;
    if(m_xaodElectrons    ) delete m_xaodElectrons;
    if(m_xaodElectronsAux ) delete m_xaodElectronsAux;
    if(m_xaodTaus         ) delete m_xaodTaus;
    if(m_xaodTausAux      ) delete m_xaodTausAux;
    if(m_xaodJets         ) delete m_xaodJets;
    if(m_xaodJetsAux      ) delete m_xaodJetsAux;
    if(m_xaodPhotons      ) delete m_xaodPhotons;
    if(m_xaodPhotonsAux   ) delete m_xaodPhotonsAux;

    if(m_metContainer     ) delete m_metContainer;
    if(m_metAuxContainer  ) delete m_metAuxContainer;
    if(m_trackMetContainer ) delete m_trackMetContainer;
    if(m_trackMetAuxContainer ) delete m_trackMetAuxContainer;

    if(m_dbg>5) cout << "Check delete shallowCopied mu " << m_xaodMuons
                     << " ele "  << m_xaodElectrons
                     << " pho "  << m_xaodPhotons
                     << " jets " << m_xaodJets
                     << " taus " << m_xaodTaus 
                     << " met "  << m_metContainer
                     << endl;
    
    if(deleteNominal){
        //m_store.print();
        m_store.clear(); // this clears m_trackMetContainer, m_xaodTruthEvent and m_xaodTruthParticles
        //and any objs recorded with TStore 

        if(m_xaodMuons_nom         ) delete m_xaodMuons_nom;
        if(m_xaodMuonsAux_nom      ) delete m_xaodMuonsAux_nom;
        if(m_xaodElectrons_nom     ) delete m_xaodElectrons_nom;
        if(m_xaodElectronsAux_nom  ) delete m_xaodElectronsAux_nom;
        if(m_xaodTaus_nom          ) delete m_xaodTaus_nom;
        if(m_xaodTausAux_nom       ) delete m_xaodTausAux_nom;
        if(m_xaodJets_nom          ) delete m_xaodJets_nom;
        if(m_xaodJetsAux_nom       ) delete m_xaodJetsAux_nom;
        if(m_xaodPhotons_nom       ) delete m_xaodPhotons_nom;
        if(m_xaodPhotonsAux_nom    ) delete m_xaodPhotonsAux_nom;

        if(m_xaodTruthParticlesAux ) delete m_xaodTruthParticlesAux;
    }
    
    clearContainerPointers(deleteNominal);

    return *this;
}
//----------------------------------------------------------
XaodAnalysis& XaodAnalysis::clearContainerPointers(bool deleteNominal)
{
    //Clear the pointer of the container that are effected by systematics
    m_xaodMuons          = NULL;
    m_xaodMuonsAux       = NULL;
    m_xaodElectrons      = NULL;
    m_xaodElectronsAux   = NULL;
    m_xaodTaus           = NULL;
    m_xaodTausAux        = NULL;
    m_xaodJets           = NULL;
    m_xaodJetsAux        = NULL;
    m_xaodPhotons        = NULL;
    m_xaodPhotonsAux     = NULL;

    m_metContainer           = NULL;
    m_metAuxContainer        = NULL;
    m_trackMetContainer      = NULL;
    m_trackMetAuxContainer   = NULL;

    if(deleteNominal){
        m_xaodMuons_nom          = NULL;
        m_xaodMuonsAux_nom       = NULL;
        m_xaodElectrons_nom      = NULL;
        m_xaodElectronsAux_nom   = NULL;
        m_xaodTaus_nom           = NULL;
        m_xaodTausAux_nom        = NULL;
        m_xaodJets_nom           = NULL;
        m_xaodJetsAux_nom        = NULL;
        m_xaodPhotons_nom        = NULL;
        m_xaodPhotonsAux_nom     = NULL;

        m_xaodTruthEvent         = NULL;
        m_xaodTruthParticles     = NULL;
        m_xaodTruthParticlesAux  = NULL; 

        m_xaodEventInfo          = NULL; 
        m_xaodVertices           = NULL; 
    }


    return *this;
}
//----------------------------------------------------------
XaodAnalysis& XaodAnalysis::retrieveCollections()
{
    if(m_dbg) cout << "XaodAnalysis::retrieveCollections " << endl;

    xaodEventInfo();
    xaodVertices();

    //Retrieve containers at nominal scale
    xaodElectrons(systInfoList[0]);
    xaodMuons(systInfoList[0]);
    xaodJets(systInfoList[0]);
    xaodTaus(systInfoList[0]);
    xaodPhotons(systInfoList[0]);
    retrieveXaodMet(systInfoList[0]);//nominal
    retrieveXaodTrackMet(systInfoList[0]);

    xaodTruthEvent();
    xaodTruthParticles();

    return *this;
}
