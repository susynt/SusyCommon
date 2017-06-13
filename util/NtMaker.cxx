//std/stl
#include <cstdlib>
#include <string>
#include <iostream>

//ROOT
#include "TChain.h"
#include "TSystem.h"

//SusyCommon
#include "SusyCommon/SusyNtMaker.h"
#include "SusyCommon/MCType.h"

//SusyNtuple
//#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"

using namespace std;
using namespace Susy;
using Susy::SusyNtMaker;

void help()
{
    cout << "--------------------------------------------------------------" << endl;
    cout << " NtMaker " << endl;
    cout << endl;
    cout << "  Options:" << endl;
    cout << "   -f|--filelist       input filelist (ROOT file, *.txt file, or dir) [REQUIRED]" << endl;
    cout << "   -i|--input          input DAOD container name [REQUIRED]" << endl;
    cout << "   -o|--output         name of output container with susyNt" << endl;
    cout << "   --outfilename       name of output ROOT file [default: susyNt.root]" << endl;
    cout << "   -d|--debug          debug level (integer values) [default: 0]" << endl;
    cout << "   -n|--nevents        number of events to process [default: -1 (all events)]" << endl;
    cout << "   -t|--tag            SusyNtuple production tag [default: ""]" << endl;
    cout << "   --mctype            MC campaign (allowed: mc15b, mc15c) [default: mc15c]" << endl;
    cout << "   --sys               run with systematics [default: false]" << endl;
    cout << "   --nLepFilter        filter on number of light leptons [default: >=1]" << endl;
    cout << "   --trigFilter        require a trigger to have fired to store event [default: false]" << endl;
    cout << "   --af2               sample is AFII reco [default: false]" << endl;
    cout << "   --nowrite           do not write out an ntuple [default: false]" << endl;
    cout << "   -s|--sys            run with systematics [default: false]" << endl;
    cout << "   --contTau           save container (loose selection) taus [default: false]" << endl;
    cout << "   --storeTruth        store truth particles in the output [default: false]" << endl;

    cout << "   -h|--help           print this help message" << endl;
    cout << "--------------------------------------------------------------" << endl;
}

int main(int argc, char** argv)
{
    string filelist = "";
    string input_container = "";
    string output_container = "";
    int dbg = 0;
    int n_events = -1;
    string nt_tag = "";
    string mc_type_str = "mc15c";
    MCType mc_type = MCType::MC15c;
    bool is_af2 = false;
    bool run_sys = false;
    bool run_lep_filter = true;
    int n_lep_filter = 1;
    bool run_trig_filter = false;
    bool write_ntuple = true;
    string outfilename = "susyNt.root";
    bool contTaus = false;
    bool storeTruth = false;

    //cout << "NtMaker" << endl;

    int optin = 1;
    while(optin < argc) {
        string in = argv[optin];
        if      (in == "-f" || in == "--filelist"   ) { filelist = argv[++optin]; }
        else if (in == "-i" || in == "--input"      ) { input_container = argv[++optin]; }
        else if (in == "-o" || in == "--output"     ) { output_container = argv[++optin]; } 
        else if (in == "--outfilename"              ) { outfilename = argv[++optin]; }
        else if (in == "-d" || in == "--debug"      ) { dbg = atoi(argv[++optin]); }
        else if (in == "-n" || in == "--nevents"    ) { n_events = atoi(argv[++optin]); }
        else if (in == "-t" || in == "--tag"        ) { nt_tag = argv[++optin]; }
        else if (in == "--mctype"                   ) { mc_type_str = argv[++optin]; } 
        else if (in == "--nLepFilter"               ) { run_lep_filter = true; n_lep_filter = atoi(argv[++optin]); }
        else if (in == "--trigFilter"               ) { run_trig_filter = true; }
        else if (in == "--af2"                      ) { is_af2 = true; }
        else if (in == "--nowrite"                  ) { write_ntuple = false; }
        else if (in == "-s" || in == "--sys"        ) { run_sys = true; }
        else if (in == "--contTau"                  ) { contTaus = true; }
        else if (in == "--storeTruth"               ) { storeTruth = true; }
        else if (in == "-h" || in == "--help"       ) { help(); return 0; }
        else {
            cout << "NtMaker    Unknown command line argument '" << in << "' provided, exiting" << endl;
            help();
            return 1;
        }
        optin++;
    } // while

    ///////////////////////////////////////////////////////////////////////////
    // provided filelist
    if(filelist=="") {
        cout << "NtMaker    No input filelist provided, exiting" << endl;
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////
    // input container name
    if(input_container=="") {
        cout << "NtMaker    Input container name not provided, exiting" << endl;
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////
    // check MCType
    vector<string> valid_mc_types = { "mc15b", "mc15c" };
    bool foundit = false;
    for(auto type : valid_mc_types) {
        if(type==mc_type_str) {
            foundit = true;
            break;
        }
    }
    if(!foundit) {
        cout << "NtMaker    Invalid MCType '" << mc_type_str << "' provided. Valid types are:";
        for(auto type : valid_mc_types) cout << " " << type;
        cout << ", exiting" << endl;
        return 1;
    }
    if(mc_type_str=="mc15b")
        mc_type = MCType::MC15b;
    else if(mc_type_str=="mc15c")
        mc_type = MCType::MC15c;

    if(n_events > 20 && dbg >=10) {
        cout << "\n---------------------------------------------------------------------" << endl;
        cout << "NtMaker    You have set more than 20 events to be processed with a high (" << dbg << ") verbosity, setting number of events to process to 20" << endl;
        cout << "\n---------------------------------------------------------------------\n" << endl;
        n_events = 20;
    }

    ///////////////////////////////////////////////////////////////////////////
    // things look ok, lets start
    cout << "-------------------------------------------------------------" << endl;
    cout << " NtMaker run conditions" << endl;
    cout << endl;
    cout << "   filelist                : " << filelist << endl;
    cout << "   input container         : " << input_container << endl;
    cout << "   output container        : " << output_container << endl;
    cout << "   output file name        : " << outfilename << endl;
    cout << "   debug level             : " << dbg << endl;
    cout << "   n events to process     : " << n_events << endl;
    cout << "   production tag          : " << nt_tag << endl;
    cout << "   MC type                 : " << mc_type_str << endl;
    cout << "   AFII                    : " << is_af2 << endl;
    cout << "   run systematics         : " << run_sys << endl;
    cout << "   run lepton filter       : " << run_lep_filter << endl;
    if(run_lep_filter)
    cout << "   n leptons to filter (>=): " << n_lep_filter << endl;
    cout << "   filter on trigger       : " << run_trig_filter << endl;
    cout << "   save container taus     : " << contTaus << endl;
    cout << "   store truth particles   : " << storeTruth << endl;
    cout << "-------------------------------------------------------------" << endl;

    if(storeTruth) {
        cout << "NtMaker    You have requsted to store TruthParticles. This is not yet implemented. Exiting." << endl;
        return 1;
    }

    TChain* chain = new TChain("CollectionTree");
    int file_err = ChainHelper::addInput(chain, filelist, dbg);
    if(file_err) return 1;
    Long64_t n_entries = chain->GetEntries();
    if(!(n_entries>0)) {
        cout << "NtMaker    TChain has zero entries, exiting" << endl;
        delete chain;
        return 0;
    }
    chain->ls();

    // build up the looper
    SusyNtMaker* susyAna = new SusyNtMaker();
    susyAna->set_debug(dbg);
    susyAna->set_chain(chain);
    if(!susyAna->set_input_container(input_container)) return 1;
    susyAna->set_output_container(output_container);
    if(!susyAna->set_mctype(mc_type)) return 1;
    susyAna->set_af2(is_af2);
    susyAna->set_write(write_ntuple);
    susyAna->run_systematics(run_sys);
    susyAna->set_output_name(outfilename);
    susyAna->set_production_tag(nt_tag);
    susyAna->set_production_command(Susy::utils::commandLineArguments(argc, argv));
    susyAna->set_nlep_filter(n_lep_filter);
    susyAna->set_trig_filter(run_trig_filter);
    susyAna->set_cont_taus(contTaus);
    

    // start up the looper
    if(n_events < 0) n_events = n_entries;
    cout << endl;
    cout << "Total entries in input chain       : " << n_entries << endl;
    cout << "Total number of entries to process : " << n_events << endl;
    chain->Process(susyAna, input_container.c_str(), n_events, 0);

    cout << endl;
    cout << "NtMaker    Job done" << endl;

    // check for expected output here?

    delete chain;
    return 0;

} // end
