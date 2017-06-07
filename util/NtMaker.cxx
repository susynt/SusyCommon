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
//using Susy::SusyNtMaker;

void help()
{
    cout << "--------------------------------------------------------------" << endl;
    cout << " NtMaker " << endl;
    cout << endl;
    cout << "  Options:" << endl;
    cout << "   -f|--filelist       input filelist (ROOT file, *.txt file, or dir) [REQUIRED]" << endl;
    cout << "   -i|--input          input DAOD container name [REQUIRED]" << endl;
    cout << "   -o|--output         name of output container with susyNt" << endl;
    cout << "   -d|--debug          debug level (integer values) [default: 0]" << endl;
    cout << "   -n|--nevents        number of events to process [default: -1 (all events)]" << endl;
    cout << "   -t|--tag            SusyNtuple production tag [default: ""]" << endl;
    cout << "   --mctype            MC campaign (allowed: mc15b, mc15c) [default: mc15c]" << endl;
    cout << "   --sys               run with systematics [default: false]" << endl;
    cout << "   --nLepFilter        filter on number of light leptons [default: 0]" << endl;
    cout << "   --trigFilter        require a trigger to have fired to store event [default: false]" << endl;

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
    bool run_sys = false;
    bool run_lep_filter = false;
    int n_lep_filter = 0;
    bool run_trig_filter = false;

    //cout << "NtMaker" << endl;

    int optin = 1;
    while(optin < argc) {
        string in = argv[optin];
        if      (in == "-f" || in == "--filelist"   ) { filelist = argv[++optin]; }
        else if (in == "-i" || in == "--input"      ) { input_container = argv[++optin]; }
        else if (in == "-o" || in == "--output"     ) { output_container = argv[++optin]; } 
        else if (in == "-d" || in == "--debug"      ) { dbg = atoi(argv[++optin]); }
        else if (in == "-n" || in == "--nevents"    ) { n_events = atoi(argv[++optin]); }
        else if (in == "-t" || in == "--tag"        ) { nt_tag = argv[++optin]; }
        else if (in == "--mctype"                   ) { mc_type_str = argv[++optin]; } 
        else if (in == "--sys"                      ) { run_sys = true; }
        else if (in == "--nLepFilter"               ) { run_lep_filter = true; n_lep_filter = atoi(argv[++optin]); }
        else if (in == "--trigFilter"               ) { run_trig_filter = true; }
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

    ///////////////////////////////////////////////////////////////////////////
    // things look ok, lets start
    cout << "-------------------------------------------------------------" << endl;
    cout << " NtMaker run conditions" << endl;
    cout << endl;
    cout << "   filelist                : " << filelist << endl;
    cout << "   input container         : " << input_container << endl;
    cout << "   output container        : " << output_container << endl;
    cout << "   debug level             : " << dbg << endl;
    cout << "   n events to process     : " << n_events << endl;
    cout << "   production tag          : " << nt_tag << endl;
    cout << "   MC type                 : " << mc_type_str << endl;
    cout << "   run systematics         : " << run_sys << endl;
    cout << "   run lepton filter       : " << run_lep_filter << endl;
    if(run_lep_filter)
    cout << "   n leptons to filter (>=): " << n_lep_filter << endl;
    cout << "   filter on trigger       : " << run_trig_filter << endl;
    cout << "-------------------------------------------------------------" << endl;

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
    susyAna->set_mctype(mc_type);

    // start up the looper
    if(n_events < 0) n_events = n_entries;
    cout << endl;
    cout << "Total entries in input chain       : " << n_entries << endl;
    cout << "Total number of entries to process : " << n_events << endl;
    chain->Process(susyAna, "", n_events);

    cout << endl;
    cout << "NtMaker    Job done" << endl;

    // check for expected output here?

    delete chain;
    return 0;

} // end
