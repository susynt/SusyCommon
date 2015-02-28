#ifndef SusyCommon_TriggerMap_h
#define SusyCommon_TriggerMap_h

#include <iostream>
#include <string>


#include "TBits.h"

using namespace std;

// --------------------------------------------------
// TriggerMap
//
// A little module to hold anything related to the 
// trigger bits and flags.
// --------------------------------------------------
/*
struct triggerbits
{
    enum bits {
        // 2012 triggers
        e7_medium1 = 0,
        HLT,
        e12Tvh_loose1,
        e12Tvh_medium1,
     //   e24vh_medium1,
     //   e24vhi_medium1,
     //   2e12Tvh_loose1,
     //   e24vh_medium1_e7_medium1,
    
        N_TRIG
    };

};
*/
namespace triggerbits {

    enum {
        // 2012 triggers
        BIT_e7_medium1 = 0,
        BIT_e12Tvh_loose1,
        BIT_e12Tvh_medium1,
        BIT_e24vh_medium1,
        BIT_e24vhi_medium1,
        BIT_2e12Tvh_loose1,
        BIT_e24vh_medium1_e7_medium1,
        
        BIT_mu8,
        BIT_mu13,
        BIT_mu18_tight,
        BIT_mu24i_tight,
        BIT_2mu13,
        BIT_mu18_tight_mu8_EFFS,
        BIT_e12Tvh_medium1_mu8,
        BIT_mu18_tight_e7_medium1,

        // photon triggers
        BIT_g20_loose,
        BIT_g40_loose,
        BIT_g60_loose,
        BIT_g80_loose,
        BIT_g100_loose,
        BIT_g120_loose,

        // Tau triggers
        BIT_tau20_medium1,
        BIT_tau20Ti_medium1,
        BIT_tau29Ti_medium1,
        BIT_tau29Ti_medium1_tau20Ti_medium1,
        BIT_tau20Ti_medium1_e18vh_medium1,
        BIT_tau20_medium1_mu15,

        // trigger flags for lep-tau  matching
        BIT_e18vh_medium1,
        BIT_mu15,

        // MET triggers
        BIT_2mu8_EFxe40wMu_tclcw,
        BIT_xe80_tclcw_loose,
        BIT_xe80T_tclcw_loose,
        
        
        
 
        N_TRIG
    };
    
    // 2012 trigger bit masks
    const long long TRIG_e7_medium1     = 1LL << BIT_e7_medium1;        
    const long long TRIG_e12Tvh_loose1 = 1LL << BIT_e12Tvh_loose1;      
    const long long TRIG_e12Tvh_medium1 = 1LL << BIT_e12Tvh_medium1;    
    const long long TRIG_e24vh_medium1  = 1LL << BIT_e24vh_medium1;
    const long long TRIG_e24vhi_medium1 = 1LL << BIT_e24vhi_medium1;
    const long long TRIG_2e12Tvh_loose1 = 1LL << BIT_2e12Tvh_loose1;
    const long long TRIG_e24vh_medium1_e7_medium1 = 1LL << BIT_e24vh_medium1_e7_medium1;

    const long long TRIG_mu8            = 1LL << BIT_mu8;
    const long long TRIG_mu13           = 1LL << BIT_mu13;
    const long long TRIG_mu18_tight     = 1LL << BIT_mu18_tight;
    const long long TRIG_mu24i_tight    = 1LL << BIT_mu24i_tight;
    const long long TRIG_2mu13          = 1LL << BIT_2mu13;
    const long long TRIG_mu18_tight_mu8_EFFS = 1LL << BIT_mu18_tight_mu8_EFFS;
    
    const long long TRIG_e12Tvh_medium1_mu8 = 1LL << BIT_e12Tvh_medium1_mu8;
    const long long TRIG_mu18_tight_e7_medium1 = 1LL << BIT_mu18_tight_e7_medium1;

    // photon trigger bit masks
    const long long TRIG_g20_loose      = 1LL << BIT_g20_loose;
    const long long TRIG_g40_loose      = 1LL << BIT_g40_loose;
    const long long TRIG_g60_loose      = 1LL << BIT_g60_loose;
    const long long TRIG_g80_loose      = 1LL << BIT_g80_loose;
    const long long TRIG_g100_loose     = 1LL << BIT_g100_loose;
    const long long TRIG_g120_loose     = 1LL << BIT_g120_loose;

    // tau trigger bit masks
    const long long TRIG_tau20_medium1  = 1LL << BIT_tau20_medium1;
    const long long TRIG_tau20Ti_medium1 = 1LL << BIT_tau20Ti_medium1;
    const long long TRIG_tau29Ti_medium1 = 1LL << BIT_tau29Ti_medium1;
    const long long TRIG_tau29Ti_medium1_tau20Ti_medium1 = 1LL << BIT_tau29Ti_medium1_tau20Ti_medium1;
    const long long TRIG_tau20Ti_medium1_e18vh_medium1 = 1LL << BIT_tau20Ti_medium1_e18vh_medium1;
    const long long TRIG_tau20_medium1_mu15 = 1LL << BIT_tau20_medium1_mu15;

    // lep-tau matching
    const long long TRIG_e18vh_medium1  = 1LL << BIT_e18vh_medium1;
    const long long TRIG_mu15           = 1LL << BIT_mu15;

    // MET triggers
    const long long TRIG_2mu8_EFxe40wMu_tclcw = 1LL << BIT_2mu8_EFxe40wMu_tclcw;
    const long long TRIG_xe80_tclcw_loose = 1LL << BIT_xe80_tclcw_loose;
    const long long TRIG_xe80T_tclcw_loose = 1LL << BIT_xe80T_tclcw_loose;
 
    
    const string trigger_names[N_TRIG] = {
        "EF_e7_medium1",
        "EF_e12Tvh_loose1",
        "EF_e12Tvh_medium1"
        "EF_e24vh_medium1",
        "EF_e24vhi_medium1",
        "EF_e24vh_medium1_e7_medium1",
        
        "EF_mu8",
        "EF_mu13",
        "EF_mu18_tight",
        "EF_mu24i_tight",
        "EF_2mu13",
        "EF_mu18_tight_mu8_EFFS",
        "EF_e12Tvh_medium1_mu8",
        "EF_mu18_tight_e7_medium1",

        // photon triggers
        "EF_g20_loose",
        "EF_g40_loose",
        "EF_g60_loose",
        "EF_g80_loose",
        "EF_g100_loose",
        "EF_g120_loose",
        
        // tau triggers
        "EF_tau20_medium1",
        "EF_tau20Ti_medium1",
        "EF_tau29Ti_medium1",
        "EF_tau29Ti_medium1_tau20Ti_medium1",
        "EF_tau20Ti_medium1_e18vh_medium1",
        "EF_tau20_medium1_mu15",

        // lep-tau matching
        "EF_e18vh_medium1",
        "EF_mu15",
        
        // MET triggers
        "EF_2mu8_EFxe40wMu_tclcw",
        "EF_xe80_tclcw_loose",
        "EF_xe80T_tclcw_loose"
      
    };

}; // end namespace triggerbits

#endif
