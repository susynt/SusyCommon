#ifndef SusyCommon_TriggerMap_h
#define SusyCommon_TriggerMap_h

#include <iostream>
#include <string>
#include <map>


#include "TBits.h"

using namespace std;

// --------------------------------------------------
// TriggerMap
//
// A little module to hold anything related to the 
// trigger bits and flags.
// --------------------------------------------------


namespace TriggerMap {

    const std::vector<std::string> triggermap = {
       // electron triggers 
                "EF_e7_medium1",                       
                "EF_e12Tvh_loose1",                    
                "EF_e12Tvh_medium1",                   
                "EF_e24vh_medium1",                    
                "EF_e24vhi_medium1",                   
                "EF_e24vh_medium1_e7_medium1",         

        // muon triggers
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
        

/*
        const std::map<std::string, int> triggermap = {
       // electron triggers 
                {"EF_e7_medium1",                      0},            
                {"EF_e12Tvh_loose1",                   1},
                {"EF_e12Tvh_medium1",                  2},
                {"EF_e24vh_medium1",                   3}, 
                {"EF_e24vhi_medium1",                  4}, 
                {"EF_e24vh_medium1_e7_medium1",        5}, 

        // muon triggers
                {"EF_mu8",                              6},
                {"EF_mu13",                             7},
                {"EF_mu18_tight",                       8},
                {"EF_mu24i_tight",                      9},
                {"EF_2mu13",                            10}, 
                {"EF_mu18_tight_mu8_EFFS",              11},
                {"EF_e12Tvh_medium1_mu8",               12},
                {"EF_mu18_tight_e7_medium1",            13},

        // photon triggers
                {"EF_g20_loose",                        14},
                {"EF_g40_loose",                        15},
                {"EF_g60_loose",                        16},
                {"EF_g80_loose",                        17},
                {"EF_g100_loose",                       18},
                {"EF_g120_loose",                       19},

        // tau triggers
                {"EF_tau20_medium1",                    20},
                {"EF_tau20Ti_medium1",                  21},
                {"EF_tau29Ti_medium1",                  22},
                {"EF_tau29Ti_medium1_tau20Ti_medium1",  23},
                {"EF_tau20Ti_medium1_e18vh_medium1",    24},
                {"EF_tau20_medium1_mu15",               25},

        // lep-tau matching
                {"EF_e18vh_medium1",                    26},
                {"EF_mu15",                             27},
                       
        // MET triggers
                {"EF_2mu8_EFxe40wMu_tclcw",             28},
                {"EF_xe80_tclcw_loose",                 29},
                {"EF_xe80T_tclcw_loose"                 30}
        };
*/

}; // end namespace TriggerMap

#endif
