#ifndef SusyCommon_TriggerMap_h
#define SusyCommon_TriggerMap_h

#include <string>


// --------------------------------------------------
// TriggerMap
//
// A little home for the triggers used in SusyNt
// 
//  TODO: dantrim Mar 17 2015 :
//     (1) implement a function that takes arguments
//          based (potentially) on metadata stored
//          in the xAOD itself or from user input 
//          and returns the appropriate vector
//          of triggers
//     (2)  implement a function that will be used
//          when reading SusyNt that will read the
//          trigger histogram to determine which
//          triggers are stored and from bin-labels
//          will then construct the appropriate
//          map between the bits actually stored
//          and the std::string of the trigger. User
//          can then call triggers by name, regardless
//          of how bits may have shifted.
//
// --------------------------------------------------


namespace susy {
//
// a pseudo-random set of triggers used for testing the code
// trigger naming found in :
// https://twiki.cern.ch/twiki/bin/view/Atlas/TriggerMenuDC14Run2#Main_physics_triggers_AN2
//
    const std::vector<std::string> triggerNames = {
        // electron triggers
                "HLT_e24_medium1_iloose",
                "HLT_e24_loose1",
                "HLT_e28_tight1_iloose",
                "HLT_e60_medium1",
                "HLT_e60_loose1",

        // muon triggers
                "HLT_mu26_imedium",
                "HLT_mu50",
                "HLT_mu60_0eta105_msonly",
                "HLT_2mu4",
                "HLT_2mu6",
                "HLT_2mu10",
                "HLT_2mu14",
                "HLT_3mu6",
                
        // photon triggers
                "HLT_g120_loose1",
                "HLT_g140_loose1",

        // tau triggers
                "HLT_e18_loose1_tau25_medium1_calo",
                "HLT_e18_lhloose1_tau25_medium1_calo",
                "HLT_mu14_tau25_medium1_calo",
                "HLT_tau35_medium1_calo_tau25_medium1_calo",

        // met triggers
                "HLT_xe100"

        };
                
                

/*
    const std::vector<std::string> triggerNames = {
       // electron triggers 
               "_e7_medium1",                       
                "_e12Tvh_loose1",                    
                "_e12Tvh_medium1",                   
                "_e24vh_medium1",                    
                "_e24vhi_medium1",                   
                "_e24vh_medium1_e7_medium1",      

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
        
*/

}; // end namespace susy

#endif
