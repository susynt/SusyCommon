#ifndef SusyCommon_TriggerMap_h
#define SusyCommon_TriggerMap_h

#include <string>


// --------------------------------------------------
// TriggerMap
//
// A little home for the triggers used in SusyNt
// 
// --------------------------------------------------


namespace susy {

    const std::vector<std::string> triggerNames = {
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
        


}; // end namespace susy

#endif
