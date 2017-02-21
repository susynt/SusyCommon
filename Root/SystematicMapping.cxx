#include "SusyCommon/SystematicMapping.h"

#include <cassert>

namespace Susy {

namespace NtSys {

//-----------------------------------------
SusyNtSys CPsys2sys(const std::string &s)
{
  SusyNtSys r = SYS_UNKNOWN;
  
  
  if     ( s== "EG_RESOLUTION_ALL__1down" )                r = EG_RESOLUTION_ALL_DN;
  else if( s== "EG_RESOLUTION_ALL__1up" )                  r = EG_RESOLUTION_ALL_UP;
  else if( s== "EG_SCALE_ALL__1down" )                     r = EG_SCALE_ALL_DN;
  else if( s== "EG_SCALE_ALL__1up" )                       r = EG_SCALE_ALL_UP;

  else if( s== "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1down" )         r = EL_EFF_ID_TOTAL_Uncorr_DN;
  else if( s== "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1up" )           r = EL_EFF_ID_TOTAL_Uncorr_UP; 
  else if( s== "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1down" )        r = EL_EFF_Iso_TOTAL_Uncorr_DN; 
  else if( s== "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR__1up" )          r = EL_EFF_Iso_TOTAL_Uncorr_UP; 
  else if( s== "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1down" )       r = EL_EFF_Reco_TOTAL_Uncorr_DN; 
  else if( s== "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1up" )         r = EL_EFF_Reco_TOTAL_Uncorr_UP; 
  else if( s== "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR__1down")  r = EL_EFF_Trigger_TOTAL_Uncorr_DN; 
  else if( s== "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR__1up")    r = EL_EFF_Trigger_TOTAL_Uncorr_UP; 

  else if( s== "FT_EFF_B_systematics__1down" )             r = FT_EFF_B_systematics_DN;
  else if( s== "FT_EFF_B_systematics__1up" )               r = FT_EFF_B_systematics_UP;
  else if( s== "FT_EFF_C_systematics__1down" )             r = FT_EFF_C_systematics_DN;
  else if( s== "FT_EFF_C_systematics__1up" )               r = FT_EFF_C_systematics_UP;
  else if( s== "FT_EFF_Light_systematics__1down" )         r = FT_EFF_Light_systematics_DN;
  else if( s== "FT_EFF_Light_systematics__1up" )           r = FT_EFF_Light_systematics_UP;
  else if( s== "FT_EFF_extrapolation__1down" )             r = FT_EFF_extrapolation_DN;
  else if( s== "FT_EFF_extrapolation__1up" )               r = FT_EFF_extrapolation_UP;
  else if( s== "FT_EFF_extrapolation_from_charm__1down" )  r = FT_EFF_extrapolation_charm_DN;
  else if( s== "FT_EFF_extrapolation_from_charm__1up" )    r = FT_EFF_extrapolation_charm_UP;

  else if( s== "JET_JER_SINGLE_NP__1up" )                  r = JER;
  else if( s== "JET_GroupedNP_1__1up" )                    r = JET_GroupedNP_1_UP;
  else if( s== "JET_GroupedNP_1__1down" )                  r = JET_GroupedNP_1_DN;
  else if( s== "JET_GroupedNP_2__1up" )                    r = JET_GroupedNP_2_UP;
  else if( s== "JET_GroupedNP_2__1down" )                  r = JET_GroupedNP_2_DN;
  else if( s== "JET_GroupedNP_3__1up" )                    r = JET_GroupedNP_3_UP;
  else if( s== "JET_GroupedNP_3__1down" )                  r = JET_GroupedNP_3_DN;
  else if( s== "JET_JvtEfficiency__1up" )                  r = JET_JVTEff_UP;
  else if( s== "JET_JvtEfficiency__1down" )                r = JET_JVTEff_DN;

  else if( s== "MET_SoftCalo_Reso" )                       r = MET_SoftCalo_Reso;
  else if( s== "MET_SoftCalo_ScaleDown" )                  r = MET_SoftCalo_ScaleDown;
  else if( s== "MET_SoftCalo_ScaleUp" )                    r = MET_SoftCalo_ScaleUp;
  else if( s== "MET_SoftTrk_ResoPara" )                    r = MET_SoftTrk_ResoPara;
  else if( s== "MET_SoftTrk_ResoPerp" )                    r = MET_SoftTrk_ResoPerp;
  else if( s== "MET_SoftTrk_ScaleDown" )                   r = MET_SoftTrk_ScaleDown;
  else if( s== "MET_SoftTrk_ScaleUp" )                     r = MET_SoftTrk_ScaleUp;

  else if( s== "MUON_EFF_STAT__1down" )                    r = MUON_EFF_STAT_DN;
  else if( s== "MUON_EFF_STAT__1up" )                      r = MUON_EFF_STAT_UP;
  else if( s== "MUON_EFF_STAT_LOWPT__1down" )              r = MUON_EFF_STAT_LOWPT_DN;
  else if( s== "MUON_EFF_STAT_LOWPT__1up" )                r = MUON_EFF_STAT_LOWPT_UP;
  else if( s== "MUON_EFF_SYS__1down" )                     r = MUON_EFF_SYS_DN;
  else if( s== "MUON_EFF_SYS__1up" )                       r = MUON_EFF_SYS_UP;
  else if( s== "MUON_EFF_SYS_LOWPT__1down" )               r = MUON_EFF_SYS_LOWPT_DN;
  else if( s== "MUON_EFF_SYS_LOWPT__1up" )                 r = MUON_EFF_SYS_LOWPT_UP;
  else if( s== "MUON_EFF_TrigStatUncertainty__1down" )     r = MUON_EFF_TRIG_STAT_DN;
  else if( s== "MUON_EFF_TrigStatUncertainty__1up" )       r = MUON_EFF_TRIG_STAT_UP;
  else if( s== "MUON_EFF_TrigSystUncertainty__1down" )     r = MUON_EFF_TRIG_SYST_DN;
  else if( s== "MUON_EFF_TrigSystUncertainty__1up" )       r = MUON_EFF_TRIG_SYST_UP;
  else if( s== "MUON_ISO_STAT__1down" )                    r = MUON_ISO_STAT_DN;
  else if( s== "MUON_ISO_STAT__1up" )                      r = MUON_ISO_STAT_UP;
  else if( s== "MUON_ISO_SYS__1down" )                     r = MUON_ISO_SYS_DN;
  else if( s== "MUON_ISO_SYS__1up" )                       r = MUON_ISO_SYS_UP;
  else if( s== "MUON_ID__1down" )                          r = MUON_ID_DN;
  else if( s== "MUON_ID__1up" )                            r = MUON_ID_UP;
  else if( s== "MUON_MS__1down" )                          r = MUON_MS_DN;
  else if( s== "MUON_MS__1up" )                            r = MUON_MS_UP;
  else if( s== "MUON_SCALE__1down" )                       r = MUON_SCALE_DN;
  else if( s== "MUON_SCALE__1up" )                         r = MUON_SCALE_UP;
  else if( s== "MUON_TTVA_STAT__1down" )                   r = MUON_TTVA_STAT_DN;
  else if( s== "MUON_TTVA_STAT__1up" )                     r = MUON_TTVA_STAT_UP;
  else if( s== "MUON_TTVA_SYS__1down" )                    r = MUON_TTVA_SYS_DN;
  else if( s== "MUON_TTVA_SYS__1up" )                      r = MUON_TTVA_SYS_UP; 
/*
  else if( s== "PH_SCALE_CONVEFFICIENCY__1down" )          r = PH_SCALE_CONVEFFICIENCY_DN;
  else if( s== "PH_SCALE_CONVEFFICIENCY__1up" )            r = PH_SCALE_CONVEFFICIENCY_UP;
  else if( s== "PH_SCALE_CONVFAKERATE__1down" )            r = PH_SCALE_CONVFAKERATE_DN;
  else if( s== "PH_SCALE_CONVFAKERATE__1up" )              r = PH_SCALE_CONVFAKERATE_UP;
  else if( s== "PH_SCALE_CONVRADIUS__1down" )              r = PH_SCALE_CONVRADIUS_DN;
  else if( s== "PH_SCALE_CONVRADIUS__1up" )                r = PH_SCALE_CONVRADIUS_UP;
  else if( s== "PH_SCALE_LEAKAGECONV__1down" )             r = PH_SCALE_LEAKAGECONV_DN;
  else if( s== "PH_SCALE_LEAKAGECONV__1up" )               r = PH_SCALE_LEAKAGECONV_UP;
  else if( s== "PH_SCALE_LEAKAGEUNCONV__1down" )           r = PH_SCALE_LEAKAGEUNCONV_DN;
  else if( s== "PH_SCALE_LEAKAGEUNCONV__1up" )             r = PH_SCALE_LEAKAGEUNCONV_UP;
*/
  else if( s== "TAUS_EFF_CONTJETID_STAT__1down" )          r = TAUS_EFF_CONTJETID_STAT_DN;
  else if( s== "TAUS_EFF_CONTJETID_STAT__1up" )            r = TAUS_EFF_CONTJETID_STAT_UP;
  else if( s== "TAUS_EFF_CONTJETID_SYST__1down" )          r = TAUS_EFF_CONTJETID_SYST_DN;
  else if( s== "TAUS_EFF_CONTJETID_SYST__1up" )            r = TAUS_EFF_CONTJETID_SYST_UP;
  else if( s== "TAUS_SME_TOTAL__1down" )                   r = TAUS_SME_TOTAL_DN;
  else if( s== "TAUS_SME_TOTAL__1up" )                     r = TAUS_SME_TOTAL_UP;
  // PileupReweighting (mu) variations
  else if( s== "PRW_DATASF__1down" )                       r = PILEUP_DN;
  else if( s== "PRW_DATASF__1up" )                         r = PILEUP_UP;
  
  return r;
}



} //Ntsys
} // susy
