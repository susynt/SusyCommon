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
  // else if( s== "EG_RESOLUTION_LASTRESOLUTIONVARIATION" )   r = EG_RESOLUTION_MATERIALCALO_DN;
  // else if( s== "EG_RESOLUTION_MATERIALCALO__1down" )       r = EG_RESOLUTION_MATERIALCALO_DN;
  // else if( s== "EG_RESOLUTION_MATERIALCALO__1up" )         r = EG_RESOLUTION_MATERIALCALO_UP;
  // else if( s== "EG_RESOLUTION_MATERIALCRYO__1down" )       r = EG_RESOLUTION_MATERIALCRYO_DN;
  // else if( s== "EG_RESOLUTION_MATERIALCRYO__1up" )         r = EG_RESOLUTION_MATERIALCRYO_UP;
  // else if( s== "EG_RESOLUTION_MATERIALGAP__1down" )        r = EG_RESOLUTION_MATERIALGAP_DN;
  // else if( s== "EG_RESOLUTION_MATERIALGAP__1up" )          r = EG_RESOLUTION_MATERIALGAP_UP;
  // else if( s== "EG_RESOLUTION_MATERIALID__1down" )         r = EG_RESOLUTION_MATERIALID_DN;
  // else if( s== "EG_RESOLUTION_MATERIALID__1up" )           r = EG_RESOLUTION_MATERIALID_UP;
  // else if( s== "EG_RESOLUTION_NOMINAL" )                   r = EG_RESOLUTION_NOMINAL;
  // else if( s== "EG_RESOLUTION_NONE" )                      r = EG_RESOLUTION_NONE;
  // else if( s== "EG_RESOLUTION_PILEUP__1down" )             r = EG_RESOLUTION_PILEUP_DN;
  // else if( s== "EG_RESOLUTION_PILEUP__1up" )               r = EG_RESOLUTION_PILEUP_UP;
  // else if( s== "EG_RESOLUTION_SAMPLINGTERM__1down" )       r = EG_RESOLUTION_SAMPLINGTERM_DN;
  // else if( s== "EG_RESOLUTION_SAMPLINGTERM__1up" )         r = EG_RESOLUTION_SAMPLINGTERM_UP;
  // else if( s== "EG_RESOLUTION_ZSMEARING__1down" )          r = EG_RESOLUTION_ZSMEARING_DN;
  // else if( s== "EG_RESOLUTION_ZSMEARING__1up" )            r = EG_RESOLUTION_ZSMEARING_UP;
  else if( s== "EG_SCALE_ALL__1down" )                     r = EG_SCALE_ALL_DN;
  else if( s== "EG_SCALE_ALL__1up" )                       r = EG_SCALE_ALL_UP;
  // else if( s== "EG_SCALE_G4__1down" )                      r = EG_SCALE_G4_DN;
  // else if( s== "EG_SCALE_G4__1up" )                        r = EG_SCALE_G4_UP;
  // else if( s== "EG_SCALE_L1GAIN__1down" )                  r = EG_SCALE_L1GAIN_DN;
  // else if( s== "EG_SCALE_L1GAIN__1up" )                    r = EG_SCALE_L1GAIN_UP;
  // else if( s== "EG_SCALE_L2GAIN__1down" )                  r = EG_SCALE_L2GAIN_DN;
  // else if( s== "EG_SCALE_L2GAIN__1up" )                    r = EG_SCALE_L2GAIN_UP;
  // else if( s== "EG_SCALE_LARCALIB__1down" )                r = EG_SCALE_LARCALIB_DN;
  // else if( s== "EG_SCALE_LARCALIB__1up" )                  r = EG_SCALE_LARCALIB_UP;
  // else if( s== "EG_SCALE_LARELECCALIB__1down" )            r = EG_SCALE_LARELECCALIB_DN;
  // else if( s== "EG_SCALE_LARELECCALIB__1up" )              r = EG_SCALE_LARELECCALIB_UP;
  // else if( s== "EG_SCALE_LARELECUNCONV__1down" )           r = EG_SCALE_LARELECUNCONV_DN;
  // else if( s== "EG_SCALE_LARELECUNCONV__1up" )             r = EG_SCALE_LARELECUNCONV_UP;
  // else if( s== "EG_SCALE_LARUNCONVCALIB__1down" )          r = EG_SCALE_LARUNCONVCALIB_DN;
  // else if( s== "EG_SCALE_LARUNCONVCALIB__1up" )            r = EG_SCALE_LARUNCONVCALIB_UP;
  // else if( s== "EG_SCALE_LASTSCALEVARIATION" )             r = EG_SCALE_LASTSCALEVARIATION;
  // else if( s== "EG_SCALE_MATCALO__1down" )                 r = EG_SCALE_MATCALO_DN;
  // else if( s== "EG_SCALE_MATCALO__1up" )                   r = EG_SCALE_MATCALO_UP;
  // else if( s== "EG_SCALE_MATCRYO__1down" )                 r = EG_SCALE_MATCRYO_DN;
  // else if( s== "EG_SCALE_MATCRYO__1up" )                   r = EG_SCALE_MATCRYO_UP;
  // else if( s== "EG_SCALE_MATID__1down" )                   r = EG_SCALE_MATID_DN;
  // else if( s== "EG_SCALE_MATID__1up" )                     r = EG_SCALE_MATID_UP;
  // else if( s== "EG_SCALE_NOMINAL" )                        r = EG_SCALE_NOMINAL;
  // else if( s== "EG_SCALE_NONE" )                           r = EG_SCALE_NONE;
  // else if( s== "EG_SCALE_PEDESTAL__1down" )                r = EG_SCALE_PEDESTAL_DN;
  // else if( s== "EG_SCALE_PEDESTAL__1up" )                  r = EG_SCALE_PEDESTAL_UP;
  // else if( s== "EG_SCALE_PS__1down" )                      r = EG_SCALE_PS_DN;
  // else if( s== "EG_SCALE_PS__1up" )                        r = EG_SCALE_PS_UP;
  // else if( s== "EG_SCALE_S12__1down" )                     r = EG_SCALE_S12_DN;
  // else if( s== "EG_SCALE_S12__1up" )                       r = EG_SCALE_S12_UP;
  // else if( s== "EG_SCALE_ZEESTAT__1down" )                 r = EG_SCALE_ZEESTAT_DN;
  // else if( s== "EG_SCALE_ZEESTAT__1up" )                   r = EG_SCALE_ZEESTAT_UP;
  // else if( s== "EG_SCALE_ZEESYST__1down" )                 r = EG_SCALE_ZEESYST_DN;
  // else if( s== "EG_SCALE_ZEESYST__1up" )                   r = EG_SCALE_ZEESYST_UP;
  
  else if( s== "EL_EFF_ID_TotalCorrUncertainty__1down" )   r = EL_EFF_ID_TotalCorrUncertainty_DN;
  else if( s== "EL_EFF_ID_TotalCorrUncertainty__1up" )     r = EL_EFF_ID_TotalCorrUncertainty_UP;
  else if( s== "EL_EFF_Iso_TotalCorrUncertainty__1down" )  r = EL_EFF_Iso_TotalCorrUncertainty_DN;
  else if( s== "EL_EFF_Iso_TotalCorrUncertainty__1up" )    r = EL_EFF_Iso_TotalCorrUncertainty_UP;
  else if( s== "EL_EFF_Reco_TotalCorrUncertainty__1down" ) r = EL_EFF_Reco_TotalCorrUncertainty_DN;
  else if( s== "EL_EFF_Reco_TotalCorrUncertainty__1up" )   r = EL_EFF_Reco_TotalCorrUncertainty_UP;
  else if( s== "EL_EFF_Trigger_TotalCorrUncertainty__1down" ) r = EL_EFF_Trigger_TotalCorrUncertainty_DN;
  else if( s== "EL_EFF_Trigger_TotalCorrUncertainty__1up" )   r = EL_EFF_Trigger_TotalCorrUncertainty_UP;

  // else if( s== "EL_SCALE_MOMENTUM__1down" )                r = EL_SCALE_MOMENTUM_DN;
  // else if( s== "EL_SCALE_MOMENTUM__1up" )                  r = EL_SCALE_MOMENTUM_UP;
  else if( s== "FT_EFF_B_systematics__1down" )             r = FT_EFF_B_systematics_DN;
  else if( s== "FT_EFF_B_systematics__1up" )               r = FT_EFF_B_systematics_UP;
  else if( s== "FT_EFF_C_systematics__1down" )             r = FT_EFF_C_systematics_DN;
  else if( s== "FT_EFF_C_systematics__1up" )               r = FT_EFF_C_systematics_UP;
  else if( s== "FT_EFF_Light_systematics__1down" )         r = FT_EFF_Light_systematics_DN;
  else if( s== "FT_EFF_Light_systematics__1up" )           r = FT_EFF_Light_systematics_UP;
  else if( s== "FT_EFF_extrapolation__1down" )             r = FT_EFF_extrapolation_DN;
  else if( s== "FT_EFF_extrapolation__1up" )               r = FT_EFF_extrapolation_UP;
  else if( s== "FT_EFF_extrapolation from charm__1down" )  r = FT_EFF_extrapolation_charm_DN;
  else if( s== "FT_EFF_extrapolation from charm__1up" )    r = FT_EFF_extrapolation_charm_UP;
/*
  else if( s== "FT_Eigen_B_0__1down" )                     r = FT_Eigen_B_0_DN;
  else if( s== "FT_Eigen_B_0__1up" )                       r = FT_Eigen_B_0_UP;
  else if( s== "FT_Eigen_B_1__1down" )                     r = FT_Eigen_B_1_DN;
  else if( s== "FT_Eigen_B_1__1up" )                       r = FT_Eigen_B_1_UP;
  else if( s== "FT_Eigen_B_2__1down" )                     r = FT_Eigen_B_2_DN;
  else if( s== "FT_Eigen_B_2__1up" )                       r = FT_Eigen_B_2_UP;
  else if( s== "FT_Eigen_B_3__1down" )                     r = FT_Eigen_B_3_DN;
  else if( s== "FT_Eigen_B_3__1up" )                       r = FT_Eigen_B_3_UP;
  else if( s== "FT_Eigen_B_4__1down" )                     r = FT_Eigen_B_4_DN;
  else if( s== "FT_Eigen_B_4__1up" )                       r = FT_Eigen_B_4_UP;
  else if( s== "FT_Eigen_B_5__1down" )                     r = FT_Eigen_B_5_DN;
  else if( s== "FT_Eigen_B_5__1up" )                       r = FT_Eigen_B_5_UP;
  else if( s== "FT_Eigen_B_6__1down" )                     r = FT_Eigen_B_6_DN;
  else if( s== "FT_Eigen_B_6__1up" )                       r = FT_Eigen_B_6_UP;
  else if( s== "FT_Eigen_B_7__1down" )                     r = FT_Eigen_B_7_DN;
  else if( s== "FT_Eigen_B_7__1up" )                       r = FT_Eigen_B_7_UP;
  else if( s== "FT_Eigen_B_8__1down" )                     r = FT_Eigen_B_8_DN;
  else if( s== "FT_Eigen_B_8__1up" )                       r = FT_Eigen_B_8_UP;
  else if( s== "FT_Eigen_B_9__1down" )                     r = FT_Eigen_B_9_DN;
  else if( s== "FT_Eigen_B_9__1up" )                       r = FT_Eigen_B_9_UP;
  else if( s== "FT_Eigen_C_0__1down" )                     r = FT_Eigen_C_0_DN;
  else if( s== "FT_Eigen_C_0__1up" )                       r = FT_Eigen_C_0_UP;
  else if( s== "FT_Eigen_C_1__1down" )                     r = FT_Eigen_C_1_DN;
  else if( s== "FT_Eigen_C_1__1up" )                       r = FT_Eigen_C_1_UP;
  else if( s== "FT_Eigen_C_2__1down" )                     r = FT_Eigen_C_2_DN;
  else if( s== "FT_Eigen_C_2__1up" )                       r = FT_Eigen_C_2_UP;
  else if( s== "FT_Eigen_C_3__1down" )                     r = FT_Eigen_C_3_DN;
  else if( s== "FT_Eigen_C_3__1up" )                       r = FT_Eigen_C_3_UP;
  else if( s== "FT_Eigen_Light_0__1down" )                 r = FT_Eigen_Light_0_DN;
  else if( s== "FT_Eigen_Light_0__1up" )                   r = FT_Eigen_Light_0_UP;
  else if( s== "FT_Eigen_Light_10__1down" )                r = FT_Eigen_Light_10_DN;
  else if( s== "FT_Eigen_Light_10__1up" )                  r = FT_Eigen_Light_10_UP;
  else if( s== "FT_Eigen_Light_11__1down" )                r = FT_Eigen_Light_11_DN;
  else if( s== "FT_Eigen_Light_11__1up" )                  r = FT_Eigen_Light_11_UP;
  else if( s== "FT_Eigen_Light_1__1down" )                 r = FT_Eigen_Light_1_DN;
  else if( s== "FT_Eigen_Light_1__1up" )                   r = FT_Eigen_Light_1_UP;
  else if( s== "FT_Eigen_Light_2__1down" )                 r = FT_Eigen_Light_2_DN;
  else if( s== "FT_Eigen_Light_2__1up" )                   r = FT_Eigen_Light_2_UP;
  else if( s== "FT_Eigen_Light_3__1down" )                 r = FT_Eigen_Light_3_DN;
  else if( s== "FT_Eigen_Light_3__1up" )                   r = FT_Eigen_Light_3_UP;
  else if( s== "FT_Eigen_Light_4__1down" )                 r = FT_Eigen_Light_4_DN;
  else if( s== "FT_Eigen_Light_4__1up" )                   r = FT_Eigen_Light_4_UP;
  else if( s== "FT_Eigen_Light_5__1down" )                 r = FT_Eigen_Light_5_DN;
  else if( s== "FT_Eigen_Light_5__1up" )                   r = FT_Eigen_Light_5_UP;
  else if( s== "FT_Eigen_Light_6__1down" )                 r = FT_Eigen_Light_6_DN;
  else if( s== "FT_Eigen_Light_6__1up" )                   r = FT_Eigen_Light_6_UP;
  else if( s== "FT_Eigen_Light_7__1down" )                 r = FT_Eigen_Light_7_DN;
  else if( s== "FT_Eigen_Light_7__1up" )                   r = FT_Eigen_Light_7_UP;
  else if( s== "FT_Eigen_Light_8__1down" )                 r = FT_Eigen_Light_8_DN;
  else if( s== "FT_Eigen_Light_8__1up" )                   r = FT_Eigen_Light_8_UP;
  else if( s== "FT_Eigen_Light_9__1down" )                 r = FT_Eigen_Light_9_DN;
  else if( s== "FT_Eigen_Light_9__1up" )                   r = FT_Eigen_Light_9_UP;
*/
  else if( s== "JET_JER_SINGLE_NP__1up" )                  r = JER;
  else if( s== "JET_GroupedNP_1__1up" )                    r = JET_GroupedNP_1_UP;
  else if( s== "JET_GroupedNP_1__1down" )                  r = JET_GroupedNP_1_DN;
  else if( s== "JET_GroupedNP_2__1up" )                    r = JET_GroupedNP_2_UP;
  else if( s== "JET_GroupedNP_2__1down" )                  r = JET_GroupedNP_2_DN;
  else if( s== "JET_GroupedNP_3__1up" )                    r = JET_GroupedNP_3_UP;
  else if( s== "JET_GroupedNP_3__1down" )                  r = JET_GroupedNP_3_DN;
/*
  else if( s== "JET_BJES_Response__1up" )                  r = JET_BJES_Response_UP;
  else if( s== "JET_BJES_Response__1down" )                r = JET_BJES_Response_DN;
  else if( s== "JET_EffectiveNP_1__1up" )                  r = JET_EffectiveNP_1_UP;
  else if( s== "JET_EffectiveNP_1__1down" )                r = JET_EffectiveNP_1_DN;
  else if( s== "JET_EffectiveNP_2__1up" )                  r = JET_EffectiveNP_2_UP;
  else if( s== "JET_EffectiveNP_2__1down" )                r = JET_EffectiveNP_2_DN;
  else if( s== "JET_EffectiveNP_3__1up" )                  r = JET_EffectiveNP_3_UP;
  else if( s== "JET_EffectiveNP_3__1down" )                r = JET_EffectiveNP_3_DN;
  else if( s== "JET_EffectiveNP_4__1up" )                  r = JET_EffectiveNP_4_UP;
  else if( s== "JET_EffectiveNP_4__1down" )                r = JET_EffectiveNP_4_DN;
  else if( s== "JET_EffectiveNP_5__1up" )                  r = JET_EffectiveNP_5_UP;
  else if( s== "JET_EffectiveNP_5__1down" )                r = JET_EffectiveNP_5_DN;
  else if( s== "JET_EffectiveNP_6restTerm__1up" )          r = JET_EffectiveNP_6restTerm_UP;
  else if( s== "JET_EffectiveNP_6restTerm__1down" )        r = JET_EffectiveNP_6restTerm_DN;
  else if( s== "JET_EtaIntercalibration_Modelling__1up" )  r = JET_EtaIntercalibration_Modelling_UP;
  else if( s== "JET_EtaIntercalibration_Modelling__1down" )r = JET_EtaIntercalibration_Modelling_DN;
  else if( s== "JET_EtaIntercalibration_TotalStat__1up" )  r = JET_EtaIntercalibration_TotalStat_UP;
  else if( s== "JET_EtaIntercalibration_TotalStat__1down" )r = JET_EtaIntercalibration_TotalStat_DN;
  else if( s== "JET_Flavor_Composition__1up" )             r = JET_Flavor_Composition_UP;
  else if( s== "JET_Flavor_Composition__1down" )           r = JET_Flavor_Composition_DN;
  else if( s== "JET_Flavor_Response__1up" )                r = JET_Flavor_Response_UP;
  else if( s== "JET_Flavor_Response__1down" )              r = JET_Flavor_Response_DN;
  else if( s== "JET_Pileup_OffsetMu__1up" )                r = JET_Pileup_OffsetMu_UP;
  else if( s== "JET_Pileup_OffsetMu__1down" )              r = JET_Pileup_OffsetMu_DN;
  else if( s== "JET_Pileup_OffsetNPV__1up" )               r = JET_Pileup_OffsetNPV_UP;
  else if( s== "JET_Pileup_OffsetNPV__1down" )             r = JET_Pileup_OffsetNPV_DN;
  else if( s== "JET_Pileup_PtTerm__1up" )                  r = JET_Pileup_PtTerm_UP;
  else if( s== "JET_Pileup_PtTerm__1down" )                r = JET_Pileup_PtTerm_DN;
  else if( s== "JET_Pileup_RhoTopology__1up" )             r = JET_Pileup_RhoTopology_UP;
  else if( s== "JET_Pileup_RhoTopology__1down" )           r = JET_Pileup_RhoTopology_DN;
  else if( s== "JET_PunchThrough_MC12__1up" )              r = JET_PunchThrough_MC12_UP;
  else if( s== "JET_PunchThrough_MC12__1down" )            r = JET_PunchThrough_MC12_DN;
  else if( s== "JET_SingleParticle_HighPt__1up" )          r = JET_SingleParticle_HighPt_UP;
  else if( s== "JET_SingleParticle_HighPt__1down" )        r = JET_SingleParticle_HighPt_DN;
  //else if( s== "JET_RelativeNonClosure_MC12__1up" )        r = JET_RelativeNonClosure_MC12_UP;
  //else if( s== "JET_RelativeNonClosure_MC12__1down" )      r = JET_RelativeNonClosure_MC12_DN;
*/
//NEW
  else if( s== "MET_SoftCalo_Reso" )                       r = MET_SoftCalo_Reso;
  else if( s== "MET_SoftCalo_ScaleDown" )                  r = MET_SoftCalo_ScaleDown;
  else if( s== "MET_SoftCalo_ScaleUp" )                    r = MET_SoftCalo_ScaleUp;
  else if( s== "MET_SoftTrk_ResoPara" )                    r = MET_SoftTrk_ResoPara;
  else if( s== "MET_SoftTrk_ResoPerp" )                    r = MET_SoftTrk_ResoPerp;
  else if( s== "MET_SoftTrk_ScaleDown" )                   r = MET_SoftTrk_ScaleDown;
  else if( s== "MET_SoftTrk_ScaleUp" )                     r = MET_SoftTrk_ScaleUp;
//NEW END
//  else if( s== "MUONSFSTAT__1down" )                       r = MUONSFSTAT_DN;
//  else if( s== "MUONSFSTAT__1up" )                         r = MUONSFSTAT_UP;
//  else if( s== "MUONSFSYS__1down" )                        r = MUONSFSYS_DN;
//  else if( s== "MUONSFSYS__1up" )                          r = MUONSFSYS_UP;
  else if( s== "MUON_EFF_STAT__1down" )                    r = MUONSFSTAT_DN;
  else if( s== "MUON_EFF_STAT__1up" )                      r = MUONSFSTAT_UP;
  else if( s== "MUON_EFF_SYS__1down" )                     r = MUONSFSYS_DN;
  else if( s== "MUON_EFF_SYS__1up" )                       r = MUONSFSYS_UP;
  else if( s== "MUONS_ID__1down" )                         r = MUONS_ID_DN;
  else if( s== "MUONS_ID__1up" )                           r = MUONS_ID_UP;
  else if( s== "MUONS_MS__1down" )                         r = MUONS_MS_DN;
  else if( s== "MUONS_MS__1up" )                           r = MUONS_MS_UP;
  else if( s== "MUONS_SCALE__1down" )                      r = MUONS_SCALE_DN;
  else if( s== "MUONS_SCALE__1up" )                        r = MUONS_SCALE_UP;
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
