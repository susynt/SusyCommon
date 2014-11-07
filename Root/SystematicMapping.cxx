#include "SusyCommon/SystematicMapping.h"

#include <cassert>

namespace Susy {

namespace NtSys {

//-----------------------------------------
SusyNtSys CPsys2sys(const std::string &s)
{
  SusyNtSys r = SYSUNKNOWN;
  
  if     ( s== "EL_SCALE_MOMENTUM__1down" )               r = EL_SCALE_MOMENTUM_DOWN;
  else if( s== "EL_SCALE_MOMENTUM__1up" )                 r = EL_SCALE_MOMENTUM_UP;
  else if( s== "EG_RESOLUTION_ALL__1down" )               r = EG_RESOLUTION_ALL_DOWN;
  else if( s== "EG_RESOLUTION_ALL__1up" )                 r = EG_RESOLUTION_ALL_UP;
  else if( s== "EG_RESOLUTION_LASTRESOLUTIONVARIATION" )  r = EG_RESOLUTION_LASTRESOLUTIONVARIATION;
  else if( s== "EG_RESOLUTION_MATERIALCALO__1down" )      r = EG_RESOLUTION_MATERIALCALO_DOWN;
  else if( s== "EG_RESOLUTION_MATERIALCALO__1up" )        r = EG_RESOLUTION_MATERIALCALO_UP;
  else if( s== "EG_RESOLUTION_MATERIALCRYO__1down" )      r = EG_RESOLUTION_MATERIALCRYO_DOWN;
  else if( s== "EG_RESOLUTION_MATERIALCRYO__1up" )        r = EG_RESOLUTION_MATERIALCRYO_UP;
  else if( s== "EG_RESOLUTION_MATERIALGAP__1down" )       r = EG_RESOLUTION_MATERIALGAP_DOWN;
  else if( s== "EG_RESOLUTION_MATERIALGAP__1up" )         r = EG_RESOLUTION_MATERIALGAP_UP;
  else if( s== "EG_RESOLUTION_MATERIALID__1down" )        r = EG_RESOLUTION_MATERIALID_DOWN;
  else if( s== "EG_RESOLUTION_MATERIALID__1up" )          r = EG_RESOLUTION_MATERIALID_UP;
  else if( s== "EG_RESOLUTION_NOMINAL" )                  r = EG_RESOLUTION_NOMINAL;
  else if( s== "EG_RESOLUTION_NONE" )                     r = EG_RESOLUTION_NONE;
  else if( s== "EG_RESOLUTION_PILEUP__1down" )            r = EG_RESOLUTION_PILEUP_DOWN;
  else if( s== "EG_RESOLUTION_PILEUP__1up" )              r = EG_RESOLUTION_PILEUP_UP;
  else if( s== "EG_RESOLUTION_SAMPLINGTERM__1down" )      r = EG_RESOLUTION_SAMPLINGTERM_DOWN;
  else if( s== "EG_RESOLUTION_SAMPLINGTERM__1up" )        r = EG_RESOLUTION_SAMPLINGTERM_UP;
  else if( s== "EG_RESOLUTION_ZSMEARING__1down" )         r = EG_RESOLUTION_ZSMEARING_DOWN;
  else if( s== "EG_RESOLUTION_ZSMEARING__1up" )           r = EG_RESOLUTION_ZSMEARING_UP;
  else if( s== "EG_SCALE_ALL__1down" )                    r = EG_SCALE_ALL_DOWN;
  else if( s== "EG_SCALE_ALL__1up" )                      r = EG_SCALE_ALL_UP;
  else if( s== "EG_SCALE_G4__1down" )                     r = EG_SCALE_G4_DOWN;
  else if( s== "EG_SCALE_G4__1up" )                       r = EG_SCALE_G4_UP;
  else if( s== "EG_SCALE_L1GAIN__1down" )                 r = EG_SCALE_L1GAIN_DOWN;
  else if( s== "EG_SCALE_L1GAIN__1up" )                   r = EG_SCALE_L1GAIN_UP;
  else if( s== "EG_SCALE_L2GAIN__1down" )                 r = EG_SCALE_L2GAIN_DOWN;
  else if( s== "EG_SCALE_L2GAIN__1up" )                   r = EG_SCALE_L2GAIN_UP;
  else if( s== "EG_SCALE_LARCALIB__1down" )               r = EG_SCALE_LARCALIB_DOWN;
  else if( s== "EG_SCALE_LARCALIB__1up" )                 r = EG_SCALE_LARCALIB_UP;
  else if( s== "EG_SCALE_LARELECCALIB__1down" )           r = EG_SCALE_LARELECCALIB_DOWN;
  else if( s== "EG_SCALE_LARELECCALIB__1up" )             r = EG_SCALE_LARELECCALIB_UP;
  else if( s== "EG_SCALE_LARELECUNCONV__1down" )          r = EG_SCALE_LARELECUNCONV_DOWN;
  else if( s== "EG_SCALE_LARELECUNCONV__1up" )            r = EG_SCALE_LARELECUNCONV_UP;
  else if( s== "EG_SCALE_LARUNCONVCALIB__1down" )         r = EG_SCALE_LARUNCONVCALIB_DOWN;
  else if( s== "EG_SCALE_LARUNCONVCALIB__1up" )           r = EG_SCALE_LARUNCONVCALIB_UP;
  else if( s== "EG_SCALE_LASTSCALEVARIATION ")            r = EG_SCALE_LASTSCALEVARIATION;
  else if( s== "EG_SCALE_MATCALO__1down" )                r = EG_SCALE_MATCALO_DOWN;
  else if( s== "EG_SCALE_MATCALO__1up" )                  r = EG_SCALE_MATCALO_UP;
  else if( s== "EG_SCALE_MATCRYO__1down" )                r = EG_SCALE_MATCRYO_DOWN;
  else if( s== "EG_SCALE_MATCRYO__1up" )                  r = EG_SCALE_MATCRYO_UP;
  else if( s== "EG_SCALE_MATID__1down" )                  r = EG_SCALE_MATID_DOWN;
  else if( s== "EG_SCALE_MATID__1up" )                    r = EG_SCALE_MATID_UP;
  else if( s== "EG_SCALE_NOMINAL")                        r = EG_SCALE_NOMINAL;
  else if( s== "EG_SCALE_NONE ")                          r = EG_SCALE_NONE;
  else if( s== "EG_SCALE_PEDESTAL__1down" )               r = EG_SCALE_PEDESTAL_DOWN;
  else if( s== "EG_SCALE_PEDESTAL__1up" )                 r = EG_SCALE_PEDESTAL_UP;
  else if( s== "EG_SCALE_PS__1down" )                     r = EG_SCALE_PS_DOWN;
  else if( s== "EG_SCALE_PS__1up" )                       r = EG_SCALE_PS_UP;
  else if( s== "EG_SCALE_S12__1down" )                    r = EG_SCALE_S12_DOWN;
  else if( s== "EG_SCALE_S12__1up" )                      r = EG_SCALE_S12_UP;
  else if( s== "EG_SCALE_ZEESTAT__1down" )                r = EG_SCALE_ZEESTAT_DOWN;
  else if( s== "EG_SCALE_ZEESTAT__1up" )                  r = EG_SCALE_ZEESTAT_UP;
  else if( s== "EG_SCALE_ZEESYST__1down" )                r = EG_SCALE_ZEESYST_DOWN;
  else if( s== "EG_SCALE_ZEESYST__1up" )                  r = EG_SCALE_ZEESYST_UP;

  else if( s== "MUONSFSTAT__1down" )                      r = MU_SF_STAT_DOWN;
  else if( s== "MUONSFSTAT__1up"   )                      r = MU_SF_STAT_UP;
  else if( s== "MUONSFSYS__1down"  )                      r = MU_SF_SYS_DOWN;
  else if( s== "MUONSFSYS__1up"    )                      r = MU_SF_SYS_UP;
  else if( s== "MUONS_ID__1down"   )                      r = MU_ID_DOWN;
  else if( s== "MUONS_ID__1up"     )                      r = MU_ID_UP;
  else if( s== "MUONS_MS__1down"   )                      r = MU_MS_DOWN;
  else if( s== "MUONS_MS__1up"     )                      r = MU_MS_UP;
  else if( s== "MUONS_SCALE__1down")                      r = MU_SCALE_DOWN;
  else if( s== "MUONS_SCALE__1up"  )                      r = MU_SCALE_UP;

  else if( s== "TAUS_STAT__1down")                        r = TAUS_STAT_DOWN;
  else if( s== "TAUS_STAT__1up")                          r = TAUS_STAT_UP; 
  else if( s== "TAUS_SYST__2down")                        r = TAUS_SYST_DOWN;
  else if( s== "TAUS_SYST__2up")                          r = TAUS_SYST_UP;
  else if( s== "TAUS_TOTAL__21down")                      r= TAUS_TOTAL_DOWN;
  else if( s== "TAUS_TOTAL__21up")                        r= TAUS_TOTAL_UP;

  else if( s== "JER__1up")                                r = JER_UP;
  else if( s== "BJES_Response__1down")                    r = BJES_Response_DOWN;
  else if( s== "BJES_Response__1up")                      r = BJES_Response_UP;

  return r;
}

//-----------------------------------------
bool isObjectSystematic(const SusyNtSys &s){
  bool b = false;
  switch(s) {
  case MU_ID_DOWN    : b=true; break;
  case MU_ID_UP      : b=true; break;
  case MU_MS_DOWN    : b=true; break;
  case MU_MS_UP      : b=true; break;
  case MU_SCALE_DOWN : b=true; break;
  case MU_SCALE_UP   : b=true; break;
  }
  return b;
}

} //Ntsys
} // susy
