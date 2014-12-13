#include "SusyCommon/SystematicMapping.h"

#include <cassert>

namespace Susy {

namespace NtSys {

//-----------------------------------------
SusyNtSys CPsys2sys(const std::string &s)
{
  SusyNtSys r = SYSUNKNOWN;
  
  
  if     ( s== "EG_RESOLUTION_ALL__1down" )               r = EG_RESOLUTION_ALL_DN;
  else if( s== "EG_RESOLUTION_ALL__1up" )                 r = EG_RESOLUTION_ALL_UP;


  return r;
}



} //Ntsys
} // susy
