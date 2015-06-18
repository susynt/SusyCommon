#ifndef SusyCommon_SystematicMapping_h
#define SusyCommon_SystematicMapping_h

/*
List of possible systematic uncertainties 
Perform the mapping between SusyNt list and the CP list.
*/

#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/SusyNtSys.h"
#include <string>

namespace Susy {

namespace NtSys {

 inline bool isValid(const SusyNtSys &s) { return s>=NOM && s<SYSUNKNOWN; }
 inline std::string syst2str(const SusyNtSys &s) { return isValid(s) ? SusyNtSysNames[s] : "unknown"; }
 
 SusyNtSys CPsys2sys(const std::string &s); //!< convert CP:Systematic string to our sys list
 

} //NtSys

} // Susy
#endif
