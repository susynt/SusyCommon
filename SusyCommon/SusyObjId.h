// Dear emacs, this is -*- c++ -*-
#ifndef SUSYCOMMON_SUSYOBJID_H
#define SUSYCOMMON_SUSYOBJID_H

#include <string>
#include <vector>

namespace Susy
{
enum SusyObjId {
    // the separate electron and muon WP
    // should be in the same order that corresponds
    // to SusyNtuple/ElectronId.h, SusyNtuple/MuonId.h!!!
    eleTightLLH=0
    ,eleMediumLLH
    ,eleLooseLLH
    ,muoLoose
    ,muoMedium
    ,Invalid
};

/// Human-readable names
std::string SusyObjId2str(const SusyObjId &id);

/// is an electron ID or not
bool isEleObj(const SusyObjId &id);

/// list of electron working points XaodAnalysis will loop over
/**
   Please use these function rather than using a range-based loop

   \code{.cpp}
   for(int i=SusyObjId::eleTightLLH; i<SusyObjId::Invalid; i++)
   \endcode
*/
inline std::vector<SusyObjId> electronIds() { return {eleTightLLH, eleMediumLLH }; }

/// list of muon working points XaodAnalysis will loop over
inline std::vector<SusyObjId> muonIds() { return { muoLoose, muoMedium }; }

/// list of electron+muon working points XaodAnalysis will loop over
inline std::vector<SusyObjId> leptonIds() { return { eleTightLLH, eleMediumLLH, muoLoose, muoMedium }; }

} // namespace Susy

#endif
