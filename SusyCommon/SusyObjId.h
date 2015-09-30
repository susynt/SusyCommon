#ifndef SUSYCOMMON_SUSYOBJID_H
#define SUSYCOMMON_SUSYOBJID_H

#include <string>

namespace Susy
{
enum SusyObjId {
    // the separate electron and muon WP
    // should be in the same order that corresponds
    // to SusyNtuple/ElectronId.h, SusyNtuple/MuonId.h!!!
    eleTightLH=0
    ,eleMediumLH
    ,eleLooseLH
    ,muoLoose
    ,muoMedium
    ,Invalid
};

/// Human-readable names
std::string SusyObjId2str(const SusyObjId &id);

/// is an electron ID or not
bool isEleObj(const SusyObjId &id);

} // namespace Susy

#endif
