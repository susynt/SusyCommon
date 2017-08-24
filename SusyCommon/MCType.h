#ifndef SUSYCOMMON_MCTYPE_H
#define SUSYCOMMON_MCTYPE_H

#include <string>

namespace Susy {

    enum MCType {
        MC15b=0,
        MC15c,
        MCInvalid
    };

    std::string MCType2str(const MCType &type);

} // namespace Susy


#endif
