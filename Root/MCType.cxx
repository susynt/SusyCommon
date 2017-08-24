#include "SusyCommon/MCType.h"

using namespace std;

namespace Susy {

    string MCType2str(const MCType &type)
    {
        string s = "Invalid";
        switch(type) {
            case MCType::MC15b          : s = "MC15b"       ; break;
            case MCType::MC15c          : s = "MC15c"       ; break;
            case MCType::MCInvalid      : s = "Invalid"     ; break;
        } // switch 
        return s;
    }


} // namespace Susy
