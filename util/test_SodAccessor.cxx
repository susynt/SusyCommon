
#include "SusyCommon/SodAccessor.h"

#include <string>
#include <vector>

using namespace std;
using namespace Susy;

int main(int argc, char **argv)
{
    string sodName = "";
    Susy::SodAccessor soda(sodName);
    cout<<"Accessing the private attributes through the accessor:"<<endl
        <<"m_dataSource : "<<soda.dataSource ()<<endl
        <<"m_doJetArea  : "<<soda.doJetArea  ()<<endl
        <<"m_doJetGSC   : "<<soda.doJetGSC   ()<<endl
        <<"m_is8TeV     : "<<soda.is8TeV     ()<<endl
        <<"m_isDerived  : "<<soda.isDerived  ()<<endl
        <<"m_jesNPset   : "<<soda.jesNPset   ()<<endl
        <<"m_eleTerm    : "<<soda.eleTerm    ()<<endl
        <<"m_gammaTerm  : "<<soda.gammaTerm  ()<<endl
        <<"m_tauTerm    : "<<soda.tauTerm    ()<<endl
        <<"m_jetTerm    : "<<soda.jetTerm    ()<<endl
        <<"m_muonTerm   : "<<soda.muonTerm   ()<<endl
        <<endl;
}
