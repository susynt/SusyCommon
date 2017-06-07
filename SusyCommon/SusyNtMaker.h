#ifndef SUSYCOMMON_SUSYNTMAKER_H
#define SUSYCOMMON_SUSYNTMAKER_H

//SusyNtuple
#include "SusyNtuple/SusyNtObject.h"
#include "SusyNtuple/TriggerTools.h"

//SusyCommon
#include "SusyCommon/XaodAnalysis.h"
#include "SusyCommon/SystematicMapping.h"

//ROOT
#include "TStopwatch.h"

//std/stl
#include <iostream>
#include <string>
#include <vector>

//Tools
namespace Root { class TElectronEfficiencyCorrectionTool; }


namespace Susy {
class SusyNtMaker : public XaodAnalysis
{
    public :

        //////////////////////////////////////////////////////////////////////
        // TSelector
        //////////////////////////////////////////////////////////////////////
        virtual void SlaveBegin(TTree* tree);
        virtual Bool_t Process(Long64_t entry);
        virtual void Terminate();

        //////////////////////////////////////////////////////////////////////
        // SusyNtMaker
        //////////////////////////////////////////////////////////////////////
        SusyNtMaker();
        virtual ~SusyNtMaker(){};

    private :

        TFile* m_file_outtree;
        TTree* m_outtree; 

        Susy::SusyNtObject m_susyNt; // SusyNt interface

        // timer
        TStopwatch m_timer;

        //bool initialize_output_tree();
        //bool save_output_tree();
        ////bool initialize_cutflow_histograms();

        //void write_metadata();


}; // clas SusyNtMaker
} // namespace Susy


#endif
