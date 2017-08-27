###################################################################
# Package: SusyCommon
###################################################################

# declare the package name
atlas_subdir( SusyCommon )

# extra deps based on the build environment:
set( extra_deps )
set( extra_libs )

# declare the SusyCommon package dependencies
atlas_depends_on_subdirs(
    PUBLIC
    # collect local SUSYTools
    SUSYTools
    #SUSYPhys/SUSYTools
    # egamma
    PhysicsAnalysis/ElectronPhotonID/ElectronEfficiencyCorrection
    PhysicsAnalysis/ElectronPhotonID/ElectronPhotonSelectorTools
    # billy butt
    PhysicsAnalysis/AnalysisCommon/PileupReweighting 
    # muon
    PhysicsAnalysis/MuonID/MuonIDAnalysis/MuonEfficiencyCorrections
    # some people worry about these
    PhysicsAnalysis/TauID/TauAnalysisTools
    # DQ
    DataQuality/GoodRunsLists
    # trigger
    Trigger/TrigAnalysis/TrigDecisionTool
    Trigger/TrigConfiguration/TrigConfxAOD
    # isolation
    PhysicsAnalysis/AnalysisCommon/IsolationSelection
    # xAOD/Core
    Event/xAOD/xAODCore
    Event/xAOD/xAODCutFlow
    Control/AthToolSupport/AsgTools
    Tools/PathResolver
    PRIVATE
    # this is us
    SusyNtuple
    ${extra_deps}
)

# external dependencies
find_package( ROOT COMPONENTS Gpad Graf Graf3d Core Tree MathCore Hist Tree RIO )

atlas_add_root_dictionary( SusyCommonLib SusyCommonCintDict 
            ROOT_HEADERS SusyCommon/*.h Root/LinkDef.h
            EXTERNAL_PACKAGES ROOT
)

# libraries
atlas_add_library( SusyCommonLib 
    SusyCommon/*.h Root/*.cxx ${SusyCommonCintDict}
    PUBLIC_HEADERS SusyCommon 
    PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
    LINK_LIBRARIES SUSYToolsLib SusyNtupleLib AsgTools xAODCore xAODCutFlow
    ElectronEfficiencyCorrectionLib ElectronPhotonSelectorToolsLib
    PileupReweightingLib
    MuonEfficiencyCorrectionsLib
    TauAnalysisToolsLib
    GoodRunsListsLib
    TrigDecisionToolLib
    TrigConfxAODLib
    IsolationSelectionLib
    PathResolver
    PRIVATE_LINK_LIBRARIES ${ROOT_LIBRARIES}
)

# executable(s) in the package
set( extra_libs )

function(SusyNtExec filename)
    set(execname)
    get_filename_component(execname ${filename} NAME_WE)
    atlas_add_executable( ${execname} "util/${execname}.cxx"
        INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
        LINK_LIBRARIES ${ROOT_LIBRARIES} SusyCommonLib ${extra_libs}
    )
endfunction(SusyNtExec)

# build everything in util/ as an executable for SusyCommon 
file(GLOB files "util/*.cxx")
foreach(file ${files})
    SusyNtExec(${file})
endforeach()

# add our data files so that they can be picked up by PathResolver's awesome abilities
atlas_install_data( data/* )

