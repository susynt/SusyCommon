// Dear emacs, this is -*- c++ -*-
#ifndef SUSY_XAODANALYSIS_TYPES_H

/**
  Typedefs used to access the xaod collections

  Note to self: because of the 'typedef' used to version each object
  (e.g. Electron_v1 etc.) I cannot get to forward declare the
  Container collections. So this header should just go after the xAOD*
  includes.

  davide.gerbaudo@gmail.com
  Sept 2014
*/

#include <utility>

namespace susy{
  typedef std::pair<xAOD::MuonContainer*,     xAOD::ShallowAuxContainer*> MuonsWithAux_t;
  typedef std::pair<xAOD::ElectronContainer*, xAOD::ShallowAuxContainer*> ElectronsWithAux_t;
  typedef std::pair<xAOD::TauJetContainer*,   xAOD::ShallowAuxContainer*> TausWithAux_t;
  typedef std::pair<xAOD::JetContainer*,      xAOD::ShallowAuxContainer*> JetsWithAux_t;
  typedef std::pair<xAOD::PhotonContainer*,   xAOD::ShallowAuxContainer*> PhotonsWithAux_t;
} // susy
#endif

