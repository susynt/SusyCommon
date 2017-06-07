//  -*- c++ -*-
#ifndef SUSYNTUPLE_SS3L_CHARGEFLIP_H
#define SUSYNTUPLE_SS3L_CHARGEFLIP_H

// version 1.1, 21/09/2015

#include "xAODEgamma/ElectronContainer.h"
#include "xAODTruth/TruthParticleContainer.h"

void fillElectronChargeFlip(xAOD::ElectronContainer* reco_electrons,
                            const xAOD::TruthParticleContainer* container_particles,
                            int mcChannelNumber);
#endif
