//==========================================================================
//  EnergyBalance.h -- per-step energy bookkeeping shared by the calibration-hit
//  sensitive action and the EnergyBalanceCheckAction debug stepping action.
//==========================================================================
#ifndef CALIBHITS_ENERGYBALANCE_H
#define CALIBHITS_ENERGYBALANCE_H

#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <G4ParticleDefinition.hh>
#include <G4Step.hh>
#include <G4Track.hh>

#include <cmath>

namespace calibhits {

  /// "Measurable" energy of a particle, following the ATLAS calibration-hit
  /// convention (CaloG4Sim): rest masses that cannot be released in the
  /// calorimeter are not counted, so that nuclear binding energy and energy
  /// carried by neutrinos show up as invisible energy.
  ///   e-, p, n, nuclei      : kinetic energy (they exist in ordinary matter)
  ///   e+, anti-p, anti-n    : kinetic + 2 m (annihilation releases both masses)
  ///   other baryons         : kinetic + m - m_p   (decay to a nucleon)
  ///   other antibaryons     : kinetic + m + m_p
  ///   neutrinos             : 0 (never subtracted: their energy is invisible
  ///                           in the cell where they are produced)
  ///   everything else       : total energy (gamma, mesons, muons, ...)
  inline double measurableEnergy(const G4ParticleDefinition* def, double ekin)  {
    const int    pdg  = def->GetPDGEncoding();
    const int    apdg = std::abs(pdg);
    const double m    = def->GetPDGMass();
    if ( apdg == 12 || apdg == 14 || apdg == 16 )   return 0.0;
    if ( pdg == 11 || pdg == 2212 || pdg == 2112 )   return ekin;
    if ( pdg == -11 || pdg == -2212 || pdg == -2112 ) return ekin + 2.0 * m;
    if ( def->GetParticleType() == "nucleus" )        return ekin;
    const double B = def->GetBaryonNumber();
    if ( B > 0 ) return ekin + m - CLHEP::proton_mass_c2;
    if ( B < 0 ) return ekin + m + CLHEP::proton_mass_c2;
    return ekin + m;
  }

  /// Per-step energy balance: measurable energy lost by the track, minus
  /// the visible deposit, minus the measurable energy given to secondaries.
  inline double invisibleEnergy(const G4Step* step)  {
    const G4Track*               trk = step->GetTrack();
    const G4ParticleDefinition*  def = trk->GetDefinition();
    const G4TrackStatus          st  = trk->GetTrackStatus();
    // Leaving the world is an escape (counted by EnergyBalanceCheckAction), not invisible energy
    const bool exits  = step->GetPostStepPoint()->GetStepStatus() == fWorldBoundary;
    const bool killed = !exits && (st == fStopAndKill || st == fKillTrackAndSecondaries);
    const double ekinPre  = step->GetPreStepPoint()->GetKineticEnergy();
    const double ekinPost = step->GetPostStepPoint()->GetKineticEnergy();
    // Alive: the mass terms cancel, so use the kinetic-energy difference (no rounding noise)
    double balance = killed ? measurableEnergy(def, ekinPre) : (ekinPre - ekinPost);
    if ( const auto* secs = step->GetSecondaryInCurrentStep() )  {
      for ( const G4Track* s : *secs )
        balance -= measurableEnergy(s->GetDefinition(), s->GetKineticEnergy());
    }
    balance -= step->GetTotalEnergyDeposit();
    // Remove floating-point noise from transport steps
    return std::abs(balance) < 1e-6 * CLHEP::keV ? 0.0 : balance;
  }

}  // namespace calibhits

#endif
