//==========================================================================
//  EnergyBalanceCheckAction
//--------------------------------------------------------------------------
//  Debug stepping action for validating the invisible-energy bookkeeping of
//  CalibHitCalorimeterSDAction. It applies the same per-step energy balance
//  (EnergyBalance.h) to EVERY step in the world and counts the measurable
//  energy of tracks leaving the world. By energy conservation, per event
//
//      sum(visible) + sum(invisible) + escaped  ==  measurable E of primaries
//
//  Any deficit means steps or tracks that the bookkeeping misses (e.g. tracks
//  killed without being stepped). Prints one line per event and a run summary.
//
//  With OutputFile set, it also writes per event the visible and invisible
//  energy in each top-level volume (the daughters of the world: one per
//  subdetector, e.g. ECalBarrel, Solenoid, OuterTracker; "world" for the world
//  volume itself), plus the primary energy and the escaped energy, as CSV:
//      event,region,visible_GeV,invisible_GeV
//  with special regions "E0" (primary energy) and "escaped".
//
//  Enable in a steering file with:
//      SIM.action.step = ["EnergyBalanceCheckAction/EBCheck"]
//  or, with the per-region CSV:
//      SIM.action.step = {"name": "EnergyBalanceCheckAction/EBCheck",
//                         "parameter": {"OutputFile": "energy_balance.csv"}}
//==========================================================================
#include <DD4hep/InstanceCount.h>
#include <DDG4/Factories.h>
#include <DDG4/Geant4Context.h>
#include <DDG4/Geant4EventAction.h>
#include <DDG4/Geant4Kernel.h>
#include <DDG4/Geant4SteppingAction.h>

#include <G4Event.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <G4Step.hh>

#include <G4VTouchable.hh>
#include <G4VPhysicalVolume.hh>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <utility>

#include "EnergyBalance.h"

namespace dd4hep {
  namespace sim {

    class EnergyBalanceCheckAction : public Geant4SteppingAction {
      double m_vis { 0 }, m_inv { 0 }, m_esc { 0 };
      double m_maxRelDev { 0 };
      long   m_nEvents { 0 };
      std::string   m_outputFile;
      std::ofstream m_out;
      /// region -> (visible, invisible)
      std::map<std::string, std::pair<double, double>> m_regions;
      std::map<const G4VPhysicalVolume*, std::string>  m_regionName;

      const std::string& regionOf(const G4Step* step)  {
        static const std::string world = "world";
        const G4VTouchable* t = step->GetPreStepPoint()->GetTouchable();
        const int depth = t->GetHistoryDepth();
        if ( depth < 1 ) return world;
        const G4VPhysicalVolume* top = t->GetVolume(depth - 1);   // daughter of the world
        auto it = m_regionName.find(top);
        if ( it == m_regionName.end() )
          it = m_regionName.emplace(top, subdetectorName(top->GetName())).first;
        return it->second;
      }

      /// "Solenoid_10" -> "Solenoid"; assembly imprints such as
      /// "AV_68!InnerTrackers_8#8!InnerTrackerBarrelModule_01_23#23" -> "InnerTrackers"
      static std::string subdetectorName(std::string n)  {
        if ( n.rfind("AV_", 0) == 0 )  {
          const auto a = n.find('!');
          const auto b = n.find('#', a + 1);
          if ( a != std::string::npos && b != std::string::npos ) n = n.substr(a + 1, b - a - 1);
        }
        const auto u = n.find_last_of('_');
        if ( u != std::string::npos && u + 1 < n.size() &&
             n.find_first_not_of("0123456789", u + 1) == std::string::npos )
          n = n.substr(0, u);
        return n;
      }

    public:
      EnergyBalanceCheckAction(Geant4Context* ctxt, const std::string& nam)
        : Geant4SteppingAction(ctxt, nam)
      {
        declareProperty("OutputFile", m_outputFile);
        context()->kernel().eventAction().callAtBegin(this, &EnergyBalanceCheckAction::beginEvent);
        context()->kernel().eventAction().callAtEnd(this, &EnergyBalanceCheckAction::endEvent);
        InstanceCount::increment(this);
      }

      virtual ~EnergyBalanceCheckAction()  {
        always("+++ EnergyBalance run summary: %ld events, max |sum/E0 - 1| = %.3e",
               m_nEvents, m_maxRelDev);
        InstanceCount::decrement(this);
      }

      void beginEvent(const G4Event*)  {
        m_vis = m_inv = m_esc = 0.0;
        m_regions.clear();
      }

      virtual void operator()(const G4Step* step, G4SteppingManager*) override  {
        const double vis = step->GetTotalEnergyDeposit();
        const double inv = calibhits::invisibleEnergy(step);
        m_vis += vis;
        m_inv += inv;
        if ( !m_outputFile.empty() && (vis != 0.0 || inv != 0.0) )  {
          auto& r = m_regions[regionOf(step)];
          r.first  += vis;
          r.second += inv;
        }
        if ( step->GetPostStepPoint()->GetStepStatus() == fWorldBoundary )
          m_esc += calibhits::measurableEnergy(step->GetTrack()->GetDefinition(),
                                               step->GetPostStepPoint()->GetKineticEnergy());
      }

      void endEvent(const G4Event* evt)  {
        double e0 = 0.0;
        for ( int iv = 0; iv < evt->GetNumberOfPrimaryVertex(); ++iv )  {
          for ( const G4PrimaryParticle* p = evt->GetPrimaryVertex(iv)->GetPrimary(); p; p = p->GetNext() )  {
            if ( p->GetG4code() ) e0 += calibhits::measurableEnergy(p->GetG4code(), p->GetKineticEnergy());
          }
        }
        const double sum = m_vis + m_inv + m_esc;
        const double rel = e0 > 0 ? sum / e0 - 1.0 : 0.0;
        m_maxRelDev = std::max(m_maxRelDev, std::abs(rel));
        ++m_nEvents;
        always("+++ EnergyBalance event %d: E0 %.4f GeV | visible %.4f  invisible %.4f  escaped %.4f GeV"
               " | sum/E0 - 1 = %+.3e", evt->GetEventID(), e0 / CLHEP::GeV, m_vis / CLHEP::GeV,
               m_inv / CLHEP::GeV, m_esc / CLHEP::GeV, rel);
        if ( !m_outputFile.empty() )  {
          if ( !m_out.is_open() )  {
            m_out.open(m_outputFile);
            m_out << "event,region,visible_GeV,invisible_GeV\n";
          }
          const int id = evt->GetEventID();
          m_out.precision(9);
          m_out << id << ",E0," << e0 / CLHEP::GeV << ",0\n";
          m_out << id << ",escaped," << m_esc / CLHEP::GeV << ",0\n";
          for ( const auto& [name, e] : m_regions )
            m_out << id << "," << name << "," << e.first / CLHEP::GeV << "," << e.second / CLHEP::GeV << "\n";
          m_out.flush();
        }
      }
    };
  }  // namespace sim
}  // namespace dd4hep

DECLARE_GEANT4ACTION(EnergyBalanceCheckAction)
