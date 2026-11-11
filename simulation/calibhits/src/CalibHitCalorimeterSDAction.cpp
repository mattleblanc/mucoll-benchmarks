//==========================================================================
//  CalibHitCalorimeterSDAction
//--------------------------------------------------------------------------
//  DDG4 sensitive action producing ATLAS-style "calibration hits" for
//  sampling calorimeters, alongside the standard simulated hits.
//
//  For each calorimeter it writes up to seven collections:
//
//   <ReadoutName>                   Standard active-layer hits, identical to
//                                   Geant4ScintillatorCalorimeterAction with
//                                   an edep>0 filter (Birks' law optional).
//
//   <ReadoutName><CalibSuffix>      Calibration hits: the raw Geant4 energy
//                                   deposited anywhere in the cell volume,
//                                   i.e. active slice + absorber + passive
//                                   slices of the same layer, keyed by the
//                                   cellID of the active cell. No Birks.
//                                   Includes invisible energy (see below)
//                                   when IncludeInvisibleEnergy is set.
//
//   <ReadoutName><InvisibleSuffix>  The invisible component alone (energy
//                                   only; can be negative, e.g. n capture).
//
//   <ReadoutName><EMSuffix>         Visible energy deposited by e+-/gamma
//                                   tracks (energy only): the EM component.
//                                   Non-EM visible = Calib - Invisible - EM.
//
//   ...<PromptSuffix>               Prompt versions of the three calibration
//                                   collections above (energy only), e.g.
//                                   <ReadoutName>CalibPrompt: only steps with
//                                   PromptWindowMin < t - |cell centre|/c < PromptWindowMax,
//                                   i.e. the same timing cut as the k4RecCalorimeter
//                                   digitisers (BaseDigitiser::standardIntegration)
//                                   apply to the active hits.
//
//  Invisible energy is computed per step from the energy balance of the
//  stepping track and its secondaries (see measurableEnergy()). It covers
//  nuclear binding energy and neutrinos (counted in the cell where they are
//  produced), as in ATLAS calibration hits.
//
//  The geometry is NOT modified at the DD4hep level. At the start of the
//  first event the action attaches its Geant4 sensitive detector to every
//  sibling volume of the active slice inside each layer (Geant4-level only,
//  so the DD4hep volume manager never sees duplicate volume IDs). Steps in
//  those passive slices are mapped onto the active slice of the same layer
//  by swapping the leaf of the navigation history, and the standard
//  segmentation then gives the cellID.
//
//  Not included: energy deposited in envelopes outside the layer volumes.
//==========================================================================

#include <DD4hep/InstanceCount.h>
#include <DD4hep/Printout.h>
#include <DD4hep/VolumeManager.h>
#include <DDG4/Factories.h>
#include <DDG4/Geant4Data.h>
#include <DDG4/Geant4HitCollection.h>
#include <DDG4/Geant4SensDetAction.h>
#include <DDG4/Geant4StepHandler.h>

#include "EnergyBalance.h"

#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4Material.hh>
#include <G4NavigationHistory.hh>
#include <G4ParticleDefinition.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4TouchableHistory.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSensitiveDetector.hh>

#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace dd4hep {
  namespace sim {

    class CalibHitCalorimeterSDAction : public Geant4Sensitive {
      using Hit          = Geant4Calorimeter::Hit;
      using Contribution = Geant4HitData::Contribution;

      // ---- Properties
      /// Suffix appended to the readout name for the calibration-hit collection
      std::string m_calibSuffix       { "Calib" };
      /// Apply Birks' law to the *active* collection (as Geant4ScintillatorCalorimeterAction)
      bool        m_applyBirks        { true };
      /// Store MC contributions for calibration hits at all (false: energy only)
      bool        m_storeCalibContrib { true };
      /// Merge calibration-hit contributions per Geant4 track within a cell
      bool        m_mergeCalibContrib { true };
      /// Add invisible energy (per-step energy balance) to the calibration hits
      bool        m_includeInvisible  { true };
      /// Suffix for the invisible-energy-only collection (energies only, can be negative)
      std::string m_invisibleSuffix   { "CalibInvisible" };
      /// Suffix for the EM (e+-/gamma) visible component (energies only; "" disables it)
      std::string m_emSuffix          { "CalibEM" };
      /// Also write prompt versions of the calibration collections
      bool        m_includePrompt     { true };
      /// Suffix appended to each calibration collection name for its prompt version
      std::string m_promptSuffix      { "Prompt" };
      /// Prompt window [ns] on t - |cell centre| / c (match the digitiser's timingWindowMin/Max)
      double      m_promptMin         { -0.5 };
      double      m_promptMax         {  1.0 };
      /// Subtract the light-speed time of flight from the IP (digitiser: timingCorrectForPropagation)
      bool        m_promptCorrectTOF  { true };

      // ---- State
      std::size_t m_activeColl { 0 };
      std::size_t m_calibColl  { 0 };
      std::size_t m_invisColl  { 0 };
      std::size_t m_emColl     { 0 };
      std::size_t m_calibPromptColl { 0 };
      std::size_t m_invisPromptColl { 0 };
      std::size_t m_emPromptColl    { 0 };
      /// cellID -> time of flight [ns] from the IP to the cell centre (cleared every event)
      std::map<VolumeID, double> m_tofCache;
      bool        m_attached   { false };
      std::size_t m_nMapErrors { 0 };
      /// Passive logical volumes we attached ourselves to
      std::set<const G4LogicalVolume*> m_passiveLVs;
      /// Layer (mother) logical volume -> active slice placement inside it
      std::map<const G4LogicalVolume*, G4VPhysicalVolume*> m_activeSibling;
      /// Subdetector volume manager (DDCore), for cell-centre positions
      VolumeManager m_volMgr;

    public:
      CalibHitCalorimeterSDAction(Geant4Context* ctxt, const std::string& nam,
                                  DetElement det, Detector& description)
        : Geant4Sensitive(ctxt, nam, det, description)
      {
        declareProperty("CalibCollectionSuffix",   m_calibSuffix);
        declareProperty("ApplyBirksLaw",           m_applyBirks);
        declareProperty("StoreCalibContributions", m_storeCalibContrib);
        declareProperty("MergeCalibContributions", m_mergeCalibContrib);
        declareProperty("IncludeInvisibleEnergy",  m_includeInvisible);
        declareProperty("InvisibleCollectionSuffix", m_invisibleSuffix);
        declareProperty("EMCollectionSuffix",      m_emSuffix);
        declareProperty("IncludePrompt",           m_includePrompt);
        declareProperty("PromptSuffix",            m_promptSuffix);
        declareProperty("PromptWindowMin",         m_promptMin);
        declareProperty("PromptWindowMax",         m_promptMax);
        declareProperty("PromptCorrectForPropagation", m_promptCorrectTOF);
        InstanceCount::increment(this);
      }

      virtual ~CalibHitCalorimeterSDAction() {
        InstanceCount::decrement(this);
      }

      /// Standard active hits, calibration hits, their invisible and EM components,
      /// and (optionally) prompt versions of the three calibration collections
      virtual void defineCollections() override {
        const std::string ro = m_readout.name();
        m_activeColl = defineCollection<Hit>(ro);
        m_calibColl  = defineCollection<Hit>(ro + m_calibSuffix);
        if ( m_includeInvisible )
          m_invisColl = defineCollection<Hit>(ro + m_invisibleSuffix);
        if ( !m_emSuffix.empty() )
          m_emColl = defineCollection<Hit>(ro + m_emSuffix);
        if ( m_includePrompt )  {
          m_calibPromptColl = defineCollection<Hit>(ro + m_calibSuffix + m_promptSuffix);
          if ( m_includeInvisible )
            m_invisPromptColl = defineCollection<Hit>(ro + m_invisibleSuffix + m_promptSuffix);
          if ( !m_emSuffix.empty() )
            m_emPromptColl = defineCollection<Hit>(ro + m_emSuffix + m_promptSuffix);
        }
      }

      virtual void begin(G4HCofThisEvent* hce) override {
        Geant4Sensitive::begin(hce);
        m_tofCache.clear();
        if ( !m_attached )  {
          attachPassiveVolumes();
          if ( m_includePrompt )
            always("+++ %s: prompt window %.3g < t - %s < %.3g ns", m_sensitive.name(), m_promptMin,
                   m_promptCorrectTOF ? "|r_cell|/c" : "0", m_promptMax);
        }
      }

      virtual bool process(const G4Step* step, G4TouchableHistory*) override {
        const G4VTouchable*     touch   = step->GetPreStepPoint()->GetTouchable();
        const G4LogicalVolume*  lv      = touch->GetVolume(0)->GetLogicalVolume();
        const double            edep    = step->GetTotalEnergyDeposit();
        const double            einv    = m_includeInvisible ? calibhits::invisibleEnergy(step) : 0.0;
        const bool              passive = m_passiveLVs.find(lv) != m_passiveLVs.end();
        if ( edep <= 0.0 && einv == 0.0 ) return true;
        Geant4StepHandler       h(step);

        VolumeID cell = 0;
        if ( !passive )  {
          // ---------------- Active slice: standard cellID
          try {
            cell = cellID(step);
          } catch ( const std::exception& e ) {
            error("+++ cellID failed for active step: %s", e.what());
            return true;
          }
          // Standard active hit: identical to Geant4ScintillatorCalorimeterAction + "edep0" filter
          if ( edep > 0.0 )  {
            Contribution c = Hit::extractContribution(step, m_applyBirks);
            addActiveHit(cell, c, h);
            mark(h.track);
          }
        }
        else  {
          // ---------------- Passive slice: map onto active slice of same layer
          auto it = m_activeSibling.find(touch->GetVolume(1)->GetLogicalVolume());
          if ( it == m_activeSibling.end() )  {
            if ( m_nMapErrors++ < 10 )
              error("+++ No active sibling for passive volume %s", lv->GetName().c_str());
            return true;
          }
          G4NavigationHistory nav(*touch->GetHistory());
          nav.BackLevel();
          nav.NewLevel(it->second, kNormal, it->second->GetCopyNo());
          G4TouchableHistory activeTouch(nav);
          const G4ThreeVector global = 0.5 * (h.prePosG4() + h.postPosG4());
          try {
            cell = cellID(&activeTouch, global);
          } catch ( const std::exception& e ) {
            if ( m_nMapErrors++ < 10 )
              error("+++ cellID failed for passive step in %s: %s", lv->GetName().c_str(), e.what());
            return true;
          }
        }

        // ---------------- Calibration hit: raw visible + invisible energy
        Contribution c = Hit::extractContribution(step, false);
        c.deposit += einv;
        addCalibHit(cell, c);
        bool isEM = false;
        if ( edep > 0.0 && !m_emSuffix.empty() )  {
          const int pdg = step->GetTrack()->GetDefinition()->GetPDGEncoding();
          isEM = (pdg == 22 || pdg == 11 || pdg == -11);
        }
        if ( einv != 0.0 ) addEnergyHit(m_invisColl, cell, einv);
        if ( isEM ) addEnergyHit(m_emColl, cell, edep);

        // ---------------- Prompt versions: the digitiser's timing cut, per step
        // (c.time is the same track time the active-hit contributions carry)
        if ( m_includePrompt && inPromptWindow(cell, c.time) )  {
          addEnergyHit(m_calibPromptColl, cell, c.deposit);
          if ( einv != 0.0 && m_includeInvisible ) addEnergyHit(m_invisPromptColl, cell, einv);
          if ( isEM ) addEnergyHit(m_emPromptColl, cell, edep);
        }
        return true;
      }

    private:
      /// Attach our G4 sensitive detector to all siblings of the active slices
      void attachPassiveVolumes()  {
        m_attached = true;
        G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
        const std::string sdName = m_sensitive.name();

        G4VSensitiveDetector*       mySD = nullptr;
        std::set<G4LogicalVolume*>  activeLVs;
        for ( G4LogicalVolume* lv : *store )  {
          G4VSensitiveDetector* sd = lv->GetSensitiveDetector();
          if ( sd && sd->GetName() == sdName )  {
            activeLVs.insert(lv);
            mySD = sd;
          }
        }
        if ( activeLVs.empty() || !mySD )  {
          except("+++ No active logical volumes found for sensitive detector %s", sdName.c_str());
        }

        std::map<std::string, std::size_t> passiveMaterials;
        std::size_t nMothers = 0, nMultiActive = 0;
        for ( G4LogicalVolume* mother : *store )  {
          G4VPhysicalVolume* activePV = nullptr;
          const std::size_t nd = mother->GetNoDaughters();
          for ( std::size_t i = 0; i < nd; ++i )  {
            G4VPhysicalVolume* d = mother->GetDaughter(i);
            if ( activeLVs.count(d->GetLogicalVolume()) )  {
              if ( activePV ) ++nMultiActive;
              else activePV = d;
            }
          }
          if ( !activePV ) continue;
          ++nMothers;
          m_activeSibling[mother] = activePV;
          for ( std::size_t i = 0; i < nd; ++i )  {
            G4LogicalVolume* dlv = mother->GetDaughter(i)->GetLogicalVolume();
            if ( activeLVs.count(dlv) ) continue;
            G4VSensitiveDetector* other = dlv->GetSensitiveDetector();
            if ( other && other != mySD )  {
              warning("+++ %s: sibling %s already has SD %s -- not attached",
                      sdName.c_str(), dlv->GetName().c_str(), other->GetName().c_str());
              continue;
            }
            if ( dlv->GetNoDaughters() > 0 )  {
              warning("+++ %s: sibling %s has %zu daughters; only its own material is counted",
                      sdName.c_str(), dlv->GetName().c_str(), dlv->GetNoDaughters());
            }
            dlv->SetSensitiveDetector(mySD);
            if ( m_passiveLVs.insert(dlv).second )
              ++passiveMaterials[dlv->GetMaterial()->GetName()];
          }
        }
        if ( nMultiActive )  {
          warning("+++ %s: %zu layer volumes contain more than one active slice; "
                  "passive energy is assigned to the first one", sdName.c_str(), nMultiActive);
        }
        std::stringstream mats;
        for ( const auto& [mat, n] : passiveMaterials ) mats << " " << mat << "(" << n << ")";
        always("+++ %s: calibration hits enabled. %zu active LVs, %zu layer LVs, "
               "%zu passive LVs attached. Passive materials:%s",
               sdName.c_str(), activeLVs.size(), nMothers, m_passiveLVs.size(), mats.str().c_str());
      }

      /// Identical to the standard DDG4 calorimeter hit handling
      void addActiveHit(VolumeID cell, const Contribution& contrib, const Geant4StepHandler& h)  {
        Geant4HitCollection* coll = collection(m_activeColl);
        Hit* hit = coll->findByKey<Hit>(cell);
        if ( !hit )  {
          DDSegmentation::Vector3D pos;
          Position global;
          if ( !useVolumeManager() || !m_segmentation.isValid() )  {
            pos    = h.avgPosition();
            global = h.localToGlobal(pos);
          }
          else if ( !m_segmentation.cellsSpanVolumes() )  {
            pos    = m_segmentation.position(cell);
            global = h.localToGlobal(pos);
          }
          else  {
            global = cellCentre(cell);
          }
          hit = new Hit(global);
          hit->cellID = cell;
          coll->add(cell, hit);
        }
        hit->truth.emplace_back(contrib);
        hit->energyDeposit += contrib.deposit;
      }

      /// Calibration hit: position is always the centre of the active cell
      void addCalibHit(VolumeID cell, const Contribution& c)  {
        Geant4HitCollection* coll = collection(m_calibColl);
        Hit* hit = coll->findByKey<Hit>(cell);
        if ( !hit )  {
          hit = new Hit(cellCentre(cell));
          hit->cellID = cell;
          coll->add(cell, hit);
        }
        hit->energyDeposit += c.deposit;
        if ( !m_storeCalibContrib ) return;
        if ( m_mergeCalibContrib )  {
          for ( auto& t : hit->truth )  {
            if ( t.trackID == c.trackID )  {
              t.deposit += c.deposit;
              t.length  += c.length;
              t.time     = std::min(t.time, c.time);
              return;
            }
          }
        }
        hit->truth.emplace_back(c);
      }

      /// Energy-only hit (invisible or EM component; no contributions)
      void addEnergyHit(std::size_t which, VolumeID cell, double e)  {
        Geant4HitCollection* coll = collection(which);
        Hit* hit = coll->findByKey<Hit>(cell);
        if ( !hit )  {
          hit = new Hit(cellCentre(cell));
          hit->cellID = cell;
          coll->add(cell, hit);
        }
        hit->energyDeposit += e;
      }

      /// Digitiser timing cut: min < t - |cell centre|/c < max (strict, as in BaseDigitiser)
      bool inPromptWindow(VolumeID cell, double t)  {
        double tof = 0.0;
        if ( m_promptCorrectTOF )  {
          auto it = m_tofCache.find(cell);
          if ( it == m_tofCache.end() )  {
            const Position p = cellCentre(cell);          // [mm]
            it = m_tofCache.emplace(cell, std::sqrt(p.Mag2()) / CLHEP::c_light).first;
          }
          tof = it->second;
        }
        const double rel = t / CLHEP::ns - tof;
        return rel > m_promptMin && rel < m_promptMax;
      }

      /// Global position [mm] of the centre of a cell, via the DDCore volume manager
      Position cellCentre(VolumeID cell)  {
        if ( !m_volMgr.isValid() )  {
          VolumeManager glob = VolumeManager::getVolumeManager(detectorDescription());
          m_volMgr = glob.subdetector(id());
        }
        VolumeManagerContext* vc = m_volMgr.lookupContext(m_segmentation.volumeID(cell));
        DDSegmentation::Vector3D pos = m_segmentation.position(cell);
        return vc->localToWorld(Position(pos)) / dd4hep::mm;
      }
    };
  }  // namespace sim
}  // namespace dd4hep

DECLARE_GEANT4SENSITIVE(CalibHitCalorimeterSDAction)
