## ddsim steering for calorimeter "calibration hits" (truth deposited energy per cell).
##
## Identical to steer_baseline.py, except that the four MAIA calorimeters use the
## CalibHitCalorimeterSDAction plugin (simulation/calibhits/). Each one then writes
##   <Readout>       e.g. ECalBarrelCollection       -- unchanged active-layer hits
##   <Readout>Calib           -- total deposited energy per cell (active + absorber +
##                               passive slices; visible + invisible)
##   <Readout>CalibInvisible  -- invisible component (nuclear binding, neutrinos)
##   <Readout>CalibEM         -- visible energy from e+-/gamma (EM component)
##   <Readout>Calib*Prompt    -- the same three, keeping only energy inside the
##                               digitiser's timing window (see PROMPT_WINDOWS)
## The geometry is the production one ($MUCOLL_GEO); digitisation and reconstruction
## are unchanged and simply carry the *Calib collections along.
##
## Requires the plugin to be built and on LD_LIBRARY_PATH:
##   simulation/calibhits/build.sh && source simulation/calibhits/setup.sh
##
## ddsim exec()s steering files without __file__, so steer_baseline.py is located
## relative to the steering-file path that ddsim was given, falling back to
## $MUCOLL_BENCHMARKS/simulation (set by setup_config.sh).
import os

_baseline = os.path.join(os.path.dirname(os.path.abspath(SIM.steeringFile or ".")), "steer_baseline.py")
if not os.path.isfile(_baseline):
    _baseline = os.path.join(os.environ.get("MUCOLL_BENCHMARKS", ""), "simulation", "steer_baseline.py")
_ns = {}
exec(compile(open(_baseline).read(), _baseline, "exec"), _ns)
SIM = _ns["SIM"]

CALIB_DETECTORS = ["ECalBarrel", "ECalEndcap", "HCalBarrel", "HCalEndcap"]
CALIB_ACTION = ("CalibHitCalorimeterSDAction", {
    "CalibCollectionSuffix": "Calib",
    ## Same as the baseline Geant4ScintillatorCalorimeterAction for the active hits
    "ApplyBirksLaw": True,
    ## MC-particle contributions for the *Calib hits. Off by default: the EM
    ## fraction comes from the *CalibEM collection instead. Test with 2 x 20 GeV
    ## pi-: baseline 0.45 MB, calib with contributions 5.4 MB, without 1.0 MB.
    "StoreCalibContributions": False,
    ## One contribution per Geant4 track and cell in the *Calib collections
    ## (otherwise one per step, which is very large in the absorber)
    "MergeCalibContributions": True,
    ## Add invisible energy (nuclear binding, neutrinos) to the *Calib hits, and
    ## write the invisible component alone to <Readout>CalibInvisible
    "IncludeInvisibleEnergy": True,
    ## Visible energy deposited by e+-/gamma per cell -> <Readout>CalibEM ("" disables)
    "EMCollectionSuffix": "CalibEM",
    ## Prompt versions of the three calibration collections (<...>Prompt), with the
    ## digitiser's timing cut applied per step: window on t - |cell centre|/c
    "IncludePrompt": True,
    "PromptSuffix": "Prompt",
    "PromptCorrectForPropagation": True,   # digitiser: timingCorrectForPropagation = 1
})

## Prompt windows [ns] per calorimeter: keep equal to timingWindowMin/Max of the
## calorimeter digitisers in configs/MAIAConfig/MAIAConfig/CaloDigi/calorimetry_{EM,HAD}.py
PROMPT_WINDOWS = {
    "ECalBarrel": (-0.5, 1.0), "ECalEndcap": (-0.5, 1.0),
    "HCalBarrel": (-0.5, 10.0), "HCalEndcap": (-0.5, 10.0),
}

## Closure of the energy bookkeeping over the whole world, per event
## (visible + invisible + escaped == primary energy), with the energy per
## top-level detector volume written to a CSV. Enabled by setting
## CALIBHITS_ENERGY_BALANCE_CSV=<file> (run_chain_calib.sh --energy-balance).
if os.environ.get("CALIBHITS_ENERGY_BALANCE_CSV"):
    SIM.action.step = {"name": "EnergyBalanceCheckAction/EBCheck",
                       "parameter": {"OutputFile": os.environ["CALIBHITS_ENERGY_BALANCE_CSV"]}}

for _det in CALIB_DETECTORS:
    _params = dict(CALIB_ACTION[1])
    _params["PromptWindowMin"], _params["PromptWindowMax"] = PROMPT_WINDOWS[_det]
    SIM.action.mapActions[_det] = (CALIB_ACTION[0], _params)
    ## The action applies edep>0 to the active hits itself. A sequence-level
    ## filter must not be used, so that later (invisible-energy) steps with
    ## zero visible deposit can still reach it.
    SIM.filter.mapDetFilter[_det] = None
