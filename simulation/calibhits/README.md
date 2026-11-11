# Calorimeter calibration hits

ATLAS-style "calibration hits" for MAIA: the true energy deposited in each
calorimeter **cell volume** (active layer + absorber + passive layers), including
invisible energy (nuclear binding, neutrinos). These
are meant as the truth reference for cluster calibration, in the spirit of
E_dep^clus in [arXiv:2412.04370](https://arxiv.org/abs/2412.04370).

**To run it on Perlmutter, start with [QUICKSTART.md](QUICKSTART.md).**

With `steer_calib.py`, each of ECalBarrel, ECalEndcap, HCalBarrel and HCalEndcap
writes seven collections:

| Collection | Content |
|---|---|
| `ECalBarrelCollection` etc. | Standard active-layer hits. **Hit-for-hit identical** to `steer_baseline.py` with the same seed. |
| `ECalBarrelCollectionCalib` etc. | Total energy in the cell volume: raw Geant4 deposits in every slice of the layer (no Birks' law) plus invisible energy. Keyed by the cellID of the active cell. |
| `ECalBarrelCollectionCalibInvisible` etc. | The invisible component alone. It can be negative, e.g. when a neutron capture releases binding energy. |
| `ECalBarrelCollectionCalibEM` etc. | Visible energy deposited by e±/γ: the EM component. Non-EM visible = Calib − Invisible − EM. |
| `…CalibPrompt`, `…CalibInvisiblePrompt`, `…CalibEMPrompt` | The same three, keeping only the steps inside the digitizer's timing window (see below). |

**Prompt collections.** The calorimeter digitizers (k4RecCalorimeter
`BaseDigitiser::standardIntegration`) only integrate active-hit contributions
with `windowMin < t − |r_cell|/c < windowMax`. In MAIA that is −0.5…1 ns for the
ECal and −0.5…10 ns for the HCal. The `…Prompt` collections apply the same cut
to every step, using the same cell-centre time of flight.
- **Windows:** set per calorimeter in `PROMPT_WINDOWS` in `steer_calib.py`.
  Keep them equal to MAIAConfig's `timingWindowMin/Max`.
- **Which to use:** the prompt calibration energy is the deposited energy over
  the same time range that the reconstructed cells integrate. The all-time
  collections are kept for the energy budget and for studying the cut.
- **Validation:** recomputing the prompt energy from per-step times agrees with
  `CalibPrompt` to 10⁻⁷ GeV per cell. For 3 π± events (10–40 GeV), 4.3% of the
  calibration energy is outside the windows.

The geometry stays the production `$MUCOLL_GEO`. Digitisation and
reconstruction are unchanged and just carry the `*Calib` collections along.
BIB is overlaid from separate files, so the `*Calib` collections contain
hard-scatter energy only.

## How it works

`src/CalibHitCalorimeterSDAction.cpp` is a DDG4 sensitive action, used in place
of `Geant4ScintillatorCalorimeterAction` for the four calorimeters.

At the first event it attaches its Geant4 sensitive detector to every sibling
slice of the active slice in each layer. For example, the ECal has 7 passive
slices per layer: W, G10, HV, PCB and air. The attachment is done at the Geant4
level only. If the slices were marked sensitive in the compact XML instead, the
DD4hep volume manager would abort on duplicate volume IDs.

A step in a passive slice is assigned to the cell of the active slice in the
same layer: the action swaps the leaf of the navigation history and applies the
standard segmentation.

**Invisible energy** is computed per step from the energy balance of the
stepping track and its secondaries (`src/EnergyBalance.h`). It follows the ATLAS
`CaloG4Sim` convention for rest masses:

- e⁻, p, n and nuclei count their kinetic energy only, so nuclear binding
  energy shows up as invisible;
- e⁺ and antinucleons count kinetic energy + 2m, since annihilation releases both
  masses;
- neutrinos are never subtracted, so their energy is invisible in the cell
  where they are produced;
- everything else counts its total energy.

`EnergyBalanceCheckAction` (debug, enabled with `SIM.action.step`) applies the
same balance to every step in the world. With it, visible + invisible + escaped
energy equals the primary energy to 10⁻¹¹ per event (5 × 20 GeV e⁻ and π⁻).

## Build (once per container image)

The plugin must be built inside the **same image** that runs ddsim. Rebuild it
after every image pull.

```bash
# inside the mucoll sim container (e.g. mucoll-sim-ubuntu24:main)
source /opt/setup_mucoll.sh
simulation/calibhits/build.sh   # -> simulation/calibhits/install/<stack>/lib
```

`<stack>` is the image's DD4hep build, e.g. `dd4hep-1.38-eqfd5l6` for `:main`
or `dd4hep-1.37-5u6u23u` for v3.1 (`stack_tag.sh`). Builds for different images
coexist, and `setup.sh` loads the one matching the image it runs in.

**Which image?** The pinned `configs/MAIAConfig` (≥ `v9-66`) uses the k4Reco
calorimeter digitizers (`DigitiserSiliconMip`, …). Those exist in the `:main`
image but not in the v3.1 release. To use v3.1, check out `configs/MAIAConfig`
at `v9`. The calibration-hit plugin works with both (validated with DD4hep 1.37
and 1.38).

## Run

```bash
source /opt/setup_mucoll.sh
source setup_config.sh /path/to/mucoll-benchmarks MAIA_v0
source simulation/calibhits/setup.sh   # puts the plugin on LD_LIBRARY_PATH
ddsim --steeringFile simulation/steer_calib.py --inputFiles gen.edm4hep.root \
      --outputFile sim_calib.edm4hep.root -N 100
```

`steer_calib.py` is `steer_baseline.py` plus the action mapping. It finds
`steer_baseline.py` next to itself. If you wrap it from another directory,
`export MUCOLL_BENCHMARKS` first (`setup_config.sh` sets this variable but does
not export it).

Options are set in `steer_calib.py`:

- `StoreCalibContributions` (default `False`): keep MC-particle contributions on
  the calibration hits. The files are about 5x larger with them. The EM
  fraction comes from `*CalibEM` instead.
- `MergeCalibContributions` (default `True`): one contribution per Geant4 track
  per cell, instead of one per step.
- `ApplyBirksLaw` (default `True`): applies to the active collection only, to
  match the baseline.
- `IncludeInvisibleEnergy` (default `True`) and `EMCollectionSuffix` (default
  `CalibEM`; `""` disables it).
- `IncludePrompt` (default `True`), `PromptSuffix` (`Prompt`),
  `PromptWindowMin`/`PromptWindowMax` [ns] (set per calorimeter from
  `PROMPT_WINDOWS`), and `PromptCorrectForPropagation` (`True`, as in the
  digitizer).

## Validate

```bash
python3 simulation/calibhits/check_calibhits.py sim_calib.edm4hep.root \
        [--baseline same_seed_baseline.edm4hep.root]
```

The script prints the active and calibration energy per calorimeter and checks
that every active cell has a calibration hit. With `--baseline`, it also checks
that the active collections are identical to the standard simulation.

Results from the first test (2 × 20 GeV π⁻, v3.1 image):

- **Active collections and `MCParticles`:** identical to the baseline.
- **Calibration/active energy ratio, ECal:** about 40. This agrees with the
  digitiser's ECal MIP→GeV factor, 0.00641/0.0001575 = 40.7.
- **Contribution positions:** all within the geometric bounds of their cells,
  at most 4.5 mm (ECal) and 29 mm (HCal) from the cell centre.
- **Time per event:** unchanged.

## Production (Perlmutter; OSCAR untested)

`production/submit_calib.py` submits one job array per particle type, with
Shifter on Perlmutter or Apptainer on OSCAR (`--site`, auto-detected; see
QUICKSTART.md).
Each task runs `production/run_chain_calib.sh` (gun → sim → digi → reco, no BIB)
in the Shifter image:

- **Queue:** shared QOS, one physical core and 4 GB per task (the chain peaks
  at 1.7 GB).
- **Code snapshot:** `<output>/code_snapshot/` is a frozen copy of the whole
  checkout: configs, steering, chain, and the plugin build for the image's
  stack. `VERSIONS.txt` next to it records the commits, the MAIAConfig pin and
  the image ID. Jobs run only from the snapshot, with the image pinned by its
  ID, so later edits, rebuilds or image pulls never affect them.
- **Output:** only the reco file is kept. It contains all the sim collections,
  `*Calib*` included.

```bash
python3 production/submit_calib.py --dry-run   # edit the EDIT ME block first
```

## Plots and event displays

Run them inside the container:

- `plot_validation.py`: material breakdown, depth and lateral mapping checks,
  per-cell and per-event comparisons. It needs samples made with
  `MergeCalibContributions=False` and `StoreCalibContributions=True`. It is not
  yet in the house style.
- `event_display.py reco.edm4hep.root -e 0 1 -o dir/` (house style, via
  `house_style.py`): per event, an r–z
  overview plus zoomed r–z and x–y maps of calibration, active, invisible and
  EM-fraction energy, and the Pandora cluster assignment of each cell. The true
  trajectory is drawn from the primary's own tracker SimHits.

### Closure with respect to the true energy

`plot_closure.py --sample e=DIR --sample pi=DIR -o closure/` plots
E_calib / E_true vs E_true and θ, with and without invisible energy, and E_calib
vs E_true.

For samples made with `--energy-balance`, it also plots the full budget vs θ:
every piece of energy the particle left in each top-level detector volume
(`EnergyBalanceCheckAction`, CSV per job). In tests, the calorimeter envelopes
contain no energy outside the calibration cells (difference < 10⁻⁵ GeV). The
non-closure of E_calib is therefore entirely material in front of the
calorimeters (solenoid, trackers, beam pipe/nozzles) plus leakage.

Use a separate production with `--energy-balance --sim-only`; its sim files and
CSVs belong together. The events differ from the main production's even with the
same seeds unless `--nevents` is also the same: the gun generator draws all
momenta first, then all angles.

### Animated shower development

`animate_shower.py sim.edm4hep.root -e 0 -o anim/` makes GIFs of the calibration
energy accumulating in time (r–z and x–y), plus key frames as pdf/png. The time
axis is linear through the prompt shower, then logarithmic to 100 µs, with a
final frame for everything later.

It needs per-step contributions (`StoreCalibContributions=True`,
`MergeCalibContributions=False`), which is too large for production: use it on a
few events.

## Using it for cluster calibration

Per event, build a map from cellID to calibration energy from the four `*Calib`
collections. The digitised `CalorimeterHit`s in a Pandora cluster carry the
same cellIDs. Then:

- E_dep^clus = sum of the map over the cluster's cells (each hit is in exactly
  one cluster, so w_geo = 1);
- R = E_EM^clus / E_dep^clus.

## Caveats and to-dos

- **Many cells have absorber energy but no active hit.** These are soft
  particles stopping in the W/steel, with nothing reaching the 0.5 mm Si or
  3 mm scintillator. In the π⁻ test this was roughly half of the calibration
  energy. Such cells are never digitised, so a strict "clustered cells only"
  E_dep^clus leaves them out. With 5 mm ECal cells, that makes the target depend
  strongly on segmentation, much more than with ATLAS cell sizes. Decide whether
  to include unhit cells inside the cluster footprint.
- **Energy outside the layer slices** is not counted: stave or module
  envelopes, the solenoid in front of the ECal, and so on.
- **Long term**, the action belongs in k4geo `plugins/`, next to the other
  custom sensitive actions there, so that it ships with the software images.
