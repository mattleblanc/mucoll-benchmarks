# Calibration hits on Perlmutter: step by step

This is how to produce MAIA particle-gun samples with calorimeter calibration
hits on Perlmutter, and how to make the validation plots. For what the
calibration hits are and how they work, see [README.md](README.md).

Everything runs in the `mucoll-sim-ubuntu24:main` container image. That image
is needed because the pinned `configs/MAIAConfig` uses calorimeter digitizers
that the v3.1 release doesn't have yet. Commands are run from a Perlmutter
login node unless a step says otherwise. **On OSCAR, follow the same steps with
the substitutions in [Running on OSCAR](#running-on-oscar) at the end.**

## 1. Get the code (once)

```bash
cd $HOME   # or wherever you keep code; it must be visible from compute nodes
git clone https://github.com/mattleblanc/mucoll-benchmarks.git
cd mucoll-benchmarks
git checkout calibhits
git submodule update --init
export BENCH=$PWD
```

## 2. Pull the image (once, and whenever you want a newer `main`)

```bash
shifterimg pull ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:main
```

On Perlmutter, pulled images are shared by all users. Re-pulling `:main`
updates it for everyone, and everyone then has to rebuild the plugin (step 3)
before running by hand. Submitted productions are unaffected because they pin
the image ID. Check with the group before re-pulling.

The pull takes 10–20 minutes. `main` is a moving tag: after every new pull you
must rebuild the plugin (step 3). Productions pin the exact image they ran
with, so a new pull never affects jobs that are already submitted.

## 3. Build the plugin inside the image (after every image pull)

```bash
shifter --image=ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:main bash -c \
    "source /opt/setup_mucoll.sh && $BENCH/simulation/calibhits/build.sh"
```

This builds into `simulation/calibhits/install/<stack>/`, where `<stack>`
names the image's DD4hep build. Builds for different images live side by side,
and `setup.sh` always picks the one for the image you are running in.

## 4. Quick test (a few minutes, on the login node)

Runs generation → simulation → digitization → reconstruction for 3 π± events:

```bash
shifter --image=ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:main bash \
    $BENCH/simulation/calibhits/production/run_chain_calib.sh \
    --job-id 0 --nevents 3 --pdg "211 -211" --pmin 1 --pmax 50 \
    --outdir $PSCRATCH/calibtest --workdir $PSCRATCH/calibtest/work
```

It should finish with `=== job 0 done` and leave
`$PSCRATCH/calibtest/reco_0.edm4hep.root`. Check the calibration hits:

```bash
shifter --image=ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:main bash -c \
    "source /opt/setup_mucoll.sh && python3 $BENCH/simulation/calibhits/check_calibhits.py \
     $PSCRATCH/calibtest/reco_0.edm4hep.root"
```

It prints the active and calibration energy per calorimeter and ends with
`RESULT: OK`.

## 5. Production (batch, shared queue)

Edit the `EDIT ME` block at the top of
`simulation/calibhits/production/submit_calib.py`:

- `SAMPLES`, `NJOBS`, `NEVENTS_PER_JOB`, and the momentum/angle ranges;
- the output directory: by default `/global/cfs/cdirs/m5197/$USER/calibhits/v1`.
  Every production needs a **new** directory: pass `--output-base …/v2` and so
  on, or edit `output_base` in `SITES["perlmutter"]`;
- the account in `SITES["perlmutter"]["sbatch"]` (`m5197`), if you run under a
  different NERSC project.

Then:

```bash
cd $BENCH/simulation/calibhits/production
python3 submit_calib.py --dry-run   # checks image + plugin, writes the batch scripts
python3 submit_calib.py             # submits one job array per sample
squeue --me                               # monitor
```

What happens:

- **Resources:** each task runs one chain on one core with 4 GB in the
  `shared` queue. 250 events take about 1.5 h.
- **Outputs:** `<OUTPUT_BASE>/<sample>/reco_<N>.edm4hep.root`. These contain
  all the simulation collections, the calibration hits included
  (`*Calib`, `*CalibInvisible`, `*CalibEM` and their `*Prompt` versions).
- **Logs:** live logs go to `$PSCRATCH/calibhits_logs/`. Gzipped copies land in
  `<sample>/logs/` when a task ends.
- **Code snapshot:** `<OUTPUT_BASE>/code_snapshot/` holds a frozen copy of the
  code the jobs actually run, and `VERSIONS.txt` records the commits and the
  image ID. So you can keep editing your checkout while jobs run.
- **Resubmitting:** running the submitter again for the same `OUTPUT_BASE`
  refuses, because the snapshot already exists. To redo only failed tasks, use
  `sbatch <sample>/logs/submit.sbatch`: tasks whose output already exists exit
  immediately.

For the energy-budget closure plots, make a small extra production that keeps
the simulation output and writes the per-subdetector energy CSVs:

```bash
python3 submit_calib.py --output-base /global/cfs/cdirs/m5197/<you>/calibhits/v1-ebal \
    --njobs 10 --nevents 200 --energy-balance --sim-only
```

## 6. Plots

Run these inside the image. Each plot is written as a pdf + png.

```bash
IMG=ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:main
CH=$BENCH/simulation/calibhits
P=/global/cfs/cdirs/m5197/<you>/calibhits

# closure with respect to the true energy; budget plots need the --energy-balance sample
shifter --image=$IMG bash -c "source /opt/setup_mucoll.sh && python3 $CH/plot_closure.py \
    --sample e=$P/v1/e_p1-50_theta10-170 --sample pi=$P/v1/pi_p1-50_theta10-170 -o closure/"

# event displays (calibration, active, invisible, EM fraction, ratio, clusters)
shifter --image=$IMG bash -c "source /opt/setup_mucoll.sh && python3 $CH/event_display.py \
    $P/v1/pi_p1-50_theta10-170/reco_0.edm4hep.root -e 0 1 2 -o displays/"
```

`displays/CAPTIONS.md` explains how to read each event display.

The animated shower displays (`animate_shower.py`) need every Geant4 step kept,
which is far too large for production. Make a few events by hand with
`StoreCalibContributions=True` and `MergeCalibContributions=False` (see
README.md).

## Things that commonly go wrong

| Symptom | Cause / fix |
|---|---|
| `ERROR: no CalibHits build for this software stack` | No plugin build for this image yet (or the image was re-pulled): redo step 3 |
| `cannot import name 'DigitiserSiliconMip'` | Running with the v3.1 image: use `:main`, or check out `configs/MAIAConfig` at `v9` for v3.1 |
| `... code_snapshot exists` when submitting | Each production needs a new `OUTPUT_BASE` |
| A job fails right at the end with `unexpected EOF` | A script was edited while running. Productions run from the snapshot, so this only affects hand-run chains |

## Running on OSCAR

> **Not yet tested on OSCAR.** The OSCAR mode was written from the
> `mucoll-slurm` conventions but has not been run there. Please report what
> breaks.

The code is the same. Only the container commands and the batch settings
differ. OSCAR uses Apptainer with an image *file* (SIF) instead of Shifter.

**Setup.** Set this once per shell, e.g. in `~/.bashrc`, as in the
`mucoll-slurm` README:

```bash
export APPTAINER_TMPDIR=/oscar/scratch/$USER/apptainer_tmp
export APPTAINER_CACHEDIR=/oscar/scratch/$USER/apptainer_cache
mkdir -p $APPTAINER_TMPDIR $APPTAINER_CACHEDIR
```

**Step 2 (image).** Pull the image to a *dated* SIF file. A SIF file never
changes, so the date in its name records exactly which `main` build you have:

```bash
SIF=/oscar/data/mleblan6/mucoll/mucoll-sim-ubuntu24-main-$(date +%Y%m%d).sif
apptainer pull $SIF docker://ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:main
```

**Steps 3, 4 and 6 (build, test, plots).** Replace every
`shifter --image=ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:main` with
`apptainer exec --bind /oscar/data/mleblan6/mucoll,$BENCH $SIF`. Apptainer
only sees the directories you bind, so bind the code checkout and every input
and output area. For example, step 3:

```bash
apptainer exec --bind /oscar/data/mleblan6/mucoll,$BENCH $SIF bash -c     "source /opt/setup_mucoll.sh && $BENCH/simulation/calibhits/build.sh"
```

In step 4, use an `--outdir` under `/oscar/data/...` or your home directory
instead of `$PSCRATCH`, and drop `--workdir`. The chain then works in `/tmp`
on the node, as the `mucoll-slurm` chains do.

**Step 5 (production).** Same submitter, OSCAR mode. Check `SITES["oscar"]`
in the `EDIT ME` block (output area, and `--partition`/`--account` if your
allocation needs them), then:

```bash
cd $BENCH/simulation/calibhits/production
python3 submit_calib.py --site oscar --image $SIF --dry-run
python3 submit_calib.py --site oscar --image $SIF
```

The site is detected automatically, so `--site oscar` is optional on OSCAR.
`--image` is always required there. The outputs, code snapshot, `VERSIONS.txt`
and logs work as on Perlmutter. Live logs go to
`/oscar/scratch/$USER/calibhits_logs/`.
