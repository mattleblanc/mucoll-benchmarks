#!/bin/bash
# Particle-gun chain with calorimeter calibration hits:
#   GEN (flat in p and theta) -> SIM (steer_calib.py) -> DIGI -> RECO   (no BIB)
# Runs INSIDE the mucoll sim container. The reco output keeps every collection
# (no BIB -> "keep *"), including the *Calib, *CalibInvisible and *CalibEM hits,
# so only the reco file is kept.
#
#   run_chain_calib.sh --job-id N --nevents N --outdir DIR --pdg "211 -211" \
#       --pmin 1 --pmax 50 --theta-min 10 --theta-max 170 [--seed S] [--workdir DIR] \
#       [--benchmarks DIR]   mucoll-benchmarks checkout (default: the one containing this script)
#       [--steering FILE] [--plugin-dir DIR]   (defaults: from the checkout) \
#       [--energy-balance]   also write energy_balance_<ID>.csv (energy per subdetector)
#       [--sim-only]         stop after SIM and keep sim_<ID>.edm4hep.root instead of reco
#
# Batch jobs run a copy of this script from the production's code_snapshot/.
# Everything is inside main() so that bash parses the whole file before running
# anything: editing a script while jobs execute it would otherwise break them.

main() {
    set -euo pipefail

    local JOB_ID="" NEVENTS=100 OUTDIR="" PDG="211" PMIN=1 PMAX=50
    local THETA_MIN=10 THETA_MAX=170 SEED="" WORKBASE="${PSCRATCH:-/tmp}/calibhits_work"
    local STEERING="" PLUGIN_DIR="" BENCH="" EBAL=0 SIM_ONLY=0
    while [ $# -gt 0 ]; do
        case "$1" in
            --job-id)     JOB_ID="$2";     shift 2 ;;
            --nevents)    NEVENTS="$2";    shift 2 ;;
            --outdir)     OUTDIR="$2";     shift 2 ;;
            --pdg)        PDG="$2";        shift 2 ;;
            --pmin)       PMIN="$2";       shift 2 ;;
            --pmax)       PMAX="$2";       shift 2 ;;
            --theta-min)  THETA_MIN="$2";  shift 2 ;;
            --theta-max)  THETA_MAX="$2";  shift 2 ;;
            --seed)       SEED="$2";       shift 2 ;;
            --workdir)    WORKBASE="$2";   shift 2 ;;
            --benchmarks) BENCH="$2";      shift 2 ;;
            --steering)   STEERING="$2";   shift 2 ;;
            --plugin-dir) PLUGIN_DIR="$2"; shift 2 ;;
            --energy-balance) EBAL=1;      shift 1 ;;
            --sim-only)   SIM_ONLY=1;      shift 1 ;;
            *) echo "Unknown argument: $1" >&2; return 1 ;;
        esac
    done
    if [ -z "$JOB_ID" ] || [ -z "$OUTDIR" ]; then
        echo "ERROR: --job-id and --outdir are required" >&2; return 1
    fi
    SEED="${SEED:-$((100000 + JOB_ID))}"

    local HERE
    HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    export MUCOLL_BENCHMARKS="${BENCH:-$(cd "$HERE/../../.." && pwd)}"
    STEERING="${STEERING:-$MUCOLL_BENCHMARKS/simulation/steer_calib.py}"
    if [ -n "$PLUGIN_DIR" ]; then export CALIBHITS_INSTALL_DIR="$PLUGIN_DIR"; fi
    # setup.sh next to this script (code snapshot) if present, else the checkout's
    local PLUGIN_SETUP="$HERE/setup.sh"
    [ -f "$PLUGIN_SETUP" ] || PLUGIN_SETUP="$MUCOLL_BENCHMARKS/simulation/calibhits/setup.sh"

    set +u
    source /opt/setup_mucoll.sh
    set +o pipefail   # setup_config.sh uses `find | head` (SIGPIPE under pipefail)
    source "$MUCOLL_BENCHMARKS/setup_config.sh" "$MUCOLL_BENCHMARKS" MAIA_v0
    set -o pipefail
    source "$PLUGIN_SETUP"
    set -u

    # Scratch on $PSCRATCH, not /tmp: /tmp on Perlmutter compute nodes is RAM-backed
    # and counts against the job's memory request.
    WORKDIR="$WORKBASE/job_${JOB_ID}_${SLURM_JOB_ID:-local}_$RANDOM"
    mkdir -p "$WORKDIR"
    cd "$WORKDIR"
    trap 'cd /; rm -rf "$WORKDIR"' EXIT
    cp -r "$MUCOLL_CONFIG/$MUCOLL_CONFIG_NAME/PandoraSettings/" ./

    echo "=== calib chain: job $JOB_ID, $NEVENTS events, PDG=[$PDG], p=[$PMIN,$PMAX] GeV," \
         "theta=[$THETA_MIN,$THETA_MAX] deg, seed=$SEED, host=$(hostname)"
    echo "=== code $MUCOLL_BENCHMARKS"
    echo "=== steering $STEERING, plugin ${LD_LIBRARY_PATH%%:*}"

    local t0 t1 t2 t3 t4
    t0=$SECONDS
    # shellcheck disable=SC2086
    python3 "$MUCOLL_BENCHMARKS/generation/pgun/pgun_edm4hep.py" \
        -p 1 -e "$NEVENTS" -s "$SEED" --pdg $PDG --p "$PMIN" "$PMAX" \
        --theta "$THETA_MIN" "$THETA_MAX" -- gen.edm4hep.root
    t1=$SECONDS
    if [ "$EBAL" -eq 1 ]; then export CALIBHITS_ENERGY_BALANCE_CSV="$WORKDIR/energy_balance.csv"; fi
    ddsim --steeringFile "$STEERING" \
        --numberOfEvents "$NEVENTS" --inputFiles gen.edm4hep.root \
        --outputFile sim.edm4hep.root --random.seed "$SEED"
    t2=$SECONDS
    mkdir -p "$OUTDIR"
    if [ "$EBAL" -eq 1 ]; then cp energy_balance.csv "$OUTDIR/energy_balance_${JOB_ID}.csv"; fi
    if [ "$SIM_ONLY" -eq 1 ]; then
        cp sim.edm4hep.root "$OUTDIR/sim_${JOB_ID}.edm4hep.root.part"
        mv "$OUTDIR/sim_${JOB_ID}.edm4hep.root.part" "$OUTDIR/sim_${JOB_ID}.edm4hep.root"
        echo "=== timing [s]: gen $((t1-t0)) sim $((t2-t1)); sim-only"
        echo "=== job $JOB_ID done"
        return 0
    fi
    k4run "$MUCOLL_CONFIG/$MUCOLL_CONFIG_NAME/digi_steer.py" -n -1 \
        --inputFiles sim.edm4hep.root --outputFile digi.edm4hep.root --RandSeed "$SEED"
    t3=$SECONDS
    k4run "$MUCOLL_CONFIG/$MUCOLL_CONFIG_NAME/reco_steer.py" -n -1 \
        --inputFiles digi.edm4hep.root --outputFile reco.edm4hep.root
    t4=$SECONDS

    cp reco.edm4hep.root "$OUTDIR/reco_${JOB_ID}.edm4hep.root.part"
    mv "$OUTDIR/reco_${JOB_ID}.edm4hep.root.part" "$OUTDIR/reco_${JOB_ID}.edm4hep.root"
    echo "=== timing [s]: gen $((t1-t0)) sim $((t2-t1)) digi $((t3-t2)) reco $((t4-t3));" \
         "output $(du -h --apparent-size "$OUTDIR/reco_${JOB_ID}.edm4hep.root" | cut -f1)"
    echo "=== job $JOB_ID done"
}

main "$@"
exit $?
