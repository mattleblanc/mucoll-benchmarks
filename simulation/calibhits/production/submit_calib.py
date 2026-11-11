#!/usr/bin/env python3
"""Submit calibration-hit particle-gun samples (SLURM job arrays) on Perlmutter or OSCAR.

Each sample is flat in momentum and polar angle (no BIB) and runs
run_chain_calib.sh (GEN -> SIM with calibration hits -> DIGI -> RECO) inside
the mucoll sim container. One SLURM job array per sample.

    python3 submit_calib.py --dry-run     # check image + plugin, write the batch scripts
    python3 submit_calib.py               # submit all samples
    python3 submit_calib.py --samples pi  # only some samples
    python3 submit_calib.py --site oscar --image /path/to/mucoll-sim-ubuntu24-main-YYYYMMDD.sif

Sites (auto-detected, or --site):
  perlmutter  Shifter image (pinned by its image ID), shared QOS, one core + 4 GB.
  oscar       Apptainer SIF file, default partition, one core + 4 GB.
              NOT YET TESTED on OSCAR (written from the mucoll-slurm conventions).

Reproducibility: every production gets <output-base>/code_snapshot/, a frozen copy
of this mucoll-benchmarks checkout (including configs/MAIAConfig and the plugin
build for the image's software stack) plus VERSIONS.txt (commits, image).
Jobs only run code from the snapshot and a pinned image (Shifter image ID / a SIF
file), so later edits, rebuilds, submodule updates or image pulls never affect them.

Run from a login node (sbatch is not available inside the container). The plugin
must have been built inside the same image (simulation/calibhits/build.sh).
"""
import argparse
import datetime
import os
import shlex
import shutil
import socket
import subprocess
import sys

# =============================== EDIT ME ====================================
SAMPLES = {
    # tag : (PDG list, seed offset)
    "e":  ("11",       100000),
    "pi": ("211 -211", 200000),
}
NJOBS = 40               # array size per sample
NEVENTS_PER_JOB = 250    # -> 10k events per sample
PMIN, PMAX = 1.0, 50.0   # GeV, flat
THETA_MIN, THETA_MAX = 10.0, 170.0  # deg, flat
TIME = "03:00:00"
MEM = "4G"               # the chain peaks at ~1.7 GB

USER = os.environ.get("USER", "user")
SITES = {
    "perlmutter": {
        # per-user default; use a new directory for every production (or --output-base)
        "output_base": f"/global/cfs/cdirs/m5197/{USER}/calibhits/v1",
        "image": "ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:main",   # Shifter image name
        # one physical core = 2 logical CPUs; 4 GB stays within one core's memory
        # share, so the shared QOS charges exactly one core
        "sbatch": ["--account=m5197", "--constraint=cpu", "--qos=shared", "--cpus-per-task=2"],
        "log_base": os.path.join(os.environ.get("PSCRATCH", "/tmp"), "calibhits_logs"),
    },
    "oscar": {
        "output_base": f"/oscar/data/mleblan6/mucoll/{USER}/calibhits/v1",
        "image": None,   # path to an Apptainer SIF file: pass --image (dated filename recommended)
        "sbatch": ["--nodes=1", "--cpus-per-task=1"],   # add e.g. "--partition=batch" if needed
        "log_base": f"/oscar/scratch/{USER}/calibhits_logs",
    },
}
# ============================================================================

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
CALIB = os.path.join(BENCH, "simulation", "calibhits")


def run(cmd, **kw):
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True, **kw)


def detect_site():
    if os.environ.get("NERSC_HOST") == "perlmutter":
        return "perlmutter"
    if os.path.isdir("/oscar") or "oscar" in socket.getfqdn():
        return "oscar"
    raise SystemExit("cannot tell the site: pass --site perlmutter or --site oscar")


# --------------------------------------------------------------------------
# Site-specific container handling
class Perlmutter:
    """Shifter: the image is pinned by its ID in #SBATCH --image."""

    def __init__(self, image):
        res = run(["shifterimg", "lookup", image])
        self.image, self.pin = image, res.stdout.strip()
        if res.returncode != 0 or not self.pin:
            raise SystemExit(f"image {image} not pulled? (shifterimg pull {image})\n{res.stderr}")

    def exec_in_image(self, bash_cmd):
        return run(["shifter", f"--image={self.image}", "bash", "-c", bash_cmd])

    def sbatch_lines(self):
        return [f"--image=id:{self.pin}"]

    def wrap(self, chain_cmd, binds):
        return f"shifter {chain_cmd}"

    def describe(self):
        return f"image:           {self.image}\nimage id:        {self.pin}\n"


class Oscar:
    """Apptainer: the SIF file itself is the pinned image."""

    def __init__(self, image):
        if not image or not os.path.isfile(image):
            raise SystemExit("on OSCAR pass --image /path/to/<image>.sif (apptainer pull ... docker://"
                             "ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:main)")
        self.image = os.path.abspath(image)
        st = os.stat(self.image)
        self.stamp = (f"{st.st_size} bytes, modified "
                      f"{datetime.datetime.fromtimestamp(st.st_mtime).isoformat(timespec='seconds')}")

    def exec_in_image(self, bash_cmd):
        return run(["apptainer", "exec", "--bind", BENCH, self.image, "bash", "-c", bash_cmd])

    def sbatch_lines(self):
        return []

    def wrap(self, chain_cmd, binds):
        return f"apptainer exec --bind {','.join(binds)} {shlex.quote(self.image)} {chain_cmd}"

    def describe(self):
        return f"image (SIF):     {self.image}\nimage file:      {self.stamp}\n"


def stack_tag(container):
    """Software-stack tag of the image (selects the matching plugin build)."""
    cmd = (f"source /opt/setup_mucoll.sh >/dev/null 2>&1; source {CALIB}/stack_tag.sh; "
           "calibhits_stack_tag")
    res = container.exec_in_image(cmd)
    tag = res.stdout.strip().splitlines()[-1] if res.stdout.strip() else ""
    if res.returncode != 0 or not tag:
        raise SystemExit(f"could not determine the software stack of the image:\n{res.stderr}")
    return tag


def git(*a):
    return run(["git", *a], cwd=BENCH).stdout.strip()


def make_snapshot(snap, tag, container, site):
    """Frozen copy of the checkout: no .git, no builds except the one for this stack."""
    if os.path.exists(snap):
        raise SystemExit(f"{snap} exists: choose a new --output-base for a new production")

    def ignore(d, names):
        skip = {n for n in names if n in (".git", "__pycache__")}
        rel = os.path.relpath(d, BENCH)
        if rel == os.path.join("simulation", "calibhits"):
            skip |= {"build"} & set(names)
        if rel == os.path.join("simulation", "calibhits", "install"):
            skip |= {n for n in names if n != tag}
        return skip
    shutil.copytree(BENCH, os.path.join(snap, "mucoll-benchmarks"), ignore=ignore, symlinks=True)
    with open(os.path.join(snap, "VERSIONS.txt"), "w") as fh:
        fh.write(f"date:            {datetime.datetime.now().isoformat(timespec='seconds')}\n")
        fh.write(f"site:            {site}\n")
        fh.write(f"command:         {' '.join(sys.argv)}\n")
        fh.write(container.describe())
        fh.write(f"software stack:  {tag}\n")
        fh.write(f"mucoll-benchmarks: {git('rev-parse', 'HEAD')} ({git('rev-parse', '--abbrev-ref', 'HEAD')})\n")
        fh.write("submodules:\n" + run(["git", "submodule", "status"], cwd=BENCH).stdout.rstrip() + "\n")
        fh.write(f"uncommitted changes (included in the snapshot):\n{git('status', '--short')}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--site", choices=list(SITES), help="default: auto-detect")
    ap.add_argument("--samples", nargs="+", default=list(SAMPLES), choices=list(SAMPLES))
    ap.add_argument("--njobs", type=int, default=NJOBS)
    ap.add_argument("--nevents", type=int, default=NEVENTS_PER_JOB)
    ap.add_argument("--image", help="Perlmutter: Shifter image name; OSCAR: SIF file (required)")
    ap.add_argument("--output-base", help="production directory (one code snapshot each)")
    ap.add_argument("--energy-balance", action="store_true",
                    help="also write per-event energy per subdetector (energy_balance_<ID>.csv)")
    ap.add_argument("--sim-only", action="store_true", help="stop after SIM, keep sim_<ID>.edm4hep.root")
    args = ap.parse_args()

    site = args.site or detect_site()
    cfg = SITES[site]
    base = os.path.abspath(args.output_base or cfg["output_base"])
    image = args.image or cfg["image"]
    container = (Perlmutter if site == "perlmutter" else Oscar)(image)
    chain_extra = (" --energy-balance" if args.energy_balance else "") + (" --sim-only" if args.sim_only else "")
    final = "sim" if args.sim_only else "reco"

    tag = stack_tag(container)
    if not os.path.isfile(os.path.join(CALIB, "install", tag, "lib", "libCalibHits.components")):
        raise SystemExit(f"no plugin build for this image (stack {tag}): run simulation/calibhits/build.sh "
                         f"inside that image first")
    print(f"site {site}, image {container.image}, stack {tag}")

    snap = os.path.join(base, "code_snapshot")
    if not args.dry_run:
        make_snapshot(snap, tag, container, site)
        print(f"code snapshot -> {snap} (see VERSIONS.txt)")
    chain = os.path.join(snap, "mucoll-benchmarks", "simulation", "calibhits", "production", "run_chain_calib.sh")

    for tag_s in args.samples:
        pdg, seed0 = SAMPLES[tag_s]
        outdir = os.path.join(base, f"{tag_s}_p{PMIN:g}-{PMAX:g}_theta{THETA_MIN:g}-{THETA_MAX:g}")
        logdir = os.path.join(outdir, "logs")
        os.makedirs(logdir, exist_ok=True)
        # Live logs on scratch (many small files are wasteful on project/data file
        # systems); one gzipped log per task is copied to <outdir>/logs/ at the end.
        tmplog = os.path.join(cfg["log_base"], os.path.basename(base.rstrip("/")), os.path.basename(outdir))
        os.makedirs(tmplog, exist_ok=True)
        chain_cmd = (f"bash {chain} --job-id \"$ID\" --nevents {args.nevents} --outdir {outdir} \\\n"
                     f"    --pdg \"{pdg}\" --pmin {PMIN:g} --pmax {PMAX:g} \\\n"
                     f"    --theta-min {THETA_MIN:g} --theta-max {THETA_MAX:g} --seed $(({seed0} + ID)){chain_extra}")
        directives = [f"--job-name=calib_{tag_s}", f"--time={TIME}", "--ntasks=1", f"--mem={MEM}",
                      f"--array=0-{args.njobs - 1}", f"--output={tmplog}/job_%a.log"]
        directives += cfg["sbatch"] + container.sbatch_lines()
        script = "#!/bin/bash\n" + "".join(f"#SBATCH {d}\n" for d in directives) + f"""
# site {site}; image {container.image}
ID=$SLURM_ARRAY_TASK_ID
if [ -f "{outdir}/{final}_$ID.edm4hep.root" ]; then echo "output exists, skipping"; exit 0; fi
{container.wrap(chain_cmd, [base])}
rc=$?
gzip -c {tmplog}/job_$ID.log > {logdir}/job_$ID.log.gz
exit $rc
"""
        # A dry run must never overwrite the record of a real submission
        path = os.path.join(logdir, "submit.sbatch" + (".dryrun" if args.dry_run else ""))
        with open(path, "w") as fh:
            fh.write(script)
        print(f"[{tag_s}] PDG [{pdg}], {args.njobs} x {args.nevents} events -> {outdir}")
        if args.dry_run:
            print(f"    wrote {path} (dry run)")
            continue
        res = run(["sbatch", path])
        print("    " + (res.stdout.strip() or res.stderr.strip()))


if __name__ == "__main__":
    main()
