#!/usr/bin/env python3
"""Closure of the calibration hits with respect to the true particle energy.

    python3 plot_closure.py --sample e=DIR_e --sample pi=DIR_pi -o closure/

Each DIR holds the outputs of run_chain_calib.sh: sim_<ID>.edm4hep.root or
reco_<ID>.edm4hep.root, and optionally energy_balance_<ID>.csv (--energy-balance).
Writes one pdf + png per plot (house style):

  closure_vs_E, closure_vs_theta  median of E_calib / E_true per bin (bars: statistical
                                  uncertainty of the median), with and without invisible energy
  closure_dist                    distribution of E_calib / E_true
  calib_vs_true_<s>               E_calib vs E_true per event
  budget_vs_theta_<s>             (needs the CSVs) where the energy goes, as a fraction of
                                  E_true: calorimeter cells (= the calibration hits),
                                  calorimeter structure outside the layers, solenoid,
                                  trackers, beam pipe + nozzles, other, escaped
  bookkeeping_closure             (needs the CSVs) sum of all components / E_true - 1

With the prompt collections (*CollectionCalib*Prompt), also:
  closure_* and budget_* get the prompt calibration energy as an extra series
  calib_vs_true_prompt_<s>        E_calib(prompt) vs E_true per event
  late_vs_E, late_vs_theta        late fraction 1 - E_prompt / E_all (energy-weighted),
                                  ECal and HCal, per sample
  late_components_<s>             late fraction per component: EM visible, non-EM
                                  visible, invisible

E_calib is summed over the four calorimeters (*CollectionCalib).
Run inside the mucoll container (uproot, numpy, matplotlib).
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import house_style as hs  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

CALOS = ["ECalBarrel", "ECalEndcap", "HCalBarrel", "HCalEndcap"]
LABEL = {"e": "e⁻", "pi": "π±", "gamma": "γ", "n": "n", "K0L": "K⁰_L"}
SAMPLE_COLOR = {"e": hs.CAT[0], "pi": hs.CAT[1]}
EXPERIMENT = "MAIA simulation"


# --------------------------------------------------------------------------
def _sum(tree, coll):
    key = f"{coll}/{coll}.energy"
    if key not in tree:
        return None
    return np.array([x.sum() for x in tree[key].array(library="np")])


def load_dir(d):
    """Per-event arrays for all files of one sample."""
    files = sorted(glob.glob(os.path.join(d, "sim_*.edm4hep.root"))) or \
        sorted(glob.glob(os.path.join(d, "reco_*.edm4hep.root")))
    if not files:
        raise SystemExit(f"no sim_/reco_ files in {d}")
    cols = defaultdict(list)
    for f in files:
        job = int(re.search(r"_(\d+)\.edm4hep\.root$", f).group(1))
        t = uproot.open(f)["events"]
        mc = {k: t[f"MCParticles/MCParticles.{k}"].array(library="np")
              for k in ("generatorStatus", "momentum.x", "momentum.y", "momentum.z", "mass")}
        n = len(mc["mass"])
        etrue, theta = np.zeros(n), np.zeros(n)
        for i in range(n):
            j = int(np.flatnonzero(mc["generatorStatus"][i] == 1)[0])
            px, py, pz, m = (mc[k][i][j] for k in ("momentum.x", "momentum.y", "momentum.z", "mass"))
            p = np.sqrt(px * px + py * py + pz * pz)
            etrue[i] = np.sqrt(p * p + m * m)
            theta[i] = np.degrees(np.arccos(pz / p))
        tot = defaultdict(lambda: np.zeros(n))
        per_calo = {}
        for c in CALOS:
            grp = c[:4]   # ECal / HCal
            for k, suf in (("calib", "Calib"), ("inv", "CalibInvisible"), ("em", "CalibEM"), ("active", ""),
                           ("calibP", "CalibPrompt"), ("invP", "CalibInvisiblePrompt"), ("emP", "CalibEMPrompt")):
                v = _sum(t, f"{c}Collection{suf}")
                if v is not None:
                    tot[k] = tot[k] + v
                    tot[f"{k}_{grp}"] = tot[f"{k}_{grp}"] + v
                    if k == "calib":
                        per_calo[c] = v
        cols["job"].append(np.full(n, job))
        cols["event"].append(np.arange(n))
        cols["etrue"].append(etrue)
        cols["theta"].append(theta)
        for k, v in tot.items():
            cols[k].append(v)
        for c, v in per_calo.items():
            cols[f"calib_{c}"].append(v)
    out = {k: np.concatenate(v) for k, v in cols.items()}
    out["has_prompt"] = "calibP" in out
    out["budget"] = load_budget(d)
    return out


REGION_GROUPS = [  # (label, regex on the subdetector name)
    ("Solenoid", r"^Solenoid"),
    ("Trackers", r"^(Vertex|InnerTrackers|OuterTrackers)"),
    ("Beam pipe, nozzles", r"^(Beampipe|Nozzle)"),
]


def load_budget(d):
    """(job, event) -> {region: visible + invisible}, plus E0 and escaped."""
    files = sorted(glob.glob(os.path.join(d, "energy_balance_*.csv")))
    if not files:
        return None
    b = defaultdict(dict)
    for f in files:
        job = int(re.search(r"_(\d+)\.csv$", f).group(1))
        with open(f) as fh:
            for row in csv.DictReader(fh):
                b[(job, int(row["event"]))][row["region"]] = float(row["visible_GeV"]) + float(row["invisible_GeV"])
    return b


# --------------------------------------------------------------------------
def median_and_err(x):
    if len(x) < 5:
        return np.nan, np.nan
    q16, q50, q84 = np.percentile(x, [16, 50, 84])
    return q50, 1.2533 * 0.5 * (q84 - q16) / np.sqrt(len(x))


def profile(xv, yv, edges):
    c = 0.5 * (edges[1:] + edges[:-1])
    m, e = np.full(len(c), np.nan), np.full(len(c), np.nan)
    idx = np.digitize(xv, edges) - 1
    for i in range(len(c)):
        m[i], e[i] = median_and_err(yv[idx == i])
    return c, m, e


def closure_profile(samples, xkey, edges, xlabel, name, outdir, lines):
    fig, ax = hs.squarefig()
    ymin, ymax = 1e9, -1e9
    for s, d in samples.items():
        x = d[xkey]
        series = [(d["calib"], "o", "full", "all times")]
        if d["has_prompt"]:
            series.append((d["calibP"], "s", "full", "prompt"))
        series.append((d["calib"] - d["inv"], "o", "none", "visible, all times"))
        for num, marker, fill, lab in series:
            c, m, e = profile(x, num / d["etrue"], edges)
            hs.mc_points(ax, c, m, e, color=SAMPLE_COLOR[s], marker=marker, fill=fill,
                         label=f"{LABEL[s]}, {lab}")
            ymin, ymax = min(ymin, np.nanmin(m - e)), max(ymax, np.nanmax(m + e))
    ax.axhline(1, color=hs.INK, lw=0.8, ls="--", zorder=0)
    span = max(ymax, 1.0) - ymin
    ax.set_ylim(ymin - 0.05 * span, max(ymax, 1.0) + 1.1 * span)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Median E_calib / E_true")
    hs.labels(ax, EXPERIMENT, lines)
    hs.legend(ax)
    hs.save(fig, os.path.join(outdir, name))


def closure_dist(samples, outdir):
    fig, ax = hs.squarefig()
    bins = np.linspace(0.5, 1.1, 61)
    top = 0
    for s, d in samples.items():
        series = [(d["calib"] / d["etrue"], "-", "all times")]
        if d["has_prompt"]:
            series.append((d["calibP"] / d["etrue"], (0, (1, 1)), "prompt"))
        series.append(((d["calib"] - d["inv"]) / d["etrue"], "--", "visible, all times"))
        for y, ls, lab in series:
            h, _ = np.histogram(np.clip(y, bins[0], bins[-1] - 1e-9), bins=bins)
            h = h / max(h.sum(), 1)
            ax.stairs(h, bins, color=SAMPLE_COLOR[s], ls=ls, lw=1.6, label=f"{LABEL[s]}, {lab}")
            top = max(top, h.max())
    ax.set_yscale("log")
    ax.set_ylim(1e-4, top * 300)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel("E_calib / E_true")
    ax.set_ylabel("Fraction of events")
    hs.labels(ax, EXPERIMENT, ["1 < p < 50 GeV, 10° < θ < 170°", "under/overflow in edge bins"])
    hs.legend(ax)
    hs.save(fig, os.path.join(outdir, "closure_dist"))


def calib_vs_true(samples, outdir):
    for s, d in samples.items():
        for key, suffix, ylab in (("calib", "", "E_calib [GeV]"), ("calibP", "_prompt", "E_calib, prompt [GeV]")):
            if key not in d:
                continue
            fig, ax, cax = hs.mapfig()
            b = np.linspace(0, 52, 53)
            h, xe, ye, im = ax.hist2d(d["etrue"], d[key], bins=[b, b], cmap=hs.density_cmap, norm=LogNorm(), cmin=1)
            ax.plot([0, 52], [0, 52], color=hs.INK, lw=0.8, ls="--")
            ax.set_xlabel("E_true [GeV]")
            ax.set_ylabel(ylab)
            hs.labels(ax, EXPERIMENT, [f"{LABEL[s]}, 10° < θ < 170°"])
            fig.colorbar(im, cax=cax, label="Events")
            hs.save(fig, os.path.join(outdir, f"calib_vs_true{suffix}_{s}"))


def budget_components(d):
    """Per-event fractions of E0 for each budget component (events with a CSV entry only)."""
    b = d["budget"]
    comps = defaultdict(list)
    for i in range(len(d["etrue"])):
        r = b.get((int(d["job"][i]), int(d["event"][i])))
        if r is None:
            comps["_missing"].append(i)
            continue
        e0 = r["E0"]
        used = {"E0", "escaped"}
        calo_env = 0.0
        for c in CALOS:
            k = f"{c}_envelope"
            calo_env += r.get(k, 0.0)
            used.add(k)
        cells = d["calib"][i]
        comps["Calorimeter cells (calib. hits)"].append(cells / e0)
        if d["has_prompt"]:
            comps["Calorimeter cells, prompt"].append(d["calibP"][i] / e0)
        comps["Calorimeter structure"].append((calo_env - cells) / e0)
        for lab, rx in REGION_GROUPS:
            v = sum(e for k, e in r.items() if re.match(rx, k))
            used |= {k for k in r if re.match(rx, k)}
            comps[lab].append(v / e0)
        comps["Other (yoke, world, ...)"].append(sum(e for k, e in r.items() if k not in used) / e0)
        comps["Escaped"].append(r["escaped"] / e0)
        comps["_sum"].append(sum(e for k, e in r.items() if k not in ("E0",)) / e0)
        comps["_idx"].append(i)
    return {k: np.array(v) for k, v in comps.items()}


BUDGET_STYLE = [  # label, colour, line style
    ("Calorimeter cells (calib. hits)", hs.CAT[0], "-"),
    ("Calorimeter structure", hs.CAT[1], "--"),
    ("Solenoid", hs.CAT[2], "-."),
    ("Trackers", hs.CAT[3], ":"),
    ("Beam pipe, nozzles", hs.CAT[4], (0, (3, 1, 1, 1, 1, 1))),
    ("Other (yoke, world, ...)", hs.CAT[5], (0, (5, 2))),
    ("Escaped", hs.CAT[6], (0, (1, 1))),
    ("Calorimeter cells, prompt", hs.CAT[7], (0, (8, 2, 1, 2))),
]


SHORT = {"Calorimeter cells (calib. hits)": "Calo. cells", "Calorimeter structure": "Calo. structure",
         "Calorimeter cells, prompt": "Calo. cells, prompt",
         "Beam pipe, nozzles": "Beam pipe + nozzles", "Other (yoke, world, ...)": "Other"}


def budget_vs_theta(samples, outdir):
    for s, d in samples.items():
        if d["budget"] is None:
            continue
        comps = budget_components(d)
        idx = comps["_idx"].astype(int)
        theta = d["theta"][idx]
        edges = np.linspace(10, 170, 33)
        c = 0.5 * (edges[1:] + edges[:-1])
        bi = np.digitize(theta, edges) - 1
        fig, ax = hs.squarefig()
        for lab, col, ls in BUDGET_STYLE:
            if lab not in comps:
                continue
            mean = np.array([comps[lab][bi == k].mean() if (bi == k).sum() else np.nan for k in range(len(c))])
            mean = np.where(mean > 0, mean, np.nan)
            short = SHORT.get(lab, lab)
            if not np.any(mean >= 1e-5):
                short += " (< 10⁻⁵)"
            ax.stairs(mean, edges, color=col, ls=ls, lw=1.6, label=short)
        ax.set_yscale("log")
        ax.set_ylim(1e-5, 3e4)
        ax.set_xlim(edges[0], edges[-1])
        ax.set_xlabel("θ [deg]")
        ax.set_ylabel("Mean fraction of E_true")
        hs.labels(ax, EXPERIMENT, [f"{LABEL[s]}, 1 < p < 50 GeV", f"{len(idx)} events"])
        hs.legend(ax)
        hs.save(fig, os.path.join(outdir, f"budget_vs_theta_{s}"))


def weighted_late(x, late, allt, edges):
    """Energy-weighted late fraction sum(late)/sum(all) per bin, with its statistical uncertainty."""
    c = 0.5 * (edges[1:] + edges[:-1])
    f, e = np.full(len(c), np.nan), np.full(len(c), np.nan)
    idx = np.digitize(x, edges) - 1
    for k in range(len(c)):
        m = (idx == k) & (allt > 0)
        if m.sum() < 5:
            continue
        w, r = allt[m], late[m] / allt[m]
        f[k] = late[m].sum() / w.sum()
        e[k] = np.sqrt(np.sum(w * w * (r - f[k]) ** 2)) / w.sum()
    return c, f, e


def late_profiles(samples, xkey, edges, xlabel, name, outdir, lines):
    fig, ax = hs.squarefig()
    for s, d in samples.items():
        if not d["has_prompt"]:
            continue
        for grp, marker in (("ECal", "o"), ("HCal", "s")):
            allt, pr = d.get(f"calib_{grp}"), d.get(f"calibP_{grp}")
            if allt is None or pr is None:
                continue
            c, f, e = weighted_late(d[xkey], allt - pr, allt, edges)
            f = np.where(f > 0, f, np.nan)
            hs.mc_points(ax, c, f, e, color=SAMPLE_COLOR[s], marker=marker, label=f"{LABEL[s]}, {grp}")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 30)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Late fraction 1 − E_prompt / E_all")
    hs.labels(ax, EXPERIMENT, lines)
    hs.legend(ax)
    hs.save(fig, os.path.join(outdir, name))


def late_components(samples, outdir, edges):
    for s, d in samples.items():
        if not d["has_prompt"]:
            continue
        fig, ax = hs.squarefig()
        comps = [("total", d["calib"], d["calibP"], hs.INK, "o"),
                 ("EM visible", d["em"], d["emP"], hs.CAT[3], "s"),
                 ("non-EM visible", d["calib"] - d["inv"] - d["em"], d["calibP"] - d["invP"] - d["emP"], hs.CAT[2], "D"),
                 ("invisible", d["inv"], d["invP"], hs.CAT[5], "^")]
        for lab, allt, pr, col, marker in comps:
            c, f, e = weighted_late(d["etrue"], allt - pr, allt, edges)
            f = np.where(f > 0, f, np.nan)
            hs.mc_points(ax, c, f, e, color=col if col != hs.INK else hs.MUTED, marker=marker, label=lab)
        ax.set_yscale("log")
        ax.set_ylim(1e-4, 30)
        ax.set_xlim(edges[0], edges[-1])
        ax.set_xlabel("E_true [GeV]")
        ax.set_ylabel("Late fraction 1 − E_prompt / E_all")
        hs.labels(ax, EXPERIMENT, [f"{LABEL[s]}, ECal + HCal", "10° < θ < 170°"])
        hs.legend(ax)
        hs.save(fig, os.path.join(outdir, f"late_components_{s}"))


def bookkeeping_closure(samples, outdir):
    if all(d["budget"] is None for d in samples.values()):
        return
    fig, ax = hs.squarefig()
    bins = np.linspace(-16, -2, 57)
    top = 0
    for s, d in samples.items():
        if d["budget"] is None:
            continue
        dev = np.abs(budget_components(d)["_sum"] - 1.0)
        h, _ = np.histogram(np.clip(np.log10(np.maximum(dev, 1e-16)), bins[0], bins[-1] - 1e-9), bins=bins)
        ax.stairs(h, bins, color=SAMPLE_COLOR[s], lw=1.6, label=f"{LABEL[s]} (max {dev.max():.0e})")
        top = max(top, h.max())
    ax.axvspan(-16, -8, color=hs.NOT_USED, alpha=0.35, lw=0, zorder=0)
    ax.annotate("CSV rounding", xy=(-12, 0.02), xycoords=("data", "axes fraction"),
                ha="center", va="bottom", fontsize=8, color=hs.MUTED)
    ax.set_ylim(0, top * 1.6)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel("log₁₀ |Σ components / E_true − 1|")
    ax.set_ylabel("Events")
    hs.labels(ax, EXPERIMENT, ["energy bookkeeping"])
    hs.legend(ax)
    hs.save(fig, os.path.join(outdir, "bookkeeping_closure"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", action="append", required=True, help="name=DIR, name in e, pi, ...")
    ap.add_argument("-o", "--outdir", default="closure")
    args = ap.parse_args()
    hs.apply()
    samples = {}
    for spec in args.sample:
        name, d = spec.split("=", 1)
        samples[name] = load_dir(d)
        n = len(samples[name]["etrue"])
        r = samples[name]["calib"] / samples[name]["etrue"]
        msg = (f"{name}: {n} events, median E_calib/E_true = {np.median(r):.4f}, "
               f"budget CSV: {'yes' if samples[name]['budget'] else 'no'}")
        if samples[name]["has_prompt"]:
            d = samples[name]
            msg += (f", prompt median {np.median(d['calibP'] / d['etrue']):.4f}, "
                    f"late fraction {1 - d['calibP'].sum() / d['calib'].sum():.4f}")
        print(msg)
    os.makedirs(args.outdir, exist_ok=True)
    closure_profile(samples, "etrue", np.array([1, 2, 3, 5, 7.5, 10, 15, 20, 30, 40, 50.]),
                    "E_true [GeV]", "closure_vs_E", args.outdir, ["10° < θ < 170°"])
    closure_profile(samples, "theta", np.linspace(10, 170, 33), "θ [deg]", "closure_vs_theta",
                    args.outdir, ["1 < p < 50 GeV"])
    closure_dist(samples, args.outdir)
    calib_vs_true(samples, args.outdir)
    budget_vs_theta(samples, args.outdir)
    bookkeeping_closure(samples, args.outdir)
    e_edges = np.array([1, 2, 3, 5, 7.5, 10, 15, 20, 30, 40, 50.])
    late_profiles(samples, "etrue", e_edges, "E_true [GeV]", "late_vs_E", args.outdir,
                  ["10° < θ < 170°", "prompt: digitiser time windows"])
    late_profiles(samples, "theta", np.linspace(10, 170, 17), "θ [deg]", "late_vs_theta", args.outdir,
                  ["1 < p < 50 GeV", "prompt: digitiser time windows"])
    late_components(samples, args.outdir, e_edges)
    print("written to", args.outdir)


if __name__ == "__main__":
    main()
