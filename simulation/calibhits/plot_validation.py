#!/usr/bin/env python3
"""Validation plots for calorimeter calibration hits (steer_calib.py output).

  python3 plot_validation.py --electron e.root --pion pi.root \
      [--calib calib.root --baseline baseline.root] --outdir plots/

--electron / --pion: endcap-pointing particle-gun samples produced with
    MergeCalibContributions=False, so every Geant4 step keeps its own contribution
    and midpoint position. This allows an exact per-slice breakdown.
--calib / --baseline: optional same-seed pair (steer_calib.py vs steer_baseline.py)
    for the "active hits unchanged" plot.

Run inside the mucoll sim container (needs uproot, numpy, matplotlib).
"""
import argparse
import os

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

# --------------------------------------------------------------------------
# Geometry constants (k4geo MuColl/MAIA/compact/MAIA_v0)
# Layer stacks: (label, category, thickness [mm]); the active slice is marked.
ECAL_STACK = [("W", "Absorber", 2.20), ("G10", "G10", 0.15), ("HV/ground", "HV/ground", 0.10),
              ("Si", "Active", 0.50), ("Air", "Air", 0.10), ("PCB", "PCB", 1.30),
              ("Air", "Air", 0.25), ("G10", "G10", 0.75)]
HCAL_STACK = [("Steel", "Absorber", 0.5), ("Steel", "Absorber", 19.0), ("Scint.", "Active", 3.0),
              ("Cu", "Copper", 0.1), ("PCB", "PCB", 0.7), ("Steel", "Absorber", 0.5),
              ("Air", "Air", 2.7)]
CELL = {"ECal": 5.1, "HCal": 30.0}                       # cell size [mm]
EM_FACTOR = {"ECal": 0.00641222630095 / 0.0001575,       # MAIAConfig CaloDigi:
             "HCal": 0.0287783798145 / 0.0004925}        # calibration_factorsMipGev / calibration_mip
LAYER_SHIFT, LAYER_MASK = 19, 0x1FF                      # "layer:9" after system:5,side:-2,module:8,stave:4

# Categorical colours: reference palette slots 1-7, fixed order (validated for adjacent pairs)
CATEGORIES = ["Absorber", "Active", "PCB", "G10", "HV/ground", "Copper", "Air"]
CAT_COLOR = dict(zip(CATEGORIES, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                                  "#e87ba4", "#008300", "#6250d6"]))
PART_COLOR = {"e⁻": "#2a78d6", "π⁻": "#eb6834"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "figure.dpi": 150, "savefig.bbox": "tight",
    "figure.facecolor": "white", "axes.facecolor": "#fcfcfb",
})


def stack_bounds(stack):
    """Slice boundaries relative to the centre of the active slice [mm]."""
    t = np.array([s[2] for s in stack])
    edges = np.concatenate([[0.0], np.cumsum(t)])
    ia = [s[1] for s in stack].index("Active")
    centre = 0.5 * (edges[ia] + edges[ia + 1])
    return edges - centre


# --------------------------------------------------------------------------
def _br(tree, coll, field):
    return tree[f"{coll}/{coll}.{field}"].array(library="np")


def load_sample(path, calos=("ECalEndcap", "HCalEndcap", "ECalBarrel", "HCalBarrel")):
    """Per calorimeter: active hits, calib hits, and (if present) calib contributions."""
    tree = uproot.open(path)["events"]
    data = {}
    for calo in calos:
        act, cal = f"{calo}Collection", f"{calo}CollectionCalib"
        d = {"act": {f: _br(tree, act, f) for f in ("cellID", "energy")},
             "cal": {f: _br(tree, cal, f) for f in ("cellID", "energy", "position.x", "position.y",
                                                    "position.z", "contributions_begin", "contributions_end")}}
        ccol = f"{cal}Contributions"
        if f"{ccol}/{ccol}.stepPosition.x" in tree:
            d["con"] = {f: _br(tree, ccol, f) for f in ("energy", "stepPosition.x", "stepPosition.y", "stepPosition.z")}
            d["rel"] = tree[f"_{cal}_contributions/_{cal}_contributions.index"].array(library="np")
        data[calo] = d
    return data


def contributions_vs_cell(d):
    """Flatten contributions with their offset from the centre of the cell they were assigned to."""
    out = {k: [] for k in ("e", "dx", "dy", "depth", "closure")}
    for iev in range(len(d["cal"]["energy"])):
        b, e = d["cal"]["contributions_begin"][iev], d["cal"]["contributions_end"][iev]
        if len(b) == 0:
            continue
        hit = np.repeat(np.arange(len(b)), e - b)
        rel = np.concatenate([np.arange(bb, ee) for bb, ee in zip(b, e)])
        ci = d["rel"][iev][rel]
        ce = d["con"]["energy"][iev][ci]
        hx, hy, hz = (d["cal"][f"position.{a}"][iev][hit] for a in "xyz")
        sx, sy, sz = (d["con"][f"stepPosition.{a}"][iev][ci] for a in "xyz")
        out["e"].append(ce)
        out["dx"].append(sx - hx)
        out["dy"].append(sy - hy)
        out["depth"].append((sz - hz) * np.sign(hz))   # endcaps: depth axis = +-z
        # closure: sum of contributions == hit energy
        out["closure"].append(np.bincount(hit, weights=ce, minlength=len(b)) - d["cal"]["energy"][iev])
    return {k: np.concatenate(v) if v else np.array([]) for k, v in out.items()}


def cell_table(d):
    """Per (event, cell): active energy (0 if none) and calib energy."""
    ea, ec, lay = [], [], []
    for iev in range(len(d["cal"]["energy"])):
        amap = dict(zip(d["act"]["cellID"][iev].tolist(), d["act"]["energy"][iev].tolist()))
        cid = d["cal"]["cellID"][iev]
        ec.append(d["cal"]["energy"][iev])
        ea.append(np.array([amap.get(c, 0.0) for c in cid.tolist()]))
        lay.append((cid.astype(np.uint64) >> np.uint64(LAYER_SHIFT)) & np.uint64(LAYER_MASK))
    cat = lambda v: np.concatenate(v) if v else np.array([])  # noqa: E731
    return {"act": cat(ea), "cal": cat(ec), "layer": cat(lay).astype(int)}


# --------------------------------------------------------------------------
def fig_material(samples, outdir):
    """100% stacked bars: what the calibration-hit energy is made of."""
    rows = []
    for pname, calo, stack in (("e⁻", "ECal", ECAL_STACK), ("π⁻", "ECal", ECAL_STACK), ("π⁻", "HCal", HCAL_STACK)):
        c = samples[pname]["contrib"][calo]
        edges = stack_bounds(stack)
        idx = np.clip(np.searchsorted(edges, c["depth"]) - 1, 0, len(stack) - 1)
        per_cat = {cat: 0.0 for cat in CATEGORIES}
        for i, s in enumerate(stack):
            per_cat[s[1]] += c["e"][idx == i].sum()
        tot = sum(per_cat.values())
        rows.append((f"{calo} endcap · {pname}", {k: v / tot for k, v in per_cat.items()}, tot))

    fig, ax = plt.subplots(figsize=(8.2, 3.1))
    ylab = [r[0] for r in rows][::-1]
    for iy, (_, frac, _) in enumerate(rows[::-1]):
        left = 0.0
        for cat in CATEGORIES:
            w = frac[cat]
            if w <= 0:
                continue
            ax.barh(iy, w - 0.002, left=left, height=0.62, color=CAT_COLOR[cat], edgecolor="white", linewidth=1.0)
            if w > 0.08:
                ax.text(left + w / 2, iy, f"{100 * w:.1f}%", ha="center", va="center",
                        fontsize=8.5, color="white" if cat in ("Absorber", "PCB") else INK)
            left += w
        ax.text(1.015, iy, f"active {100 * frac['Active']:.2f}%\n= 1 / {1 / frac['Active']:.1f}",
                va="center", ha="left", fontsize=8.5, color=INK, clip_on=False)
    ax.set_yticks(range(len(rows)), ylab)
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax.set_xlabel("Share of calibration-hit energy (20 GeV, θ = 15–35°)")
    ax.grid(axis="y", visible=False)
    handles = [Rectangle((0, 0), 1, 1, color=CAT_COLOR[c]) for c in CATEGORIES]
    ax.legend(handles, ["Absorber (W / steel)", "Active (Si / scint.)", "PCB", "G10", "HV/ground", "Copper", "Air"],
              ncol=7, loc="lower left", bbox_to_anchor=(0, 1.02), fontsize=8, handlelength=1.0, columnspacing=1.0)
    ax.set_title("What a calibration hit contains: energy by slice material", loc="left", pad=28)
    fig.savefig(os.path.join(outdir, "1_material_breakdown.png"))
    plt.close(fig)
    return rows


def fig_depth(samples, outdir):
    """Energy-weighted depth of deposits relative to the active slice of their cell."""
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.2))
    for ax, (pname, calo, stack, bw) in zip(axes, (("e⁻", "ECal", ECAL_STACK, 0.025), ("π⁻", "HCal", HCAL_STACK, 0.1))):
        c = samples[pname]["contrib"][calo]
        edges = stack_bounds(stack)
        bins = np.arange(edges[0] - 3 * bw, edges[-1] + 3 * bw + 1e-9, bw)
        h, _ = np.histogram(c["depth"], bins=bins, weights=c["e"])
        h = h / (c["e"].sum() * bw)
        centres = 0.5 * (bins[1:] + bins[:-1])
        idx = np.searchsorted(edges, centres) - 1
        for i, s in enumerate(stack):
            m = idx == i
            ax.bar(centres[m], h[m], width=bw, color=CAT_COLOR[s[1]], linewidth=0)
            ax.axvline(edges[i], color=INK2, lw=0.4, alpha=0.5)
        ax.axvline(edges[-1], color=INK2, lw=0.4, alpha=0.5)
        outside = c["e"][(c["depth"] < edges[0] - 1e-3) | (c["depth"] > edges[-1] + 1e-3)].sum() / c["e"].sum()
        ax.set_yscale("log")
        ax.set_xlim(bins[0], bins[-1])
        ax.set_xlabel("Depth of energy deposit relative to centre of active slice [mm]")
        ax.set_ylabel("Energy fraction per mm")
        ax.set_title(f"{calo} endcap · {pname}: all layers folded onto one", loc="left")
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi * 40)
        trans = matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
        for i, s in enumerate(stack):
            if edges[i + 1] - edges[i] >= (0.2 if calo == "ECal" else 1.5):
                ax.text(0.5 * (edges[i] + edges[i + 1]), 0.97, s[0], transform=trans,
                        ha="center", va="top", fontsize=8, color=INK2)
        ax.text(0.99, 0.84, f"energy outside layer: {100 * outside:.3f}%", transform=ax.transAxes,
                ha="right", fontsize=8, color=INK2)
    fig.suptitle("Every deposit lies inside the layer of the cell it was assigned to "
                 "(spikes = midpoints of tracks crossing a whole slice)", x=0.01, ha="left", y=1.02, fontsize=11)
    fig.savefig(os.path.join(outdir, "2_depth_profile.png"))
    plt.close(fig)


def fit_grid_rotation(c, half):
    """Rotation [deg] of the segmentation grid w.r.t. global x-y: minimises energy outside the cell."""
    sel = np.random.default_rng(1).choice(len(c["e"]), size=min(len(c["e"]), 200000), replace=False)
    dx, dy, e = c["dx"][sel], c["dy"][sel], c["e"][sel]

    def outside(a):
        ca, sa = np.cos(np.radians(a)), np.sin(np.radians(a))
        return e[(np.abs(ca * dx + sa * dy) > half) | (np.abs(-sa * dx + ca * dy) > half)].sum()
    coarse = np.arange(-45, 45, 0.5)
    a0 = coarse[np.argmin([outside(a) for a in coarse])]
    fine = np.arange(a0 - 0.5, a0 + 0.5, 0.01)
    return fine[np.argmin([outside(a) for a in fine])]


def fig_lateral(samples, outdir):
    """Transverse offset of deposits from the cell centre: must fill exactly one cell."""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"wspace": 0.45})
    for ax, (pname, calo) in zip(axes, (("e⁻", "ECal"), ("π⁻", "HCal"))):
        c = samples[pname]["contrib"][calo]
        half = CELL[calo] / 2
        rng = 1.6 * half
        ang = fit_grid_rotation(c, half)
        ca, sa = np.cos(np.radians(ang)), np.sin(np.radians(ang))
        u, v = ca * c["dx"] + sa * c["dy"], -sa * c["dx"] + ca * c["dy"]
        h, xe, ye, im = ax.hist2d(u, v, bins=80, range=[[-rng, rng], [-rng, rng]],
                                  weights=c["e"], cmap="Blues", norm=LogNorm())
        ax.add_patch(Rectangle((-half, -half), 2 * half, 2 * half, fill=False, ec=INK, ls="--", lw=1.0))
        out = c["e"][(np.abs(u) > half + 1e-3) | (np.abs(v) > half + 1e-3)].sum() / c["e"].sum()
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_xlabel("Δu  (deposit − cell centre, cell frame) [mm]")
        ax.set_ylabel("Δv [mm]")
        ax.set_title(f"{calo} endcap · {pname}  ({CELL[calo]:g} mm cells; grid rotated {ang:.2f}° in x-y)",
                     loc="left", fontsize=10)
        ax.text(0.02, 0.02, f"energy outside cell: {100 * out:.3f}%", transform=ax.transAxes, fontsize=8, color=INK2)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("Energy [GeV]")
    fig.suptitle("Deposits in absorber and passive slices fill exactly one cell (dashed = cell boundary)",
                 x=0.01, ha="left", y=1.0, fontsize=11)
    fig.savefig(os.path.join(outdir, "3_lateral_in_cell.png"))
    plt.close(fig)


def fig_cells(samples, outdir):
    """Per cell: calibration energy vs active energy, plus what is in cells without an active hit."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    notes = {}
    for ax, calo in zip(axes, ("ECal", "HCal")):
        lo, hi = 1e-7, 1.0
        for pname in ("e⁻", "π⁻"):
            t = samples[pname]["cells"][calo]
            m = t["act"] > 0
            ax.scatter(t["act"][m], t["cal"][m], s=3, alpha=0.25, color=PART_COLOR[pname], linewidths=0,
                       label=pname, rasterized=True)
            frac_dark = t["cal"][~m].sum() / t["cal"].sum() if t["cal"].sum() > 0 else np.nan
            notes.setdefault(calo, []).append(
                f"{pname}: {100 * frac_dark:.0f}% of calib energy in {(~m).sum()} cells with no active hit"
                f" (vs {m.sum()} with)")
        x = np.logspace(np.log10(lo), np.log10(hi), 10)
        ax.plot(x, EM_FACTOR[calo] * x, color=INK, lw=1.2, ls="--")
        ax.text(1.2e-3, 2.2, f"y = {EM_FACTOR[calo]:.1f} x\n(digitiser\nEM-scale factor)",
                fontsize=8, color=INK, ha="right", va="top")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, 1e-1)
        ax.set_ylim(1e-6, 5)
        ax.set_xlabel("Active-layer energy in cell [GeV]  (*Collection)")
        ax.set_ylabel("Calibration energy in cell [GeV]  (*CollectionCalib)")
        ax.set_title(f"{calo} endcap, cells with an active hit", loc="left")
        leg = ax.legend(loc="upper left", markerscale=4)
        for lh in leg.legend_handles:
            lh.set_alpha(1)
    fig.suptitle("Per cell, calibration energy tracks the active energy × the sampling factor", x=0.01, ha="left",
                 y=1.0, fontsize=11)
    fig.savefig(os.path.join(outdir, "4_cells_calib_vs_active.png"))
    plt.close(fig)
    return notes


def fig_dark_cells(samples, outdir):
    """Share of calibration energy in cells that do / do not have an active hit."""
    rows = [("ECal endcap · e⁻", "e⁻", "ECal"), ("ECal endcap · π⁻", "π⁻", "ECal"), ("HCal endcap · π⁻", "π⁻", "HCal")]
    fig, ax = plt.subplots(figsize=(8.2, 2.7))
    cols = {"with": "#2a78d6", "without": "#eb6834"}
    for iy, (lab, pname, calo) in enumerate(rows[::-1]):
        t = samples[pname]["cells"][calo]
        m = t["act"] > 0
        tot = t["cal"].sum()
        fw = t["cal"][m].sum() / tot
        ax.barh(iy, fw - 0.002, height=0.6, color=cols["with"], edgecolor="white", lw=1)
        ax.barh(iy, 1 - fw, left=fw, height=0.6, color=cols["without"], edgecolor="white", lw=1)
        ax.text(fw / 2, iy, f"{100 * fw:.0f}%  ({m.sum()} cells)", ha="center", va="center", fontsize=8.5, color="white")
        ax.text(fw + (1 - fw) / 2, iy, f"{100 * (1 - fw):.0f}%  ({(~m).sum()} cells)", ha="center", va="center",
                fontsize=8.5, color=INK)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows][::-1])
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Share of calibration-hit energy (all events)")
    ax.legend([Rectangle((0, 0), 1, 1, color=cols[k]) for k in ("with", "without")],
              ["cells with an active hit (can be clustered)", "cells with no active hit (never digitised)"],
              ncol=2, loc="lower left", bbox_to_anchor=(0, 1.02), fontsize=8.5)
    ax.set_title("Much of the deposited energy sits in cells the reconstruction never sees", loc="left", pad=24)
    fig.savefig(os.path.join(outdir, "7_cells_without_active_hit.png"))
    plt.close(fig)


def fig_longitudinal(samples, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.9), gridspec_kw={"width_ratios": [50, 75]})
    for ax, calo, nl in zip(axes, ("ECal", "HCal"), (50, 75)):
        for pname in ("e⁻", "π⁻"):
            t = samples[pname]["cells"][calo]
            nev = samples[pname]["nev"]
            cal = np.bincount(t["layer"], weights=t["cal"], minlength=nl)[:nl] / nev
            act = np.bincount(t["layer"], weights=t["act"], minlength=nl)[:nl] / nev * EM_FACTOR[calo]
            ax.plot(np.arange(nl), cal, color=PART_COLOR[pname], lw=2, label=f"{pname} calibration hits")
            ax.plot(np.arange(nl), act, color=PART_COLOR[pname], lw=1.2, ls=(0, (3, 2)),
                    label=f"{pname} active × EM-scale factor ({EM_FACTOR[calo]:.1f})")
        ax.set_yscale("log")
        ax.set_xlabel("Layer")
        ax.set_ylabel("Mean energy per event [GeV]")
        ax.set_title(f"{calo} endcap", loc="left")
    for ax, loc in zip(axes, ("lower left", "upper right")):
        ax.legend(fontsize=7.5, loc=loc)
    fig.suptitle("Longitudinal profiles: calibration hits vs EM-scale active energy", x=0.01, ha="left", y=1.0, fontsize=11)
    fig.savefig(os.path.join(outdir, "5_longitudinal.png"))
    plt.close(fig)


def fig_totals(samples, outdir):
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    for pname in ("e⁻", "π⁻"):
        s = samples[pname]
        ax.scatter(s["tot_em"], s["tot_cal"], s=36, color=PART_COLOR[pname], edgecolor="white", lw=1.0, label=pname, zorder=3)
    ax.plot([0, 25], [0, 25], color=INK2, lw=1, ls="--")
    ax.text(23.2, 22.0, "y = x", fontsize=8, color=INK2, rotation=45)
    ax.axhline(20, color=INK2, lw=0.8, ls=":")
    ax.text(0.5, 20.3, "beam energy 20 GeV", fontsize=8, color=INK2)
    ax.set_xlim(0, 25)
    ax.set_ylim(0, 25)
    ax.set_aspect("equal")
    ax.set_xlabel("Active energy × EM-scale factors, ECal + HCal [GeV]")
    ax.set_ylabel("Calibration-hit energy, ECal + HCal [GeV]")
    ax.set_title("Per event totals (10 events each)", loc="left")
    ax.legend(loc="lower right")
    fig.savefig(os.path.join(outdir, "6_event_totals.png"))
    plt.close(fig)


def fig_baseline(calib_path, base_path, outdir):
    calos = ("ECalBarrel", "ECalEndcap", "HCalBarrel", "HCalEndcap")
    tc, tb = uproot.open(calib_path)["events"], uproot.open(base_path)["events"]
    xs, ys, nmiss, nhits = [], [], 0, 0
    for calo in calos:
        c = f"{calo}Collection"
        cid_c, e_c = _br(tc, c, "cellID"), _br(tc, c, "energy")
        cid_b, e_b = _br(tb, c, "cellID"), _br(tb, c, "energy")
        for iev in range(len(e_b)):
            mc = dict(zip(cid_c[iev].tolist(), e_c[iev].tolist()))
            nmiss += len(set(mc) ^ set(cid_b[iev].tolist()))
            for cid, eb in zip(cid_b[iev].tolist(), e_b[iev].tolist()):
                if cid in mc:
                    xs.append(eb)
                    ys.append(mc[cid])
                    nhits += 1
    xs, ys = np.array(xs), np.array(ys)
    maxrel = np.max(np.abs(ys - xs) / np.maximum(xs, 1e-30)) if len(xs) else np.nan
    fig, ax = plt.subplots(figsize=(4.8, 4.4))
    ax.scatter(xs, ys, s=4, color=PART_COLOR["e⁻"], alpha=0.5, linewidths=0, rasterized=True)
    ax.plot([1e-9, 1], [1e-9, 1], color=INK2, lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    lo = max(xs.min() * 0.5, 1e-9)
    ax.set_xlim(lo, xs.max() * 2)
    ax.set_ylim(lo, xs.max() * 2)
    ax.set_xlabel("Hit energy, steer_baseline.py [GeV]")
    ax.set_ylabel("Hit energy, steer_calib.py [GeV]")
    ax.set_title("Active hits unchanged (same seed)", loc="left")
    ax.text(0.04, 0.96, f"{nhits} hits, all 4 calorimeters\ncells only in one run: {nmiss}\n"
            f"max |ΔE|/E: {maxrel:.1e}", transform=ax.transAxes, va="top", fontsize=8, color=INK2)
    fig.savefig(os.path.join(outdir, "0_baseline_identity.png"))
    plt.close(fig)
    return nhits, nmiss, maxrel


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--electron", required=True)
    ap.add_argument("--pion", required=True)
    ap.add_argument("--calib")
    ap.add_argument("--baseline")
    ap.add_argument("--outdir", default="calibhits_validation")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    samples = {}
    for pname, path in (("e⁻", args.electron), ("π⁻", args.pion)):
        raw = load_sample(path)
        s = {"contrib": {}, "cells": {}}
        s["contrib"]["ECal"] = contributions_vs_cell(raw["ECalEndcap"])
        s["contrib"]["HCal"] = contributions_vs_cell(raw["HCalEndcap"])
        s["cells"]["ECal"] = cell_table(raw["ECalEndcap"])
        s["cells"]["HCal"] = cell_table(raw["HCalEndcap"])
        s["nev"] = len(raw["ECalEndcap"]["cal"]["energy"])
        tot_cal = np.zeros(s["nev"])
        tot_em = np.zeros(s["nev"])
        for calo, d in raw.items():
            f = EM_FACTOR["ECal" if calo.startswith("ECal") else "HCal"]
            tot_cal += np.array([x.sum() for x in d["cal"]["energy"]])
            tot_em += f * np.array([x.sum() for x in d["act"]["energy"]])
        s["tot_cal"], s["tot_em"] = tot_cal, tot_em
        for calo in ("ECal", "HCal"):
            cl = s["contrib"][calo]["closure"]
            print(f"{pname} {calo}: {len(s['contrib'][calo]['e'])} contributions; "
                  f"max |sum(contrib) - hit E| = {np.abs(cl).max() if len(cl) else 0:.2e} GeV")
        samples[pname] = s

    rows = fig_material(samples, args.outdir)
    for name, frac, tot in rows:
        print(f"{name}: " + ", ".join(f"{k} {100 * v:.2f}%" for k, v in frac.items() if v > 0))
    fig_depth(samples, args.outdir)
    fig_lateral(samples, args.outdir)
    for calo, lines in fig_cells(samples, args.outdir).items():
        print(calo, "|", " ; ".join(lines))
    fig_dark_cells(samples, args.outdir)
    fig_longitudinal(samples, args.outdir)
    fig_totals(samples, args.outdir)
    for pname in samples:
        print(pname, "event totals: calib", np.round(samples[pname]["tot_cal"], 2),
              "active*EM", np.round(samples[pname]["tot_em"], 2))
    if args.calib and args.baseline:
        print("baseline identity: hits=%d, cells only in one run=%d, max rel diff=%.1e"
              % fig_baseline(args.calib, args.baseline, args.outdir))
    print("plots written to", args.outdir)


if __name__ == "__main__":
    main()
