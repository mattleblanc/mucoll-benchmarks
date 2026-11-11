#!/usr/bin/env python3
"""Event displays for calibration-hit validation (2D projections, house style).

    python3 event_display.py reco.edm4hep.root -e 0 1 2 -o displays/

For every event, writes one pdf + png per plot:
  evt<N>_overview_rz        whole detector, r-z: calibration energy + true trajectory
  evt<N>_<view>_calib       calibration-hit energy (all slices, visible + invisible)
  evt<N>_<view>_active      active-layer (Si / scintillator) energy, SimCalorimeterHits
  evt<N>_<view>_invisible   invisible energy (positive part; net value in the label)
  evt<N>_<view>_emfrac      EM fraction of the visible deposited energy (CalibEM / visible)
  evt<N>_<view>_ratio       calibration / active energy per pixel (pixels with active energy):
                            the local factor a cell-level calibration has to learn
  evt<N>_<view>_clusters    cells by Pandora cluster / active hit but unclustered / no active hit
With the prompt collections (*CalibPrompt) also:
  evt<N>_<view>_calib_prompt  calibration energy inside the digitiser time windows
  evt<N>_<view>_late          calibration energy outside them (all - prompt)
  evt<N>_<view>_clusters_prompt  the cluster view with prompt calibration energy
with <view> = rz and xy, all zoomed on the shower with identical binning.

The true trajectory of the primary joins the generator vertex, the primary's own
tracker SimHits (hits made by secondaries are skipped) and its end point.
Thin lines are the other saved MCParticles (vertex -> end point). Their
parent is the nearest *saved* ancestor, so albedo (e.g. neutrons
re-entering the tracker) can show up far from the shower.

Run inside the mucoll sim container (podio, numpy, matplotlib).
"""
import argparse
import math
import os
import sys

import numpy as np
from podio.root_io import Reader

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import house_style as hs  # noqa: E402
import matplotlib  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

CALOS = ["ECalBarrel", "ECalEndcap", "HCalBarrel", "HCalEndcap"]
TRACKER = ["VertexBarrelCollection", "VertexEndcapCollection", "InnerTrackerBarrelCollection",
           "InnerTrackerEndcapCollection", "OuterTrackerBarrelCollection", "OuterTrackerEndcapCollection"]
# Calorimeter envelopes (MAIA_v0 compact) [mm]: (zmin, zmax, rmin, rmax); barrels are dodecagons
ENVELOPES = [(-2307, 2307, 1857, 2124.5), (2307, 2574.5, 310, 2124.5), (-2574.5, -2307, 310, 2124.5),
             (-2574.5, 2574.5, 2126, 4113.5), (2575.5, 4600, 310, 4113.5), (-4600, -2575.5, 310, 4113.5)]
SOLENOID_R = (1500, 1857)
PDG_NAME = {11: "e⁻", -11: "e⁺", 22: "γ", 211: "π⁺", -211: "π⁻", 2112: "n", 130: "K⁰_L", 13: "μ⁻", -13: "μ⁺"}

TRAJ_COLOR = hs.CAT[1]      # rose
SEC_COLOR = hs.CAT[6]       # teal
GREY_ACTIVE, GREY_DARK = "#8C8B87", "#CFCECA"


# --------------------------------------------------------------------------
def cells(coll):
    """cellID -> (x, y, z, E) for a SimCalorimeterHit collection."""
    out = {}
    for h in coll:
        p = h.getPosition()
        out[h.getCellID()] = (p.x, p.y, p.z, h.getEnergy())
    return out


def extract(ev):
    names = set(str(n) for n in ev.getAvailableCollections())
    d = {k: {} for k in ("active", "calib", "invisible", "em", "calibP")}
    for c in CALOS:
        for key, suffix in (("active", ""), ("calib", "Calib"), ("invisible", "CalibInvisible"), ("em", "CalibEM"),
                            ("calibP", "CalibPrompt")):
            name = f"{c}Collection{suffix}"
            if name in names:
                d[key].update(cells(ev.get(name)))
    # late energy per cell = all-time - prompt calibration energy
    d["late"] = {cid: (x, y, z, e - d["calibP"].get(cid, (0, 0, 0, 0.0))[3])
                 for cid, (x, y, z, e) in d["calib"].items()} if d["calibP"] else {}
    mcs = ev.get("MCParticles")
    prim = [p for p in mcs if p.getGeneratorStatus() == 1]
    p0 = prim[0]
    pts = []
    for cname in TRACKER:
        if cname not in names:
            continue
        for h in ev.get(cname):
            if h.isProducedBySecondary():
                continue
            if h.getParticle().getObjectID().index != p0.getObjectID().index:
                continue
            x = h.getPosition()
            pts.append((h.getTime(), x.x, x.y, x.z))
    v, e = p0.getVertex(), p0.getEndpoint()
    traj = [(v.x, v.y, v.z)] + [p[1:] for p in sorted(pts)] + [(e.x, e.y, e.z)]
    secs = []
    for p in mcs:
        if p.getGeneratorStatus() == 1:
            continue
        a, b = p.getVertex(), p.getEndpoint()
        if p.getEnergy() > 0.01 or math.dist((a.x, a.y, a.z), (b.x, b.y, b.z)) > 20:
            secs.append(((a.x, a.y, a.z), (b.x, b.y, b.z)))
    m = p0.getMomentum()
    pmag = math.sqrt(m.x ** 2 + m.y ** 2 + m.z ** 2)
    d["truth"] = {"pdg": p0.getPDG(), "E": p0.getEnergy(), "theta": math.degrees(math.acos(m.z / pmag)),
                  "phi": math.degrees(math.atan2(m.y, m.x)), "traj": np.array(traj), "secs": secs}
    # Pandora clusters -> cellIDs, ranked by energy
    clus = []
    if "PandoraClusters" in names:
        for cl in ev.get("PandoraClusters"):
            clus.append((cl.getEnergy(), set(h.getCellID() for h in cl.getHits())))
    d["clusters"] = sorted(clus, key=lambda t: -t[0])
    return d


# --------------------------------------------------------------------------
def proj(xyz, view):
    xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)
    if view == "rz":
        return xyz[:, 2], np.hypot(xyz[:, 0], xyz[:, 1])
    if view == "rzs":   # signed r: + for y >= 0, - for y < 0 (full r-z plane)
        return xyz[:, 2], np.where(xyz[:, 1] >= 0, 1.0, -1.0) * np.hypot(xyz[:, 0], xyz[:, 1])
    return xyz[:, 0], xyz[:, 1]


def arrays(cmap):
    if not cmap:
        return np.zeros((0, 3)), np.zeros(0)
    a = np.array(list(cmap.values()))
    return a[:, :3], a[:, 3]


def zoom_window(d, view):
    """Square window around the calib energy (1-99% energy-weighted quantiles) incl. the shower start."""
    xyz, e = arrays(d["calib"])
    u, v = proj(xyz, view)
    w = np.clip(e, 0, None)

    def wq(a, q):
        o = np.argsort(a)
        c = np.cumsum(w[o]) / w.sum()
        return a[o][np.searchsorted(c, q)]
    ulo, uhi, vlo, vhi = wq(u, 0.01), wq(u, 0.99), wq(v, 0.01), wq(v, 0.99)
    tu, tv = proj(d["truth"]["traj"][-1:], view)          # primary end point
    ulo, uhi = min(ulo, tu[0]), max(uhi, tu[0])
    vlo, vhi = min(vlo, tv[0]), max(vhi, tv[0])
    half = max(uhi - ulo, vhi - vlo, 250.0) * 0.5 * 1.15
    uc, vc = 0.5 * (ulo + uhi), 0.5 * (vlo + vhi)
    head = 0.75 * half          # head space for the label lines and the legend
    half_new = half + 0.5 * head
    return (uc - half_new, uc + half_new), (vc - half, vc + half + head)


def draw_envelopes(ax, view):
    kw = dict(fill=False, ec=hs.MUTED, lw=0.6, ls="-", zorder=3)
    if view in ("rz", "rzs"):
        for z0, z1, r0, r1 in ENVELOPES:
            ax.add_patch(Rectangle((z0, r0), z1 - z0, r1 - r0, **kw))
            if view == "rzs":
                ax.add_patch(Rectangle((z0, -r1), z1 - z0, r1 - r0, **kw))
    else:
        t = np.linspace(0, 2 * np.pi, 361)
        for r in (1857, 2124.5, 4113.5):
            ax.plot(r * np.cos(t), r * np.sin(t), color=hs.MUTED, lw=0.6, zorder=3)


def draw_truth(ax, d, view, lim=None):
    for a, b in d["truth"]["secs"]:
        u, v = proj([a, b], view)
        ax.plot(u, v, color=SEC_COLOR, lw=0.7, alpha=0.8, zorder=4)
    u, v = proj(d["truth"]["traj"], view)
    ax.plot(u, v, color=TRAJ_COLOR, lw=1.8, zorder=5)
    ax.plot(u[-1:], v[-1:], marker="o", ms=6, color=TRAJ_COLOR, mec="white", mew=0.8, ls="none", zorder=6)
    return [Line2D([], [], color=TRAJ_COLOR, lw=1.8, marker="o", ms=5, mec="white", label="true track"),
            Line2D([], [], color=SEC_COLOR, lw=0.8, label="saved MC sec.")]


def label_lines(d, iev, extra=()):
    t = d["truth"]
    name = PDG_NAME.get(t["pdg"], str(t["pdg"]))
    return [f"{name}, E = {t['E']:.1f} GeV", f"θ = {t['theta']:.0f}°, φ = {t['phi']:.0f}°, event {iev}"] + list(extra)


def axis_labels(ax, view):
    if view in ("rz", "rzs"):
        ax.set_xlabel("z [mm]")
        ax.set_ylabel("r [mm]" if view == "rz" else "r × sign(y) [mm]")
    else:
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")


def finish(ax, ulim, vlim, handles, d, iev, extra=()):
    ax.set_xlim(*ulim)
    ax.set_ylim(*vlim)
    ax.set_aspect("equal", adjustable="box")
    hs.labels(ax, "MAIA simulation", label_lines(d, iev, extra))
    hs.legend(ax, handles, [h.get_label() for h in handles])


def energy_map(d, iev, view, key, outdir, cblabel, ulim, vlim, nb=110, tag=None, overview=False):
    fig, ax, cax = hs.mapfig()
    xyz, e = arrays(d[key])
    extra = []
    if key in ("invisible", "late"):
        extra = [f"net {'invisible' if key == 'invisible' else 'late'}: {e.sum():.2f} GeV"]
        e = np.clip(e, 0, None)
    u, v = proj(xyz, view)
    ub, vb = np.linspace(*ulim, nb + 1), np.linspace(*vlim, nb + 1)
    h, _, _ = np.histogram2d(u, v, bins=[ub, vb], weights=e)
    h[h <= 0] = np.nan
    vmax = np.nanmax(h) if np.isfinite(h).any() else 1.0
    im = ax.pcolormesh(ub, vb, h.T, cmap=hs.density_cmap, norm=LogNorm(vmin=vmax * 1e-4, vmax=vmax),
                       rasterized=True, zorder=2)
    draw_envelopes(ax, view)
    handles = draw_truth(ax, d, view)
    if key != "active" and not overview:
        extra = [f"total: {e.sum():.2f} GeV"] + extra if key not in ("invisible", "late") else extra
    finish(ax, ulim, vlim, handles, d, iev, extra)
    axis_labels(ax, view)
    fig.colorbar(im, cax=cax, label=cblabel)
    name = f"evt{iev}_overview_rz" if overview else f"evt{iev}_{view}_{tag or key}"
    hs.save(fig, os.path.join(outdir, name))


def emfrac_map(d, iev, view, outdir, ulim, vlim, nb=110):
    fig, ax, cax = hs.mapfig()
    xyz_c, e_c = arrays(d["calib"])
    xyz_i, e_i = arrays(d["invisible"])
    xyz_e, e_e = arrays(d["em"])
    ub, vb = np.linspace(*ulim, nb + 1), np.linspace(*vlim, nb + 1)
    vis = np.histogram2d(*proj(xyz_c, view), bins=[ub, vb], weights=e_c)[0] \
        - np.histogram2d(*proj(xyz_i, view), bins=[ub, vb], weights=e_i)[0]
    em = np.histogram2d(*proj(xyz_e, view), bins=[ub, vb], weights=e_e)[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        f = np.where(vis > 1e-6, em / vis, np.nan)
    ramp = matplotlib.colors.LinearSegmentedColormap.from_list("ordered", hs.RAMP)
    im = ax.pcolormesh(ub, vb, f.T, cmap=ramp, norm=Normalize(0, 1), rasterized=True, zorder=2)
    draw_envelopes(ax, view)
    handles = draw_truth(ax, d, view)
    fem = e_e.sum() / max(e_c.sum() - e_i.sum(), 1e-9)
    finish(ax, ulim, vlim, handles, d, iev, [f"EM fraction (all cells): {fem:.2f}"])
    axis_labels(ax, view)
    fig.colorbar(im, cax=cax, label="EM fraction of visible energy")
    hs.save(fig, os.path.join(outdir, f"evt{iev}_{view}_emfrac"))


def ratio_map(d, iev, view, outdir, ulim, vlim, nb=110):
    """Calibration energy / active energy per pixel, i.e. the local inverse sampling fraction."""
    fig, ax, cax = hs.mapfig()
    xyz_c, e_c = arrays(d["calib"])
    xyz_a, e_a = arrays(d["active"])
    ub, vb = np.linspace(*ulim, nb + 1), np.linspace(*vlim, nb + 1)
    cal = np.histogram2d(*proj(xyz_c, view), bins=[ub, vb], weights=e_c)[0]
    act = np.histogram2d(*proj(xyz_a, view), bins=[ub, vb], weights=e_a)[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(act > 0, cal / act, np.nan)
    ramp = matplotlib.colors.LinearSegmentedColormap.from_list("ordered", hs.RAMP)
    im = ax.pcolormesh(ub, vb, r.T, cmap=ramp, norm=LogNorm(5, 500), rasterized=True, zorder=2)
    draw_envelopes(ax, view)
    handles = draw_truth(ax, d, view)
    tot = e_c.sum() / max(e_a.sum(), 1e-12)
    finish(ax, ulim, vlim, handles, d, iev, [f"event total: {tot:.0f}"])
    axis_labels(ax, view)
    fig.colorbar(im, cax=cax, label="Calibration / active energy")
    hs.save(fig, os.path.join(outdir, f"evt{iev}_{view}_ratio"))


def cluster_map(d, iev, view, outdir, ulim, vlim, key="calib"):
    fig, ax = hs.squarefig()
    xyz, e = arrays(d[key])
    ids = np.array(list(d[key].keys()))
    u, v = proj(xyz, view)
    emax = max(e.max(), 1e-9)
    size = 2 + 60 * np.sqrt(np.clip(e, 0, None) / emax)
    in_clus = np.full(len(ids), -1)
    for k, (_, cset) in enumerate(d["clusters"]):
        in_clus[np.isin(ids, list(cset))] = k
    has_active = np.isin(ids, list(d["active"].keys()))
    handles = []
    m = (in_clus < 0) & ~has_active
    ax.scatter(u[m], v[m], s=size[m], color=GREY_DARK, lw=0, zorder=2)
    handles.append(Line2D([], [], ls="none", marker="o", color=GREY_DARK,
                          label=f"no active hit: {e[m].sum():.1f} GeV"))
    m = (in_clus < 0) & has_active
    ax.scatter(u[m], v[m], s=size[m], color=GREY_ACTIVE, lw=0, zorder=3)
    handles.append(Line2D([], [], ls="none", marker="o", color=GREY_ACTIVE,
                          label=f"unclustered: {e[m].sum():.1f} GeV"))
    for k in range(min(len(d["clusters"]), 7)):
        m = in_clus == k
        if not m.any():
            continue
        ax.scatter(u[m], v[m], s=size[m], color=hs.CAT[[0, 3, 2, 4, 5, 7, 6][k]], lw=0, zorder=4 + k)
        handles.append(Line2D([], [], ls="none", marker="o", color=hs.CAT[[0, 3, 2, 4, 5, 7, 6][k]],
                              label=f"cluster {k + 1}: {e[m].sum():.1f} GeV"))
    draw_envelopes(ax, view)
    handles = draw_truth(ax, d, view) + handles
    finish(ax, ulim, vlim, handles, d, iev,
           ["area ∝ calib. energy" + (", prompt" if key == "calibP" else "")])
    axis_labels(ax, view)
    hs.save(fig, os.path.join(outdir, f"evt{iev}_{view}_clusters" + ("_prompt" if key == "calibP" else "")))


CAPTIONS = """# How to read the calibration-hit event displays

Every plot is a 2D projection of one event. r-z: z along the beam, r = transverse
radius (all phi collapsed). x-y: looking along the beam (all z collapsed). Each
pixel sums the cells whose centres fall in it. The pink line is the true track of
the gun particle (vertex, its own tracker SimHits, end point = where it first
interacted or stopped); teal lines are other saved MC particles. Grey outlines are
the calorimeter envelopes (ECal 1857-2124.5 mm, HCal 2126-4113.5 mm in r; endcaps
beyond |z| = 2307 mm).

- **overview_rz**: whole detector, r signed by the sign of y. Where the energy went.
- **calib**: calibration-hit energy = the true energy deposited in the full cell
  volume (sensor + absorber + PCB/G10/air), visible + invisible. This is the truth
  reference: the quantity a calibrated cluster should reproduce. "total" is the sum
  over all four calorimeters.
- **active**: the energy deposited in the sensors only (0.5 mm Si / 3 mm
  scintillator, Birks-quenched in the HCal) = what the detector measures before
  digitisation. Same pixels as calib, but about 40x (ECal) to 30-60x (HCal) less
  energy, and many calib pixels have no active energy at all.
- **ratio**: calib / active per pixel = the local factor that turns measured energy
  into true deposited energy. EM showers sit close to a constant (about 40 in the
  ECal; the digitiser uses 40.7 for the ECal and 58.4 for the HCal). Hadronic
  deposits scatter widely around it (invisible energy, Birks quenching, slow
  protons/fragments). That spread is what a hadronic cell/cluster calibration has
  to learn. Pixels without active energy are blank: they have calib energy that no
  calibration of measured cells can recover.
- **invisible**: the invisible part of the calibration energy (nuclear binding,
  neutrinos); positive part shown, net value in the label. It marks where nuclear
  interactions happened and is zero in pure EM showers.
- **emfrac**: fraction of the visible deposited energy made by e+-/gamma (from the
  CalibEM collection). Dark = EM sub-showers (pi0 -> gamma gamma), light = hadronic
  (charged hadrons, protons, nuclear fragments). This is the ATLAS-style EM/hadronic
  classification input.
- **calib_prompt / late**: the calibration energy inside / outside the digitiser
  time windows (t - |r_cell|/c within -0.5..1 ns in the ECal, -0.5..10 ns in the
  HCal), i.e. what the reconstructed cells can integrate vs what they cannot.
  Late energy comes from slow neutrons, captures and decays, mostly in hadronic
  showers. **clusters_prompt** is the cluster view with prompt energy, the
  target for a calibration of the reconstructed clusters.
- **clusters**: every cell with calibration energy (marker area ~ energy), coloured
  by its fate in reconstruction: in a Pandora cluster (contributes to E_dep^clus of
  that cluster), active but not clustered (below threshold or not attached), or no
  active hit (never digitised). The legend gives the calibration energy in each
  category; together they sum to the "total" of the calib plot.
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("-e", "--events", type=int, nargs="+", default=[0])
    ap.add_argument("-o", "--outdir", default="event_displays")
    args = ap.parse_args()
    hs.apply()
    events = Reader(args.input).get("events")
    for iev in args.events:
        ev = events[iev]
        d = extract(ev)
        energy_map(d, iev, "rzs", "calib", args.outdir, "Calibration-hit energy [GeV]",
                   (-5600, 5600), (-4800, 6400), nb=200, tag="overview", overview=True)
        for view in ("rz", "xy"):
            ulim, vlim = zoom_window(d, view)
            energy_map(d, iev, view, "calib", args.outdir, "Calibration-hit energy [GeV]", ulim, vlim)
            energy_map(d, iev, view, "active", args.outdir, "Active-layer energy [GeV]", ulim, vlim)
            energy_map(d, iev, view, "invisible", args.outdir, "Invisible energy [GeV]", ulim, vlim)
            emfrac_map(d, iev, view, args.outdir, ulim, vlim)
            ratio_map(d, iev, view, args.outdir, ulim, vlim)
            cluster_map(d, iev, view, args.outdir, ulim, vlim)
            if d["calibP"]:
                energy_map(d, iev, view, "calibP", args.outdir, "Prompt calibration energy [GeV]", ulim, vlim,
                           tag="calib_prompt")
                energy_map(d, iev, view, "late", args.outdir, "Late calibration energy [GeV]", ulim, vlim)
                cluster_map(d, iev, view, args.outdir, ulim, vlim, key="calibP")
        print(f"event {iev}: {label_lines(d, iev)[0]}; {len(d['calib'])} calib cells, "
              f"{len(d['active'])} active cells, {len(d['clusters'])} clusters")
    with open(os.path.join(args.outdir, "CAPTIONS.md"), "w") as fh:
        fh.write(CAPTIONS)
    print("written to", args.outdir)


if __name__ == "__main__":
    main()
