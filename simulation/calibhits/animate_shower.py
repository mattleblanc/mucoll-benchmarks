#!/usr/bin/env python3
"""Animated event display: the shower developing in time (house style, GIF).

    python3 animate_shower.py sim.edm4hep.root -e 0 -o animations/

Needs a sample simulated with per-step calibration-hit contributions
(StoreCalibContributions=True, MergeCalibContributions=False): every Geant4 step
in the calorimeters then keeps its own time and position, so the energy can be
accumulated in time order.

Per event and view (rz, xy) it writes
  evt<N>_<view>_shower.gif          the animation
  evt<N>_<view>_shower_t<k>.pdf/png a few key frames (one plot per file)

Each frame shows the calibration-hit energy (visible + invisible, at the step
positions, so finer than the cells) deposited up to time t, the true track up to
t, and the deposits made since the previous frame (the shower front).
Time runs linearly through the prompt shower, then logarithmically to 100 us;
the last frame adds everything later (e.g. neutron capture, radioactive decay).
"""
import argparse
import math
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import house_style as hs  # noqa: E402
import event_display as ed  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from podio.root_io import Reader  # noqa: E402

C_MM_PER_NS = 299.792458
NEW_COLOR = hs.CAT[2]   # gold


def deposits(tree, iev):
    """All calorimeter step contributions of one event: x, y, z [mm], t [ns], E [GeV]."""
    cols = {k: [] for k in ("x", "y", "z", "t", "e")}
    for c in ed.CALOS:
        k = f"{c}CollectionCalibContributions"
        if f"{k}/{k}.time" not in tree:
            continue
        get = lambda f: tree[f"{k}/{k}.{f}"].array(library="np", entry_start=iev, entry_stop=iev + 1)[0]  # noqa: E731
        cols["x"].append(get("stepPosition.x"))
        cols["y"].append(get("stepPosition.y"))
        cols["z"].append(get("stepPosition.z"))
        cols["t"].append(get("time"))
        cols["e"].append(get("energy"))
    d = {k: np.concatenate(v) for k, v in cols.items()}
    if len(d["t"]) == 0:
        raise SystemExit("no calibration-hit contributions: simulate with StoreCalibContributions=True")
    o = np.argsort(d["t"])
    return {k: v[o] for k, v in d.items()}


def track_times(traj):
    """Time along the true track: assume v = c (fine for the energies here)."""
    seg = np.linalg.norm(np.diff(traj, axis=0), axis=1)
    return np.concatenate([[0.0], np.cumsum(seg)]) / C_MM_PER_NS


def frame_times(t, e):
    # Linear through the prompt shower (until 80% of the energy, at most 15 ns), then
    # logarithmic to 100 us (slow neutrons, captures); a last frame adds everything.
    t0 = t[0] - 0.2
    c = np.cumsum(np.clip(e, 0, None))
    t80 = t[np.searchsorted(c, 0.80 * c[-1])]
    prompt = min(max(t80 - t0, 0.5) + 0.5, 15.0)
    lin = np.linspace(t0, t0 + prompt, 70)
    log = np.geomspace(t0 + prompt, 1e5, 31)[1:]
    return np.concatenate([lin, log, [np.inf]])


def window(dep, view, traj):
    u, v = ed.proj(np.c_[dep["x"], dep["y"], dep["z"]], view)
    w = np.clip(dep["e"], 0, None)

    def wq(a, q):
        o = np.argsort(a)
        cc = np.cumsum(w[o]) / w.sum()
        return a[o][np.searchsorted(cc, q)]
    ulo, uhi, vlo, vhi = wq(u, 0.02), wq(u, 0.98), wq(v, 0.02), wq(v, 0.98)
    tu, tv = ed.proj(traj[-1:], view)
    ulo, uhi, vlo, vhi = min(ulo, tu[0]), max(uhi, tu[0]), min(vlo, tv[0]), max(vhi, tv[0])
    half = max(uhi - ulo, vhi - vlo, 200.0) * 0.5 * 1.15
    uc, vc = 0.5 * (ulo + uhi), 0.5 * (vlo + vhi)
    head = 0.75 * half
    hn = half + 0.5 * head
    return (uc - hn, uc + hn), (vc - half, vc + half + head)


def animate(dep, truth, iev, view, outdir, nbins=140, keyframes=(0.25, 0.5, 0.75, 1.0)):
    traj = truth["traj"]
    ttraj = track_times(traj)
    ulim, vlim = window(dep, view, traj)
    ub, vb = np.linspace(*ulim, nbins + 1), np.linspace(*vlim, nbins + 1)
    u, v = ed.proj(np.c_[dep["x"], dep["y"], dep["z"]], view)
    iu, iv = np.digitize(u, ub) - 1, np.digitize(v, vb) - 1
    inside = (iu >= 0) & (iu < nbins) & (iv >= 0) & (iv < nbins)
    e = np.clip(dep["e"], 0, None)
    times = frame_times(dep["t"], dep["e"])
    stops = np.searchsorted(dep["t"], times, side="right")
    final = np.zeros((nbins, nbins))
    np.add.at(final, (iu[inside], iv[inside]), e[inside])
    vmax = final.max()

    fig, ax, cax = hs.mapfig()
    grid = np.full((nbins, nbins), np.nan)
    mesh = ax.pcolormesh(ub, vb, grid.T, cmap=hs.density_cmap, norm=LogNorm(vmin=vmax * 1e-4, vmax=vmax),
                         rasterized=True, zorder=2)
    fig.colorbar(mesh, cax=cax, label="Calibration-hit energy [GeV]")
    ed.draw_envelopes(ax, view)
    tu, tv = ed.proj(traj, view)
    (trk,) = ax.plot([], [], color=ed.TRAJ_COLOR, lw=1.8, zorder=5)
    (head,) = ax.plot([], [], marker="o", ms=6, color=ed.TRAJ_COLOR, mec="white", mew=0.8, ls="none", zorder=6)
    new = ax.scatter([], [], s=3, color=NEW_COLOR, lw=0, zorder=4)
    ax.set_xlim(*ulim)
    ax.set_ylim(*vlim)
    ax.set_aspect("equal", adjustable="box")
    ed.axis_labels(ax, view)
    hs.legend(ax, [Line2D([], [], color=ed.TRAJ_COLOR, lw=1.8, marker="o", ms=5, mec="white"),
                   Line2D([], [], ls="none", marker="o", ms=4, color=NEW_COLOR)],
              ["true track", "new deposits"])
    base = ed.label_lines({"truth": truth}, iev)
    texts = []
    acc = np.zeros((nbins, nbins))
    state = {"last": -1}

    def draw(fi):
        lo, hi = (stops[fi - 1] if fi > 0 else 0), stops[fi]
        if fi == 0:
            acc[:] = 0.0
        sel = slice(lo, hi)
        if fi != state["last"]:          # held frames repeat fi: don't add the slice twice
            m = inside[sel]
            np.add.at(acc, (iu[sel][m], iv[sel][m]), e[sel][m])
        mesh.set_array(np.where(acc > 0, acc, np.nan).T.ravel())
        new.set_offsets(np.c_[u[sel], v[sel]] if hi > lo else np.zeros((0, 2)))
        t = times[fi]
        if np.isfinite(t):
            k = np.searchsorted(ttraj, t)
            if k >= len(ttraj):
                pu, pv = tu, tv
            else:
                f = (t - ttraj[k - 1]) / max(ttraj[k] - ttraj[k - 1], 1e-12) if k > 0 else 0.0
                pu = np.r_[tu[:k], tu[k - 1] + f * (tu[k] - tu[k - 1])] if k > 0 else tu[:1]
                pv = np.r_[tv[:k], tv[k - 1] + f * (tv[k] - tv[k - 1])] if k > 0 else tv[:1]
            tlabel = (f"t = {t:.2f} ns" if t < 100 else f"t = {t:.0f} ns" if t < 1e4
                      else f"t = {t / 1e3:.0f} μs")
        else:
            pu, pv, tlabel = tu, tv, "all times (incl. > 100 μs)"
        trk.set_data(pu, pv)
        head.set_data(pu[-1:], pv[-1:])
        for tx in texts:
            tx.remove()
        texts.clear()
        before = set(ax.texts)
        # signed sum (negative invisible energy, e.g. n capture, is only clipped for the colour map)
        hs.labels(ax, "MAIA simulation", base + [tlabel, f"deposited: {dep['e'][:hi].sum():.2f} GeV"])
        texts.extend(t_ for t_ in ax.texts if t_ not in before)
        state["last"] = fi
        return [mesh, trk, head, new] + texts

    nfr = len(times)
    frames = list(range(nfr)) + [nfr - 1] * 15     # hold the last frame
    anim = FuncAnimation(fig, draw, frames=frames, blit=False)
    os.makedirs(outdir, exist_ok=True)
    gif = os.path.join(outdir, f"evt{iev}_{view}_shower.gif")
    anim.save(gif, writer=PillowWriter(fps=12), dpi=90)
    # key frames as stand-alone plots
    cum = np.cumsum(e)
    dep_at = np.where(stops > 0, cum[np.maximum(stops, 1) - 1], 0.0)   # energy deposited by each frame
    for j, q in enumerate(keyframes):
        fi = nfr - 1 if q >= 1.0 else min(int(np.searchsorted(dep_at, q * cum[-1])), nfr - 1)
        for k in range(fi + 1):
            draw(k)
        hs.save(fig, os.path.join(outdir, f"evt{iev}_{view}_shower_t{j}"), close=False)
    import matplotlib.pyplot as plt
    plt.close(fig)
    return gif


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("-e", "--events", type=int, nargs="+", default=[0])
    ap.add_argument("-o", "--outdir", default="animations")
    ap.add_argument("--views", nargs="+", default=["rz", "xy"])
    args = ap.parse_args()
    hs.apply()
    tree = uproot.open(args.input)["events"]
    events = Reader(args.input).get("events")
    for iev in args.events:
        truth = ed.extract(events[iev])["truth"]
        dep = deposits(tree, iev)
        for view in args.views:
            gif = animate(dep, truth, iev, view, args.outdir)
            print(f"event {iev} {view}: {len(dep['t'])} steps, t = {dep['t'][0]:.2f} ... "
                  f"{dep['t'][-1]:.3g} ns -> {gif}")


if __name__ == "__main__":
    main()
