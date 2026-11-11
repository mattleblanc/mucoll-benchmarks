"""Matt LeBlanc's house plot style for matplotlib.

Copy of ~/.claude/skills/house-plot-style/house_style.py, vendored so the calibhits
scripts are self-contained.

Usage
-----
    import house_style as hs
    hs.apply()
    fig, ax = hs.squarefig()
    hs.data_points(ax, x, y, yerr, label="Data")
    ax.stairs(h, edges, color=hs.CAT[0], label="Pythia 8")
    hs.labels(ax, "MAIA simulation, 20 GeV π⁻", ["ECal endcap"])
    hs.legend(ax)
    hs.save(fig, "outdir/name")          # writes name.pdf and name.png

Note: this file was written from the rules in SKILL.md. The Julia/Makie implementation
(Alpenglow.jl + config/plot.toml in the isotropy repo) is the source of truth.
"""
import os

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm  # noqa: F401  (re-exported)
from matplotlib.ticker import AutoMinorLocator

# ---------------------------------------------------------------- colours
CAT = ["#4C3FB0", "#D6497E", "#C08512", "#3682EF", "#B34700", "#8D267F", "#00A1AA", "#177937"]
CAT_NAMES = ["violet", "rose", "gold", "blue", "orange", "plum", "teal", "green"]
LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1, 1, 1)), (0, (5, 2)), (0, (1, 1)), (0, (8, 2, 1, 2))]
RAMP = ["#86B6EF", "#5598E7", "#2A78D6", "#1C5CAB", "#104281"]
ALPENGLOW = ["#44308F", "#783997", "#B13E86", "#D6497E", "#E26A61", "#D98A3D", "#C08512"]
DENSITY = ["#FCFCFB", "#CDE2FB", "#9EC5F4", "#5598E7", "#2A78D6", "#184F95", "#0D366B"]
INK, MUTED, SURFACE = "#000000", "#52514E", "#FFFFFF"
NOT_USED = "#BDBDBD"

GENERATORS = {  # family colour + marker per member
    "Pythia 8":               ("#654EA3", "s"),
    "Dire":                   ("#654EA3", "P"),
    "Vincia":                 ("#654EA3", "X"),
    "Herwig 7 Angle-Ordered": ("#5FBBA1", (4, 1, 0)),
    "Herwig 7 Dipole":        ("#5FBBA1", (4, 1, 45)),
    "Sherpa 3":               ("#FA7268", "*"),
    "PanGlobal":              ("#D946D6", ">"),
    "PanLocal":               ("#D946D6", "<"),
    "Pythia 8, true":         ("#9E9E9E", "s"),
}

density_cmap = LinearSegmentedColormap.from_list("house_density", DENSITY)
density_cmap.set_bad(SURFACE)
density_cmap.set_under(SURFACE)

# ---------------------------------------------------------------- geometry [inches]
PANEL = 4.4          # square panel side
GAP = 0.08           # gap between stacked panels
MARGIN = dict(left=0.85, right=0.2, bottom=0.65, top=0.2)
CBAR_W, CBAR_PAD, CBAR_LABEL = 0.15, 0.1, 0.8
LEGEND_SIZE = 11
LABEL_SIZE = 12
SUBLABEL_SIZE = 10


def apply():
    """Set the house rcParams (boxed frame, inward mirrored ticks, no grid, TrueType fonts)."""
    matplotlib.rcParams.update({
        "font.size": 11, "axes.labelsize": 12, "xtick.labelsize": 11, "ytick.labelsize": 11,
        "legend.fontsize": LEGEND_SIZE, "legend.frameon": False,
        "axes.grid": False, "axes.linewidth": 0.9, "axes.edgecolor": INK,
        "axes.spines.top": True, "axes.spines.right": True,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.major.size": 6, "ytick.major.size": 6, "xtick.minor.size": 3, "ytick.minor.size": 3,
        "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "lines.linewidth": 1.6, "errorbar.capsize": 0,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.prop_cycle": matplotlib.cycler(color=CAT),
        "figure.dpi": 100, "savefig.dpi": 200,
    })


def _minor(ax):
    for axis in (ax.xaxis, ax.yaxis):
        if axis.get_scale() == "linear":
            axis.set_minor_locator(AutoMinorLocator(4))   # 4 intervals between majors
    ax.tick_params(which="both", direction="in", top=True, right=True)


def _fig(panel_heights, width=PANEL, extra_right=0.0):
    """Stack panels of the given heights [in] with fixed margins; returns fig, [axes]."""
    h = MARGIN["bottom"] + sum(panel_heights) + GAP * (len(panel_heights) - 1) + MARGIN["top"]
    w = MARGIN["left"] + width + extra_right + MARGIN["right"]
    fig = plt.figure(figsize=(w, h))
    axes, y = [], h - MARGIN["top"]
    for i, ph in enumerate(panel_heights):
        y -= ph
        ax = fig.add_axes([MARGIN["left"] / w, y / h, width / w, ph / h],
                          sharex=axes[0] if axes else None)
        axes.append(ax)
        y -= GAP
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)
    for ax in axes:
        _minor(ax)
    return fig, axes


def squarefig():
    fig, (ax,) = _fig([PANEL])
    return fig, ax


def ratiofig():
    main = 0.76 * PANEL
    fig, (ax, rax) = _fig([main, 0.32 * main])
    return fig, ax, rax


def ratiofig3():
    fig, (ax, rax, uax) = _fig([0.64 * PANEL, 0.32 * PANEL, 0.24 * PANEL])
    return fig, ax, rax, uax


def mapfig():
    """Square panel + thin colour-bar axis to its right. Returns fig, ax, cax."""
    fig, (ax,) = _fig([PANEL], extra_right=CBAR_PAD + CBAR_W + CBAR_LABEL)
    w, h = fig.get_size_inches()
    x0 = (MARGIN["left"] + PANEL + CBAR_PAD) / w
    cax = fig.add_axes([x0, MARGIN["bottom"] / h, CBAR_W / w, PANEL / h])
    cax.tick_params(direction="in", which="both")
    return fig, ax, cax


def ratio_line(rax, y=1.0):
    rax.axhline(y, color=INK, lw=0.8, ls="--", zorder=0)


def labels(ax, experiment, lines=(), x_pt=10, y_pt=-10):
    """Experiment label at top left, further lines below it, spaced in points."""
    ax.annotate(experiment, xy=(0, 1), xycoords="axes fraction", xytext=(x_pt, y_pt),
                textcoords="offset points", ha="left", va="top", fontsize=LABEL_SIZE, color=INK)
    off = y_pt - LABEL_SIZE * 1.45
    for line in lines:
        ax.annotate(line, xy=(0, 1), xycoords="axes fraction", xytext=(x_pt, off),
                    textcoords="offset points", ha="left", va="top", fontsize=SUBLABEL_SIZE, color=INK)
        off -= SUBLABEL_SIZE * 1.4


def legend(ax, handles=None, labels_=None, loc="upper right", **kw):
    """One column, top right, frameless, 11 pt. Data entries (drawn with data_points) go first.

    If it collides with the experiment label, pass e.g. bbox_to_anchor=(1, 0.85) to drop it below.
    """
    if handles is None:
        handles, labels_ = ax.get_legend_handles_labels()
        is_data = [getattr(h, "_house_data", False) for h in handles]
        order = sorted(range(len(handles)), key=lambda i: not is_data[i])
        handles, labels_ = [handles[i] for i in order], [labels_[i] for i in order]
    args = (handles, labels_)
    kw.setdefault("borderaxespad", 0.7)
    kw.setdefault("handlelength", 1.6)
    return ax.legend(*args, loc=loc, ncol=1, frameon=False, fontsize=LEGEND_SIZE, **kw)


def data_points(ax, x, y, yerr=None, xerr=None, label="Data", **kw):
    """Data: solid black circles, error bars without caps, no connecting line."""
    cont = ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt="o", color=INK, ms=4.5, lw=1.0,
                       capsize=0, label=label, zorder=5, **kw)
    cont._house_data = True
    return cont


def mc_points(ax, x, y, yerr=None, color=CAT[0], marker="s", label=None, fill="full", **kw):
    """Simulation as coloured markers with a thin white ring. fill: full | none | left."""
    mfc = {"full": color, "none": "none", "left": color}[fill]   # open markers are see-through
    return ax.errorbar(x, y, yerr=yerr, fmt=marker, color=color, mfc=mfc, mec=color if fill == "none" else "white",
                       mew=0.8 if fill != "none" else 1.2, fillstyle=fill if fill == "left" else "full",
                       ms=6, lw=1.0, capsize=0, label=label, zorder=4, **kw)


def not_used(ax, x0, x1, text="not used"):
    ax.axvspan(x0, x1, color=NOT_USED, alpha=0.35, lw=0, zorder=0)
    ax.annotate(text, xy=((x0 + x1) / 2, 0.02), xycoords=("data", "axes fraction"),
                ha="center", va="bottom", fontsize=8, color=MUTED)


def save(fig, path_noext, close=True):
    """Write path.pdf and path.png. (One plot per file: save each panel of a multipanel too.)"""
    d = os.path.dirname(path_noext)
    if d:
        os.makedirs(d, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{path_noext}.{ext}")
    if close:
        plt.close(fig)
