#!/usr/bin/env python3
"""
Analysis of the experiment grid: for each scenario of SCENARIOS (the same list,
in the same order, as run_all.sh) the 42 subsets of at least 3 maps with the per-face estimator
(output/SCENARIO/per_face/sub_*/) and the 20 triplets with method 2 of
Roney et al. 2019 (output/SCENARIO/roney/sub_i-j-k/), each with estimate.ply
(format D) and estimate.csv.

Usage, from the repository root:
    python analysis/analysis.py
Writes everything to output/figures/.

load() is the only place where the errors are computed, the same for the two
methods: estimate.ply against the ground truth data/reference/SCENARIO/truth.ply,
both in format D (CV_L = exp(mu + alpha), r = CV_L / CV_T = exp(2 alpha)).
Errors per estimable face (mu of the estimate finite); AREA-WEIGHTED statistics:
  angle  acos(min(1, |f_est . f_true|)), degrees
  cv_l   100 (CV_L_est - CV_L_true) / CV_L_true, %
  ratio  100 (r_est - r_true) / r_true, %
q50, q90 and q99 (the values below which 50, 90, 99% of the area lies) are
computed on the absolute value; for cv_l and ratio there is also the signed
median. cond = sigma_max / sigma_min of the system of a face (per-face
estimator only, from estimate.csv); its quantiles, area-weighted, are over the
faces with computed cond (estimable, and discarded for lambda2 <= 1e-8 lambda1).

per_face.csv, one row per scenario and subset of the per-face method: non-estimable
area, errors on the estimable faces of the subset and q50, q90 and q99
of cond.

comparison.csv, one row per scenario, triplet and method (per_face, roney): non-estimable
area of the method alone, paired area (estimable by both
methods with that triplet) and errors on the paired faces. The rows with sub =
pooled gather the paired faces of the 20 triplets of the scenario; their
areas are the means over the triplets.

Boxes of plots 0a and 0b: from q25 to q75, line at q50, whiskers from q10 to q90,
area-weighted. Plot 0 only on SCENARIO_COND (uniform): cond depends on r, so
it is compared only at equal r.
plot0_cond_errors.png, faces with computed cond of the 42
subsets pooled, in 10 equal-area cond bins, each at its
median cond. Three panels angle, cv_l, ratio, one box per bin on its
estimable faces, log scale on both axes, in the title the Spearman
correlation between cond and error, unweighted; in a 2 x 2 square with the
area panel: the % of area of the bin where the estimate is not an ellipse, i.e. the
faces with cond but with estimate NaN (the smaller eigenvalue lambda2 <= 1e-8 lambda1:
the estimated conic is a hyperbola, or an ellipse flat to numerical zero). The faces
discarded because the smallest singular value is below threshold have no cond and do
not appear.
plot0_maps_cond.png, one box of cond per k, on the faces
with computed cond of all the subsets with k maps, log scale.
plot0_cond_area.png, one point per subset from per_face.csv:
cond_q90 (log) on the horizontal axis, non-estimable area on the vertical axis, color and
marker by k, in the title Spearman's rho over the 42 points.

plot5_ratio.png, only on SCENARIO_RATIO (painted): errors and non-estimable area as a
function of the true r = CV_L / CV_T (column r_true of load()). The faces are split
once, on the ground truth, into 10 equal-area bins of true r, the same for every series,
each at its area-weighted median r. Three series, one Okabe-Ito color each: per-face
with the six maps on its estimable faces; per-face and Roney on the 20 triplets pooled
(each face once per triplet, as the pooled rows of comparison.csv), on the faces paired
with that triplet. 2 x 2 panels as plot 0a: angle, cv_l, ratio with q50 (solid line)
and q90 (dashed) of the absolute value per bin, log vertical scale; and the
non-estimable area, % of the area of the bin, of each series alone (on the triplets
each method with its own lost area, not the paired one), linear from 0. The numbers
are in ratio_bins.csv, one row per bin and series.

The other figures, each in four PNGs: _area (non-estimable area, one panel) and
_angle, _cv_l, _ratio. In plot 1 the latter have a single
panel with q50, q90 and q99 (circle, square, triangle, Okabe-Ito colors),
points only without lines, with the vertical axis on a log scale and labels in
plain numbers; in plot 4 three panels q50, q90, q99, each with its
own scale. The other vertical axes, including the _area files, are linear from 0;
in plot 4 also the horizontal axis, with the same limits.
  plot1_maps_SCENARIO  one per scenario, as a function of k: one point per
                       subset and quantile, staggered within each k
  plot4_roney          one point per triplet and scenario, Roney on the horizontal axis and
                       per-face on the vertical axis, same limits on both axes,
                       one color (Okabe-Ito) and one marker per scenario, in the order of
                       the scenarios; on the axes only 0, minimum, centroid (mean of all
                       the points, cross) and maximum of the cloud, with dotted
                       projections and without grid
"""

import functools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.ticker import FormatStrFormatter, LogLocator, NullFormatter
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]                # repository root
sys.path.insert(0, str(ROOT / "estimator"))
from roney import read_ply                                # the PLY reader

REF = ROOT / "data" / "reference"
OUT = ROOT / "output"
FIG = OUT / "figures"
SCENARIOS = ["uniform", "painted"]                        # as in run_all.sh: order of tables and plots 1, 4
SCENARIO_COND = "uniform"                                 # plot 0
SCENARIO_RATIO = "painted"                                # plot 5
ERRORS = ("angle", "cv_l", "ratio")
QUANTILES = (50, 90, 99)
YLABEL = {"area": "non-estimable area [%]",
          "angle": "fiber angle error [deg]",
          "cv_l": "|relative error of CV_L| [%]",
          "ratio": "|relative error of CV_L / CV_T| [%]"}
NAMES = {"angle": "fiber angle", "cv_l": "CV_L", "ratio": "CV_L / CV_T"}   # panel titles of plot 0a
METHODS = {"per_face": "per-face", "roney": "Roney"}
BLUE = "#2a78d6"
QSTYLE = {f"q{q}": (c, m) for q, c, m in zip(QUANTILES, ("#0072B2", "#D55E00", "#009E73"), "os^")}  # Okabe-Ito
OKABE_ITO = ("#000000", "#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7")  # scenarios of plot 4
MARKERS = "osD^v<>"                                       # scenarios of plot 4


def wquantile(values, weights, q):
    """Weighted quantile. Ignores NaNs."""
    ok = np.isfinite(values) & np.isfinite(weights)
    v, w = np.asarray(values)[ok], np.asarray(weights)[ok]
    if v.size == 0:
        return np.nan
    order = np.argsort(v)
    v, w = v[order], w[order]
    cw = np.cumsum(w)
    return np.interp(q * cw[-1], cw - 0.5 * w, v)


def field(face):
    """Fiber, CV_L and r = CV_L / CV_T per face, from the face element of a format D PLY."""
    f = np.stack([face[k] for k in ("fx", "fy", "fz")], axis=1)
    return f, np.exp(face["mu"] + face["alpha"]), np.exp(2 * face["alpha"])


@functools.cache
def truth(scenario):
    """Area of each face (mm^2) and true field of a scenario, from data/reference/SCENARIO/truth.ply."""
    ply = read_ply(REF / scenario / "truth.ply")
    V = np.stack([ply["vertex"][k] for k in "xyz"], axis=1)
    F = np.stack([ply["face"][k] for k in "abc"], axis=1)
    area = np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1) / 2
    return area, field(ply["face"])


def load(scenario, method, sub):
    """One run, output/SCENARIO/METHOD/SUB/estimate.{ply,csv}, against the ground truth: DataFrame
    indexed by face_id with area, estimable, cond (per-face) or k_used (Roney) and the errors angle,
    cv_l, ratio (degrees, %, %), NaN on the non-estimable faces. The only place with the
    formulas of the errors."""
    area, (f_true, cvl_true, r_true) = truth(scenario)
    case = OUT / scenario / method / sub
    face = read_ply(case / "estimate.ply")["face"]
    f, cvl, r = field(face)
    df = pd.read_csv(case / "estimate.csv").set_index("face_id")
    df["area"] = area
    df["estimable"] = np.isfinite(face["mu"])
    df["r_true"] = r_true
    # |a.b|: angle between directions, unoriented
    df["angle"] = np.degrees(np.arccos(np.minimum(1, np.abs(np.sum(f * f_true, axis=1)))))
    df["cv_l"] = 100 * (cvl - cvl_true) / cvl_true
    df["ratio"] = 100 * (r - r_true) / r_true
    return df


def non_estimable(df):
    """% of non-estimable area."""
    w = df["area"]
    return 100 * w[~df["estimable"]].sum() / w.sum()


def quantiles(df):
    """q50, q90 and q99 of the absolute value of each error on the estimable faces, then the signed
    median of cv_l and ratio."""
    est = df[df["estimable"]]
    w = est["area"]
    out = {f"{e}_q{q}": wquantile(est[e].abs(), w, q / 100) for e in ERRORS for q in QUANTILES}
    for e in ("cv_l", "ratio"):
        out[f"{e}_median"] = wquantile(est[e], w, 0.50)
    return out


def cond_quantiles(df):
    """q50, q90 and q99 of cond on the faces with computed cond."""
    c = df[df["cond"].notna()]
    return {f"cond_q{q}": wquantile(c["cond"], c["area"], q / 100) for q in QUANTILES}


def subs(scenario, method):
    """The subsets of a run, by number of maps and then by name."""
    return sorted((d.name for d in (OUT / scenario / method).glob("sub_*")), key=lambda s: (s.count("-"), s))


def tables():
    """per_face.csv and comparison.csv, as DataFrames."""
    per_face_table, comparison = [], []
    for sc in SCENARIOS:
        per_face = {sub: load(sc, "per_face", sub) for sub in subs(sc, "per_face")}
        roney = {sub: load(sc, "roney", sub) for sub in subs(sc, "roney")}

        for sub, df in per_face.items():
            per_face_table.append({"scenario": sc, "sub": sub, "k": sub.count("-") + 1,
                                   "non_estimable_area_%": non_estimable(df), **quantiles(df), **cond_quantiles(df)})

        for method, runs in (("per_face", per_face), ("roney", roney)):
            rows, paired = [], []
            for sub in roney:
                df = runs[sub]
                both = roney[sub]["estimable"] & per_face[sub]["estimable"]
                rows.append({"scenario": sc, "sub": sub, "method": method,
                             "non_estimable_area_%": non_estimable(df),
                             "paired_area_%": 100 * df["area"][both].sum() / df["area"].sum(), **quantiles(df[both])})
                paired.append(df[both])
            comparison += rows
            comparison.append({"scenario": sc, "sub": "pooled", "method": method,
                               "non_estimable_area_%": np.mean([r["non_estimable_area_%"] for r in rows]),
                               "paired_area_%": np.mean([r["paired_area_%"] for r in rows]),
                               **quantiles(pd.concat(paired))})
    return pd.DataFrame(per_face_table), pd.DataFrame(comparison)


def figure(name, title, draw, together=False, after=None):
    """Four PNGs, output/figures/NAME_{area,angle,cv_l,ratio}.png. draw(ax, cols, ylabel) draws
    in one panel the columns of cols, {label: column}. With together q50, q90 and q99 are in a
    single panel, otherwise one panel per quantile. after(fig), if given, is called before saving."""
    for q, ylabel in YLABEL.items():
        if q == "area":
            panels = {"non-estimable area": {"non-estimable area": "non_estimable_area_%"}}
        elif together:
            panels = {"q50, q90 and q99": {f"q{p}": f"{q}_q{p}" for p in QUANTILES}}
        else:
            panels = {f"q{p}": {f"q{p}": f"{q}_q{p}"} for p in QUANTILES}
        fig, axes = plt.subplots(1, len(panels), figsize=(max(4.6 * len(panels), 7.5), 4.6), squeeze=False,
                                 layout="constrained")
        for ax, (panel, cols) in zip(axes[0], panels.items()):
            ax.grid(alpha=0.3)
            draw(ax, cols, ylabel)
            ax.set_title(panel)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside lower center", ncol=min(len(labels), 6 if len(panels) > 1 else 4),
                   frameon=False)
        fig.suptitle(title)
        if after:
            after(fig)
        fig.savefig(FIG / f"{name}_{q}.png", dpi=140, bbox_inches="tight")
        plt.close(fig)


def log_scale(ax, which):
    """Log scale on the given axes ("x", "y"), with labels in plain numbers (1, 2, 5, 10, 20) instead
    of powers of 10: ticks 1-2-5 per decade, or 1-9 if with those there are fewer than four labels.
    Must be called after drawing the data."""
    for a in which:
        getattr(ax, f"set_{a}scale")("log")
        axis = getattr(ax, f"{a}axis")
        lo, hi = getattr(ax, f"get_{a}lim")()
        for subs in ((1, 2, 5), (1, 2, 3, 4, 5, 6, 7, 8, 9)):
            ticks = [t for t in LogLocator(subs=subs).tick_values(lo, hi) if lo <= t <= hi]
            if len(ticks) >= 4:
                break
        axis.set_major_locator(LogLocator(subs=subs))
        axis.set_major_formatter(FormatStrFormatter("%g"))
        axis.set_minor_formatter(NullFormatter())
        ax.grid(axis=a, which="minor", alpha=0.15)


def separate_tick_labels(fig):
    """Moves overlapping tick labels onto several rows: downward on the x axis, to the
    left on the y axis."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px_pt = 72 / fig.dpi
    for ax in fig.axes:
        for axis, side in ((ax.xaxis, 0), (ax.yaxis, 1)):
            ticks = sorted(axis.get_major_ticks(), key=lambda t: t.get_loc())
            boxes = [t.label1.get_window_extent(renderer) for t in ticks]
            step = max((b.height if side == 0 else b.width) for b in boxes) + 3
            ends = []                                     # per row, where the last label ends
            for t, b in zip(ticks, boxes):
                a, z = (b.x0, b.x1) if side == 0 else (b.y0, b.y1)
                r = next((i for i, f in enumerate(ends) if a > f + 3), len(ends))
                ends[r:r + 1] = [z]
                t.set_pad(t.get_pad() + r * step * px_pt)


def y_axis(ax, ylabel, log=False):
    """Vertical axis: on a log scale, or linear from 0."""
    if log:
        ax.set_ylabel(f"{ylabel} (log scale)")
        log_scale(ax, "y")
    else:
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)


def points(ax, x, df, cols):
    """One point per row and column of cols, color and marker of the quantile."""
    for label, col in cols.items():
        color, marker = QSTYLE.get(label, (BLUE, "o"))
        ax.scatter(x, df[col], s=24, color=color, marker=marker, label=label, zorder=3)


def boxes(ax, x, x0, x1, groups, col):
    """One box per group, in x from x0 to x1: from q25 to q75 of the absolute value of col, line at q50 and whiskers
    from q10 to q90, area-weighted."""
    q10, q25, q50, q75, q90 = (groups.apply(lambda s: wquantile(s[col].abs(), s["area"], q)).to_numpy()
                              for q in (0.10, 0.25, 0.50, 0.75, 0.90))
    blue = QSTYLE["q50"][0]
    ax.vlines(x, q10, q25, color=blue, lw=1)
    ax.vlines(x, q75, q90, color=blue, lw=1)
    ax.bar(x0, q75 - q25, width=x1 - x0, bottom=q25, align="edge", facecolor=to_rgba(blue, 0.25), edgecolor=blue,
           lw=1, zorder=3)
    ax.hlines(q50, x0, x1, color=blue, lw=2, zorder=4)


def boxes_note(fig):
    """The legend of the boxes, below the figure."""
    fig.supxlabel("box: q25–q75 · line: q50 · whiskers: q10–q90 (area-weighted)", fontsize="medium")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    per_face_table, comparison = tables()
    per_face_table.to_csv(FIG / "per_face.csv", index=False, float_format="%.6g")
    comparison.to_csv(FIG / "comparison.csv", index=False, float_format="%.6g")

    # plot 0a: SCENARIO_COND, all subsets, error as a function of cond
    df = pd.concat([load(SCENARIO_COND, "per_face", sub).assign(k=sub.count("-") + 1)
                    for sub in subs(SCENARIO_COND, "per_face")])
    c = df[df["cond"].notna()].sort_values("cond")
    w = c["area"]
    c["bin"] = np.minimum((10 * (w.cumsum() - 0.5 * w) / w.sum()).astype(int), 9)   # 10 equal-area bins
    est = c[c["estimable"]]
    xc = c.groupby("bin").apply(lambda s: wquantile(s["cond"], s["area"], 0.50)).to_numpy()
    h = 10 ** (0.4 * np.diff(np.log10(xc)).min())          # half width of the boxes, constant in log
    fig, axes = plt.subplots(2, 2, figsize=(10, 9), layout="constrained")
    axes = axes.flat
    for ax, e in zip(axes, ERRORS):
        boxes(ax, xc, xc / h, xc * h, est.groupby("bin"), e)
        ax.grid(alpha=0.3)
        ax.set_xlabel("cond (log scale)")
        ax.set_ylabel(f"{YLABEL[e]} (log scale)")
        log_scale(ax, "xy")
        ax.set_title(f"{NAMES[e]}, Spearman's ρ = {spearmanr(est['cond'], est[e].abs()).statistic:.2f}")
    ax = axes[3]
    ax.plot(xc, c.groupby("bin").apply(non_estimable), color=QSTYLE["q50"][0], lw=1.5, marker="o", ms=5)
    ax.grid(alpha=0.3)
    ax.set_xlabel("cond (log scale)")
    log_scale(ax, "x")
    y_axis(ax, "area where the estimate is not an ellipse [%]")
    ax.set_title("estimate is not an ellipse")
    boxes_note(fig)
    fig.suptitle(f"{SCENARIO_COND}, all map subsets: error as a function of cond,\nequal-area bins")
    fig.savefig(FIG / "plot0_cond_errors.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # plot 0b: SCENARIO_COND, one box of cond per k
    ks = np.array(sorted(c["k"].unique()))
    fig, ax = plt.subplots(figsize=(7.5, 4.6), layout="constrained")
    boxes(ax, ks, ks - 0.25, ks + 0.25, c.groupby("k"), "cond")
    ax.grid(alpha=0.3)
    ax.set_xticks(ks)
    ax.set_xlabel("maps used, k")
    y_axis(ax, "cond", log=True)
    boxes_note(fig)
    fig.suptitle(f"{SCENARIO_COND}: distribution of cond by number of maps")
    fig.savefig(FIG / "plot0_maps_cond.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # plot 0c: SCENARIO_COND, non-estimable area vs cond_q90, one point per subset
    g0 = per_face_table[per_face_table["scenario"] == SCENARIO_COND]
    rho = spearmanr(g0["cond_q90"], g0["non_estimable_area_%"]).statistic
    fig, ax = plt.subplots(figsize=(7.5, 4.6), layout="constrained")
    for (k, s), (color, marker) in zip(g0.groupby("k"), (("#000000", "o"), ("#0072B2", "s"), ("#D55E00", "^"),
                                                         ("#009E73", "D"))):
        ax.scatter(s["cond_q90"], s["non_estimable_area_%"], s=28, color=color, marker=marker, label=f"k = {k}",
                   zorder=3)
    ax.grid(alpha=0.3)
    ax.set_xlabel("cond, subset q90 (log scale)")
    log_scale(ax, "x")
    y_axis(ax, YLABEL["area"])
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center", ncol=4, frameon=False)
    fig.suptitle(f"{SCENARIO_COND}: non-estimable area and tail of cond, one point per map subset "
                 f"(Spearman's ρ = {rho:.2f})")
    fig.savefig(FIG / "plot0_cond_area.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # plot 1: for each scenario, as a function of k, points staggered within each k
    for sc in SCENARIOS:
        g1 = per_face_table[per_face_table["scenario"] == sc].copy()
        i, n = g1.groupby("k").cumcount(), g1.groupby("k")["k"].transform("size")
        g1["x"] = g1["k"] + 0.6 * (i - (n - 1) / 2) / n

        def draw1(ax, cols, ylabel):
            points(ax, g1["x"], g1, cols)
            ax.set_xticks(sorted(g1["k"].unique()))
            ax.set_xlabel("maps used, k")
            y_axis(ax, ylabel, log=len(cols) > 1)

        figure(f"plot1_maps_{sc}", f"{sc}: one point per map subset", draw1, together=True)

    # plot 4: per-face vs Roney, one point per triplet and scenario
    triplets = comparison[comparison["sub"] != "pooled"]
    g4 = triplets[triplets["method"] == "roney"].merge(triplets[triplets["method"] == "per_face"],
                                                       on=["scenario", "sub"],
                                                       suffixes=("_roney", "_per_face"))

    def draw4(ax, cols, ylabel):
        (col,) = cols.values()
        for sc, color, marker in zip(SCENARIOS, OKABE_ITO, MARKERS):
            s = g4[g4["scenario"] == sc]
            ax.scatter(s[f"{col}_roney"], s[f"{col}_per_face"], s=22, color=color, marker=marker, label=sc,
                       zorder=3)
        x, y = g4[f"{col}_roney"], g4[f"{col}_per_face"]
        lim = (0, 1.05 * max(x.max(), y.max()))
        ax.plot(lim, lim, color="#555", lw=1, label="y = x")
        ax.set_xlim(*lim)
        ax.set_ylim(*lim)
        # on the axes only 0, minimum, centroid and maximum of the cloud, with the projections
        mx, my = x.mean(), y.mean()
        projection = dict(color="#444", lw=0.8, ls=":", zorder=2)
        for i in (x.idxmin(), x.idxmax()):
            ax.plot([x[i], x[i]], [0, y[i]], **projection)
        for i in (y.idxmin(), y.idxmax()):
            ax.plot([0, x[i]], [y[i], y[i]], **projection)
        ax.plot([mx, mx, 0], [0, my, my], **projection)
        ax.scatter(mx, my, s=170, marker="+", linewidths=2, color="#444", label="centroid", zorder=5)
        for setter, v in ((ax.set_xticks, x), (ax.set_yticks, y)):
            vals = [0, v.min(), v.mean(), v.max()]
            setter(vals, ["0"] + [f"{t:#.3g}".rstrip(".") for t in vals[1:]])   # 19.0, not 19
        ax.grid(False)
        ax.set_aspect("equal")
        ax.set_xlabel(f"{METHODS['roney']}: {ylabel}")
        ax.set_ylabel(f"{METHODS['per_face']}: {ylabel}")

    figure("plot4_roney", "map triplets: per-face vs Roney, paired faces (area: each method alone)\n"
           "on the axes: min, centroid (mean of all the points) and max", draw4, after=separate_tick_labels)

    plot5()


def plot5():
    """Plot 5 and ratio_bins.csv: errors and non-estimable area on SCENARIO_RATIO as a function of the
    true r, in 10 equal-area bins of true r computed once on the ground truth."""
    sc = SCENARIO_RATIO
    area, (_, _, r_true) = truth(sc)
    order = np.argsort(r_true)
    w = area[order]
    bins = np.empty(len(area), dtype=int)
    bins[order] = np.minimum((10 * (np.cumsum(w) - 0.5 * w) / w.sum()).astype(int), 9)   # 10 equal-area bins
    xc = np.array([wquantile(r_true[bins == b], area[bins == b], 0.50) for b in range(10)])

    # each series: (label, color, all the faces for the area, faces for the errors)
    all_maps = load(sc, "per_face", subs(sc, "per_face")[-1])
    per_face, roney, paired_pf, paired_ro = [], [], [], []
    for sub in subs(sc, "roney"):
        pf, ro = load(sc, "per_face", sub), load(sc, "roney", sub)
        both = pf["estimable"] & ro["estimable"]
        per_face.append(pf)
        roney.append(ro)
        paired_pf.append(pf[both])
        paired_ro.append(ro[both])
    series = [("per-face, 6 maps", "#0072B2", all_maps, all_maps[all_maps["estimable"]]),
              ("per-face, 20 triplets", "#009E73", pd.concat(per_face), pd.concat(paired_pf)),
              ("Roney, 20 triplets", "#D55E00", pd.concat(roney), pd.concat(paired_ro))]

    rows = []
    for label, _, everything, errors in series:
        b_all, b_err = bins[everything.index], bins[errors.index]
        for b in range(10):
            row = {"scenario": sc, "series": label, "bin": b, "r_median": xc[b],
                   "non_estimable_area_%": non_estimable(everything[b_all == b])}
            e = errors[b_err == b]
            for err in ERRORS:
                for q in (50, 90):
                    row[f"{err}_q{q}"] = wquantile(e[err].abs(), e["area"], q / 100)
            rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(FIG / "ratio_bins.csv", index=False, float_format="%.6g")

    fig, axes = plt.subplots(2, 2, figsize=(10, 9), layout="constrained")
    axes = axes.flat
    for ax, e in zip(axes, ERRORS):
        for label, color, _, _ in series:
            t = table[table["series"] == label]
            ax.plot(xc, t[f"{e}_q50"], color=color, lw=1.5, marker="o", ms=4, label=f"{label}, q50")
            ax.plot(xc, t[f"{e}_q90"], color=color, lw=1.5, ls="--", marker="o", ms=4, mfc="white",
                    label=f"{label}, q90")
        ax.grid(alpha=0.3)
        ax.set_xlabel("true r = CV_L / CV_T")
        y_axis(ax, YLABEL[e], log=True)
        ax.set_title(NAMES[e])
    ax = axes[3]
    for label, color, _, _ in series:
        t = table[table["series"] == label]
        ax.plot(xc, t["non_estimable_area_%"], color=color, lw=1.5, marker="o", ms=4, label=label)
    ax.grid(alpha=0.3)
    ax.set_xlabel("true r = CV_L / CV_T")
    y_axis(ax, "non-estimable area of the bin [%]")
    ax.set_title("non-estimable area (each method alone)")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle(f"{sc}: errors as a function of the true r, 10 equal-area bins of r\n"
                 "per-face with 6 maps (its estimable faces); per-face and Roney on the 20 triplets pooled "
                 "(paired faces)\nq50 solid, q90 dashed, of |error|, area-weighted")
    fig.savefig(FIG / "plot5_ratio.png", dpi=140, bbox_inches="tight")
    plt.close(fig)



if __name__ == "__main__":
    main()
