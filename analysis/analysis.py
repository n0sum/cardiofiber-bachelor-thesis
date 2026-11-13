#!/usr/bin/env python3
"""
Analysis of the experiment grid: for each scenario of SCENARIOS (the same list,
in the same order, as run_all.sh: MAIN_SCENARIOS, then the rest of SWEEP) the 42
subsets of at least 3 maps with the per-face estimator
(output/SCENARIO/per_face/sub_*/) and the 20 triplets with method 2 of
Roney et al. 2019 (output/SCENARIO/roney/sub_i-j-k/), each with estimate.ply
(format D) and estimate.csv.

Usage, from the repository root:
    python analysis/analysis.py
Writes everything to output/figures/.

load() is the only place where the errors are computed, the same for the two
methods: estimate.ply against the ground truth data/reference/SCENARIO/truth.ply,
both in format D (mean speed e^mu = sqrt(CV_L CV_T), r = CV_L / CV_T = exp(2 alpha)).
Errors per estimable face (mu of the estimate finite); AREA-WEIGHTED statistics:
  angle    acos(min(1, |f_est . f_true|)), degrees
  mean_cv  100 (e^mu_est - e^mu_true) / e^mu_true, %
  ratio    100 (r_est - r_true) / r_true, %
q50, q90 and q99 (the values below which 50, 90, 99% of the area lies) are
computed on the absolute value; for mean_cv and ratio there is also the signed
median. cond = sigma_max / sigma_min of the system of a face (per-face
estimator only, from estimate.csv); its quantiles, area-weighted, are over the
faces with computed cond (estimable, and discarded for lambda2 <= 1e-8 lambda1).

per_face.csv, one row per scenario and subset of the per-face method: non-estimable
area, errors on the estimable faces of the subset and q50, q90 and q99
of cond.

comparison.csv, one row per scenario, triplet and method (per_face, roney): non-estimable
area of the method alone, paired area (estimable by both methods with that triplet)
and errors on the paired faces. For the Roney rows the non-estimable area is also
split in two: fit_area_% (k_used < 3, lost to the planar fit) and conic_area_%
(k_used = 3 and estimate NaN: the conic is not an ellipse or the system is nearly
singular); empty on the per-face rows. The rows with sub = pooled gather the
paired faces of the 20 triplets of the scenario; their areas are the means over
the triplets.

Boxes (plots 0b, 0d, 3): from q25 to q75, line at q50, whiskers from q10 to q90,
area-weighted. Plots 0a, 0b, 0c only on SCENARIO_COND (uniform): cond depends on
r, so there it is compared only at equal r; plots 0d, 2 and 3 are on the scenarios
of SWEEP (uniform fields differing only in r: 1.6, 2 (uniform), 2.4), with all
the 42 subsets.
plot0_cond_errors.png (0a), faces with computed cond of the 42 subsets pooled, in
10 equal-area cond bins, each at its median cond. Three panels angle, mean_cv,
ratio: per bin q50, q90 and q99 of |error| on its estimable faces, area-weighted,
points joined by lines (colors and markers of QSTYLE), log scale on both axes,
in the title the Spearman correlation between cond and |error|, unweighted; in a
2 x 2 square with the area panel: the % of area of the bin where the estimate is
not an ellipse, i.e. the faces with cond but with estimate NaN (the smaller
eigenvalue lambda2 <= 1e-8 lambda1: the estimated conic is a hyperbola, or an
ellipse flat to numerical zero). The faces discarded because the smallest
singular value is below threshold have no cond and do not appear.
plot0_maps_cond.png (0b), one box of cond per k, on the faces with computed cond
of all the subsets with k maps, log scale.
plot0_cond_area.png (0c), one point per subset from per_face.csv: cond_q90 (log)
on the horizontal axis, non-estimable area on the vertical axis, color and marker
by k, in the title Spearman's rho over the 42 points.
plot0_ratio_cond.png (0d), one box of cond per scenario of SWEEP, r on the
horizontal axis, on the faces with computed cond of all its 42 subsets (the same
mix of k in each box), log scale.
plot3_cond_k_ratio.png, boxes of cond grouped by k (3 to 6), in each group one box
per scenario of SWEEP, one color per r (those of plot 4), on the faces with
computed cond of the subsets with that k; log scale; in the title Spearman's rho
of cond with k and with r over the faces with computed cond of the three
scenarios, unweighted. The numbers in cond_k_ratio.csv, one row per (k, r):
cond_q50, cond_q90 (area-weighted) and the number of faces.

The other figures, each in four PNGs: _area (non-estimable area, one panel) and
_angle, _mean_cv, _ratio. In plots 1 and 2 the latter have a single
panel with q50, q90 and q99 (circle, square, triangle, Okabe-Ito colors),
points only without lines, with the vertical axis on a log scale and labels in
plain numbers; in plot 4 three panels q50, q90, q99, each with its
own scale. The other vertical axes, including the _area files, are linear from 0;
in plot 4 also the horizontal axis, with the same limits.
  plot1_maps_SCENARIO  one per scenario of PLOT1_SCENARIOS (painted), as a function of k: one point per
                       subset and quantile, staggered within each k
  plot2_ratio          per-face only, the scenarios of SWEEP, as a function of r: one point
                       per subset and quantile, the 42 subsets of each r staggered
                       within it, ordered by k
  plot4_roney          one point per triplet and scenario of SCENARIOS (all four, 80
                       points), Roney on the horizontal axis and per-face on the
                       vertical axis, same limits on both axes, one color (Okabe-Ito)
                       and one marker per scenario, in the order of run_all.sh, labels
                       of SCENARIO_LABELS; limits, minimum, maximum and centroid over
                       all the points; on the axes only 0, minimum, centroid (mean of
                       all the points, cross) and maximum of the cloud, with dotted
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
MAIN_SCENARIOS = ["uniform", "painted"]
PLOT1_SCENARIOS = ["painted"]                             # plot 1 (uniform is in plot 2, at r = 2)
SWEEP = [("r1.6", 1.6), ("uniform", 2.0), ("r2.4", 2.4)]  # plots 0d, 2, 3: uniform fields differing only in r, increasing r
SCENARIOS = MAIN_SCENARIOS + [sc for sc, _ in SWEEP if sc not in MAIN_SCENARIOS]   # tables and plot 4: run_all.sh
SCENARIO_LABELS = {"uniform": "uniform (r = 2)", "painted": "painted", "r1.6": "r = 1.6", "r2.4": "r = 2.4"}
SCENARIO_COND = "uniform"                                 # plots 0a, 0b, 0c
ERRORS = ("angle", "mean_cv", "ratio")
QUANTILES = (50, 90, 99)
YLABEL = {"area": "non-estimable area [%]",
          "angle": "fiber angle error [deg]",
          "mean_cv": "|relative error of mean speed e^μ| [%]",
          "ratio": "|relative error of CV_L / CV_T| [%]"}
NAMES = {"angle": "fiber angle", "mean_cv": "mean speed e^μ", "ratio": "CV_L / CV_T"}   # panel titles of plot 0a
METHODS = {"per_face": "per-face", "roney": "Roney"}
BLUE = "#2a78d6"
QSTYLE = {f"q{q}": (c, m) for q, c, m in zip(QUANTILES, ("#0072B2", "#D55E00", "#009E73"), "os^")}  # Okabe-Ito
OKABE_ITO = ("#000000", "#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7")  # scenarios of plots 3, 4
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
    """Fiber, mean speed e^mu = sqrt(CV_L CV_T) and r = CV_L / CV_T per face, from the face element of
    a format D PLY."""
    f = np.stack([face[k] for k in ("fx", "fy", "fz")], axis=1)
    return f, np.exp(face["mu"]), np.exp(2 * face["alpha"])


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
    mean_cv, ratio (degrees, %, %), NaN on the non-estimable faces. The only place with the
    formulas of the errors."""
    area, (f_true, speed_true, r_true) = truth(scenario)
    case = OUT / scenario / method / sub
    face = read_ply(case / "estimate.ply")["face"]
    f, speed, r = field(face)
    df = pd.read_csv(case / "estimate.csv").set_index("face_id")
    df["area"] = area
    df["estimable"] = np.isfinite(face["mu"])
    # |a.b|: angle between directions, unoriented
    df["angle"] = np.degrees(np.arccos(np.minimum(1, np.abs(np.sum(f * f_true, axis=1)))))
    df["mean_cv"] = 100 * (speed - speed_true) / speed_true
    df["ratio"] = 100 * (r - r_true) / r_true
    return df


def non_estimable(df):
    """% of non-estimable area."""
    w = df["area"]
    return 100 * w[~df["estimable"]].sum() / w.sum()


def quantiles(df):
    """q50, q90 and q99 of the absolute value of each error on the estimable faces, then the signed
    median of mean_cv and ratio."""
    est = df[df["estimable"]]
    w = est["area"]
    out = {f"{e}_q{q}": wquantile(est[e].abs(), w, q / 100) for e in ERRORS for q in QUANTILES}
    for e in ("mean_cv", "ratio"):
        out[f"{e}_median"] = wquantile(est[e], w, 0.50)
    return out


def cond_quantiles(df):
    """q50, q90 and q99 of cond on the faces with computed cond."""
    c = df[df["cond"].notna()]
    return {f"cond_q{q}": wquantile(c["cond"], c["area"], q / 100) for q in QUANTILES}


def roney_split(df, method):
    """Roney: the non-estimable area, % of the total, split into fit_area_% (k_used < 3, lost to the
    planar fit) and conic_area_% (k_used = 3 and estimate NaN). NaN for the per-face method."""
    if method != "roney":
        return {"fit_area_%": np.nan, "conic_area_%": np.nan}
    lost, w = ~df["estimable"], df["area"]
    return {"fit_area_%": 100 * w[lost & (df["k_used"] < 3)].sum() / w.sum(),
            "conic_area_%": 100 * w[lost & (df["k_used"] == 3)].sum() / w.sum()}


@functools.cache
def all_subsets(scenario):
    """The 42 per-face runs of a scenario concatenated, with the columns k (maps of the subset) and
    sub. Read-only: shared between the plots."""
    return pd.concat([load(scenario, "per_face", sub).assign(k=sub.count("-") + 1, sub=sub)
                      for sub in subs(scenario, "per_face")])


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
                             "non_estimable_area_%": non_estimable(df), **roney_split(df, method),
                             "paired_area_%": 100 * df["area"][both].sum() / df["area"].sum(), **quantiles(df[both])})
                paired.append(df[both])
            comparison += rows
            comparison.append({"scenario": sc, "sub": "pooled", "method": method,
                               **{col: np.mean([r[col] for r in rows])
                                  for col in ("non_estimable_area_%", "fit_area_%", "conic_area_%")},
                               "paired_area_%": np.mean([r["paired_area_%"] for r in rows]),
                               **quantiles(pd.concat(paired))})
    return pd.DataFrame(per_face_table), pd.DataFrame(comparison)


def figure(name, title, draw, together=False, after=None):
    """Four PNGs, output/figures/NAME_{area,angle,mean_cv,ratio}.png. draw(ax, cols, ylabel) draws
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


def boxes(ax, x, x0, x1, groups, col, color=QSTYLE["q50"][0], label=None):
    """One box per group, in x from x0 to x1: from q25 to q75 of the absolute value of col, line at q50 and whiskers
    from q10 to q90, area-weighted."""
    q10, q25, q50, q75, q90 = (groups.apply(lambda s: wquantile(s[col].abs(), s["area"], q)).to_numpy()
                              for q in (0.10, 0.25, 0.50, 0.75, 0.90))
    ax.vlines(x, q10, q25, color=color, lw=1)
    ax.vlines(x, q75, q90, color=color, lw=1)
    ax.bar(x0, q75 - q25, width=x1 - x0, bottom=q25, align="edge", facecolor=to_rgba(color, 0.25), edgecolor=color,
           lw=1, zorder=3, label=label)
    ax.hlines(q50, x0, x1, color=color, lw=2, zorder=4)


def boxes_note(fig):
    """The legend of the boxes, below the figure."""
    fig.supxlabel("box: q25–q75 · line: q50 · whiskers: q10–q90 (area-weighted)", fontsize="medium")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    per_face_table, comparison = tables()
    per_face_table.to_csv(FIG / "per_face.csv", index=False, float_format="%.6g")
    comparison.to_csv(FIG / "comparison.csv", index=False, float_format="%.6g")

    # plot 0a: SCENARIO_COND, all subsets, q50, q90, q99 of |error| as a function of cond
    c = all_subsets(SCENARIO_COND)
    c = c[c["cond"].notna()].sort_values("cond")
    w = c["area"]
    c = c.assign(bin=np.minimum((10 * (w.cumsum() - 0.5 * w) / w.sum()).astype(int), 9))   # 10 equal-area bins
    est = c[c["estimable"]]
    xc = c.groupby("bin").apply(lambda s: wquantile(s["cond"], s["area"], 0.50)).to_numpy()
    fig, axes = plt.subplots(2, 2, figsize=(10, 9), layout="constrained")
    axes = axes.flat
    for ax, e in zip(axes, ERRORS):
        for q in QUANTILES:
            color, marker = QSTYLE[f"q{q}"]
            y = est.groupby("bin").apply(lambda s: wquantile(s[e].abs(), s["area"], q / 100)).to_numpy()
            ax.plot(xc, y, color=color, marker=marker, ms=5, lw=1.5, label=f"q{q}", zorder=3)
        ax.grid(alpha=0.3)
        ax.set_xlabel("cond, median of the bin (log scale)")
        ax.set_ylabel(f"{YLABEL[e]} (log scale)")
        log_scale(ax, "xy")
        ax.set_title(f"{NAMES[e]}, Spearman's ρ = {spearmanr(est['cond'], est[e].abs()).statistic:.2f}")
    ax = axes[3]
    ax.plot(xc, c.groupby("bin").apply(non_estimable), color=QSTYLE["q50"][0], lw=1.5, marker="o", ms=5)
    ax.grid(alpha=0.3)
    ax.set_xlabel("cond, median of the bin (log scale)")
    log_scale(ax, "x")
    y_axis(ax, "area where the estimate is not an ellipse [%]")
    ax.set_title("estimate is not an ellipse")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center", ncol=3, frameon=False,
               title="quantiles of |error| in the bin, area-weighted")
    fig.suptitle(f"{SCENARIO_COND}, all map subsets: error as a function of cond, 10 equal-area bins")
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

    # plot 0d: one box of cond per r, all the subsets of each scenario of SWEEP
    rs = [r for _, r in SWEEP]
    sweep = pd.concat([all_subsets(sc).assign(r=r) for sc, r in SWEEP])
    sweep = sweep[sweep["cond"].notna()]
    fig, ax = plt.subplots(figsize=(7.5, 4.6), layout="constrained")
    boxes(ax, np.array(rs), np.array(rs) - 0.1, np.array(rs) + 0.1, sweep.groupby("r"), "cond")
    ax.grid(alpha=0.3)
    ax.set_xticks(rs, [f"{r:g}" for r in rs])
    ax.set_xlabel("r = CV_L / CV_T")
    y_axis(ax, "cond", log=True)
    boxes_note(fig)
    fig.suptitle("uniform fields differing only in r: distribution of cond by r, all map subsets")
    fig.savefig(FIG / "plot0_ratio_cond.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # plot 1: for each scenario of PLOT1_SCENARIOS, as a function of k, points staggered within each k
    for sc in PLOT1_SCENARIOS:
        g1 = per_face_table[per_face_table["scenario"] == sc].copy()
        i, n = g1.groupby("k").cumcount(), g1.groupby("k")["k"].transform("size")
        g1["x"] = g1["k"] + 0.6 * (i - (n - 1) / 2) / n

        def draw1(ax, cols, ylabel):
            points(ax, g1["x"], g1, cols)
            ax.set_xticks(sorted(g1["k"].unique()))
            ax.set_xlabel("maps used, k")
            y_axis(ax, ylabel, log=len(cols) > 1)

        figure(f"plot1_maps_{sc}", f"{sc}: one point per map subset", draw1, together=True)

    # plot 2: per-face, as a function of r, the 42 subsets of each r staggered within it, ordered by k
    g2 = per_face_table.merge(pd.DataFrame(SWEEP, columns=["scenario", "r"]), on="scenario")
    g2 = g2.sort_values(["r", "k", "sub"], kind="stable")
    i, n = g2.groupby("r").cumcount(), g2.groupby("r")["r"].transform("size")
    g2["x"] = g2["r"] + 0.25 * (i - (n - 1) / 2) / n

    def draw2(ax, cols, ylabel):
        points(ax, g2["x"], g2, cols)
        ax.set_xticks(rs, [f"{r:g}" for r in rs])
        ax.set_xlabel("r = CV_L / CV_T (within each r the 42 subsets, k increasing from left to right)")
        y_axis(ax, ylabel, log=len(cols) > 1)

    figure("plot2_ratio", "uniform fields differing only in r, per-face: one point per map subset", draw2,
           together=True)

    # plot 3: cond by k and by r, boxes grouped by k, one per r
    ks = sorted(sweep["k"].unique())
    colors = dict(zip(SCENARIOS, OKABE_ITO))              # the colors of the scenarios in plot 4
    fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
    rows = []
    for j, (sc, r) in enumerate(SWEEP):
        s3 = sweep[sweep["r"] == r]
        x = np.array(ks) + 0.27 * (j - 1)
        boxes(ax, x, x - 0.11, x + 0.11, s3.groupby("k"), "cond", color=colors[sc], label=f"r = {r:g}")
        for k, g in s3.groupby("k"):
            rows.append({"k": k, "r": r, "scenario": sc, "cond_q50": wquantile(g["cond"], g["area"], 0.50),
                         "cond_q90": wquantile(g["cond"], g["area"], 0.90), "faces": len(g)})
    pd.DataFrame(rows).sort_values(["k", "r"]).to_csv(FIG / "cond_k_ratio.csv", index=False, float_format="%.6g")
    rho_k = spearmanr(sweep["cond"], sweep["k"]).statistic
    rho_r = spearmanr(sweep["cond"], sweep["r"]).statistic
    ax.grid(alpha=0.3)
    ax.set_xticks(ks)
    ax.set_xlabel("maps used, k")
    y_axis(ax, "cond", log=True)
    fig.legend(*ax.get_legend_handles_labels(), loc="outside right center", frameon=False)
    boxes_note(fig)
    fig.suptitle(f"uniform fields differing only in r, all map subsets: cond by k and by r\n"
                 f"Spearman's ρ(cond, k) = {rho_k:.2f}, ρ(cond, r) = {rho_r:.2f} (faces with cond, unweighted)")
    fig.savefig(FIG / "plot3_cond_k_ratio.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # plot 4: per-face vs Roney, one point per triplet and scenario
    triplets = comparison[comparison["sub"] != "pooled"]
    g4 = triplets[triplets["method"] == "roney"].merge(triplets[triplets["method"] == "per_face"],
                                                       on=["scenario", "sub"],
                                                       suffixes=("_roney", "_per_face"))

    def draw4(ax, cols, ylabel):
        (col,) = cols.values()
        for sc, color, marker in zip(SCENARIOS, OKABE_ITO, MARKERS):
            s = g4[g4["scenario"] == sc]
            ax.scatter(s[f"{col}_roney"], s[f"{col}_per_face"], s=22, color=color, marker=marker,
                       label=SCENARIO_LABELS[sc], zorder=3)
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


if __name__ == "__main__":
    main()
