#!/usr/bin/env python3
"""
Analysis of the experiment grid: for each scenario aniso_zone_N the 42
subsets of at least 3 maps with the per-face estimator
(output/SCENARIO/per_face/sub_*/errors.csv) and the 20 triplets with method 2 of
Roney et al. 2019 (output/SCENARIO/roney/sub_i-j-k/errors.csv).

Usage, from the repository root:
    python analysis/analysis.py
Writes everything to output/figures/.

Errors per estimable face, read from errors.csv, which computes them (exportErrorsAsCsv in
estimator/src/GroundTruth.cpp and estimator/roney.py); AREA-WEIGHTED statistics:
  angle  err_angle, degrees
  cv_l   err_cv_l = 100 (CV_L_est - CV_L_true) / CV_L_true, %
  ratio  err_ratio = 100 (r_est - r_true) / r_true, %, r = CV_L / CV_T
q50, q90 and q99 are computed on the absolute value; for cv_l and ratio there is
also the signed median. cond = sigma_max / sigma_min of the system of a
face; its quantiles, area-weighted, are over the faces with computed cond
(estimable and discarded for lambda2 <= 1e-8 lambda1). d = sqrt(4 A / (pi N)) is the mean
zone diameter in mm, A the total area.

per_face.csv, one row per scenario and subset of the per-face method: non-estimable
area, errors on the estimable faces of the subset and q50, q90 and q99
of cond.

comparison.csv, one row per scenario, triplet and method (per_face, roney): non-estimable
area of the method alone, paired area (estimable by both
methods with that triplet) and errors on the paired faces. The rows with sub =
pooled gather the paired faces of the 20 triplets of the scenario; their
areas are the means over the triplets.

Boxes of plots 0a and 0b: from q25 to q75, line at q50, whiskers from q10 to q90,
area-weighted.
plot0_cond_errors.png, aniso_zone_1, faces with computed cond of the 42
subsets pooled, in 10 equal-area cond bins, each at its
median cond. Three panels angle, cv_l, ratio, one box per bin on its
estimable faces, log scale on both axes, in the title the Spearman
correlation between cond and error, unweighted; in a 2 x 2 square with the
area panel: the % of area of the bin where the estimate is not an ellipse, i.e. the
faces with cond but with the smaller eigenvalue lambda2 <= 1e-8 lambda1 (the estimated
conic is a hyperbola, or an ellipse flat to numerical zero). The faces discarded because
the smallest singular value is below threshold have no cond and do not appear.
plot0_maps_cond.png, aniso_zone_1, one box of cond per k, on the faces
with computed cond of all the subsets with k maps, log scale.
plot0_cond_area.png, aniso_zone_1, one point per subset from per_face.csv:
cond_q90 (log) on the horizontal axis, non-estimable area on the vertical axis, color and
marker by k, in the title Spearman's rho over the 42 points.

The other figures, each in four PNGs: _area (non-estimable area, one panel) and
_angle, _cv_l, _ratio. In plots 1 and 2 the latter have a single
panel with q50, q90 and q99 (circle, square, triangle, Okabe-Ito colors),
points only without lines, with the vertical axis on a log scale and labels in
plain numbers; in plots 3 and 4 three panels q50, q90, q99, each with its
own scale. The other vertical axes, including the _area files, are linear from 0;
in plot 4 also the horizontal axis, with the same limits.
  plot1_maps        aniso_zone_1, as a function of k: one point per
                    subset and quantile, staggered within each k
  plot2_zones       all maps, as a function of d: one point per
                    scenario and quantile
  plot3_zones_maps  as a function of d, for each k the median over the
                    subsets and the band from minimum to maximum
  plot4_roney       one point per triplet and scenario, Roney on the horizontal axis and
                    per-face on the vertical axis, same limits on both axes,
                    one color and one marker per scenario; on the axes only 0,
                    minimum, centroid (mean of the 80 points, cross) and
                    maximum of the cloud, with dotted projections and
                    without grid
Median, minimum and maximum over subsets count each subset once.
In plots 2 and 3 the dashed line is the diameter of Roney's disc,
2 patch.radius_mm of experiments/roney.conf.
"""

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
OUT = ROOT / "output"
FIG = OUT / "figures"
SCENARIO_MAPS = "aniso_zone_1"                            # plot 1
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
MARKERS = "osD^v"                                         # scenarios of plot 4
COLORS_N = {1: "#000000", 8: "#0072B2", 70: "#D55E00", 600: "#009E73"}   # plot 4, by N, Okabe-Ito
RONEY_DISC = 2 * float(next(l.split("=")[1] for l in open(ROOT / "experiments" / "roney.conf")
                            if l.startswith("patch.radius_mm")))


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


def load(path):
    """errors.csv with the error columns, as written in it; relative ones in percent, NaN on the
    non-estimable faces."""
    df = pd.read_csv(path).set_index("face_id")
    df["angle"] = df["err_angle"]
    df["cv_l"] = df["err_cv_l"]
    df["ratio"] = df["err_ratio"]
    return df


def non_estimable(df):
    """% of non-estimable area."""
    w = df["area"]
    return 100 * w[df["estimable"] == 0].sum() / w.sum()


def quantiles(df):
    """q50, q90 and q99 of the absolute value of each error on the estimable faces, then the signed
    median of cv_l and ratio."""
    est = df[df["estimable"] == 1]
    w = est["area"]
    out = {f"{e}_q{q}": wquantile(est[e].abs(), w, q / 100) for e in ERRORS for q in QUANTILES}
    for e in ("cv_l", "ratio"):
        out[f"{e}_median"] = wquantile(est[e], w, 0.50)
    return out


def cond_quantiles(df):
    """q50, q90 and q99 of cond on the faces with computed cond."""
    c = df[df["cond"].notna()]
    return {f"cond_q{q}": wquantile(c["cond"], c["area"], q / 100) for q in QUANTILES}


def tables():
    """per_face.csv and comparison.csv, as DataFrames."""
    per_face_table, comparison = [], []
    for sc_dir in sorted(OUT.glob("aniso_zone_*"), key=lambda p: int(p.name.removeprefix("aniso_zone_"))):
        sc = sc_dir.name
        n = int(sc.removeprefix("aniso_zone_"))
        per_face = {d.name: load(d / "errors.csv") for d in (sc_dir / "per_face").glob("sub_*")}
        roney = {d.name: load(d / "errors.csv") for d in (sc_dir / "roney").glob("sub_*")}
        area = next(iter(per_face.values()))["area"]
        head = {"scenario": sc, "N": n, "d": np.sqrt(4 * area.sum() / (np.pi * n))}

        for sub in sorted(per_face, key=lambda s: (s.count("-"), s)):
            df = per_face[sub]
            per_face_table.append({**head, "sub": sub, "k": sub.count("-") + 1,
                                   "non_estimable_area_%": non_estimable(df), **quantiles(df), **cond_quantiles(df)})

        for method, runs in (("per_face", per_face), ("roney", roney)):
            rows, paired = [], []
            for sub in sorted(roney):
                df = runs[sub]
                both = (roney[sub]["estimable"] == 1) & (per_face[sub]["estimable"] == 1)
                rows.append({**head, "sub": sub, "method": method, "non_estimable_area_%": non_estimable(df),
                             "paired_area_%": 100 * area[both].sum() / area.sum(), **quantiles(df[both])})
                paired.append(df[both])
            comparison += rows
            comparison.append({**head, "sub": "pooled", "method": method,
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


def d_axis(ax):
    """Horizontal axis d on a log scale, with the line of the diameter of Roney's disc."""
    ax.axvline(RONEY_DISC, color="#555", lw=1, ls="--", label=f"Roney's disc diameter ({RONEY_DISC:g} mm)")
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.set_xlabel("mean zone diameter d [mm]")


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

    # plot 0a: aniso_zone_1, all subsets, error as a function of cond
    df = pd.concat([load(d / "errors.csv").assign(k=d.name.count("-") + 1)
                    for d in (OUT / SCENARIO_MAPS / "per_face").glob("sub_*")])
    c = df[df["cond"].notna()].sort_values("cond")
    w = c["area"]
    c["bin"] = np.minimum((10 * (w.cumsum() - 0.5 * w) / w.sum()).astype(int), 9)   # 10 equal-area bins
    est = c[c["estimable"] == 1]
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
    fig.suptitle(f"{SCENARIO_MAPS}, all map subsets: error as a function of cond,\nequal-area bins")
    fig.savefig(FIG / "plot0_cond_errors.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # plot 0b: aniso_zone_1, one box of cond per k
    ks = np.array(sorted(c["k"].unique()))
    fig, ax = plt.subplots(figsize=(7.5, 4.6), layout="constrained")
    boxes(ax, ks, ks - 0.25, ks + 0.25, c.groupby("k"), "cond")
    ax.grid(alpha=0.3)
    ax.set_xticks(ks)
    ax.set_xlabel("maps used, k")
    y_axis(ax, "cond", log=True)
    boxes_note(fig)
    fig.suptitle(f"{SCENARIO_MAPS}: distribution of cond by number of maps")
    fig.savefig(FIG / "plot0_maps_cond.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # plot 0c: aniso_zone_1, non-estimable area vs cond_q90, one point per subset
    g0 = per_face_table[per_face_table["scenario"] == SCENARIO_MAPS]
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
    fig.suptitle(f"{SCENARIO_MAPS}: non-estimable area and tail of cond, one point per map subset "
                 f"(Spearman's ρ = {rho:.2f})")
    fig.savefig(FIG / "plot0_cond_area.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # plot 1: aniso_zone_1, as a function of k, points staggered within each k
    g1 = per_face_table[per_face_table["scenario"] == SCENARIO_MAPS].copy()
    i, n = g1.groupby("k").cumcount(), g1.groupby("k")["k"].transform("size")
    g1["x"] = g1["k"] + 0.6 * (i - (n - 1) / 2) / n


    def draw1(ax, cols, ylabel):
        points(ax, g1["x"], g1, cols)
        ax.set_xticks(sorted(g1["k"].unique()))
        ax.set_xlabel("maps used, k")
        y_axis(ax, ylabel, log=len(cols) > 1)

    figure("plot1_maps", f"{SCENARIO_MAPS}: one point per map subset", draw1, together=True)

    # plot 2: all maps, as a function of d
    g2 = per_face_table[per_face_table["k"] == per_face_table["k"].max()].sort_values("d")

    def draw2(ax, cols, ylabel):
        points(ax, g2["d"], g2, cols)
        d_axis(ax)
        y_axis(ax, ylabel, log=len(cols) > 1)

    figure("plot2_zones", f"aniso_zone_N, per-face, all {g2['k'].iloc[0]} maps: one point per scenario", draw2,
           together=True)

    # plot 3: as a function of d, one line and one band per k
    ks = sorted(per_face_table["k"].unique())
    colors = dict(zip(ks, plt.cm.viridis(np.linspace(0, 0.55, len(ks)))))

    def draw3(ax, cols, ylabel):
        (col,) = cols.values()
        for k in ks:
            s = per_face_table[per_face_table["k"] == k].groupby("d")[col]
            ax.fill_between(s.min().index, s.min(), s.max(), color=colors[k], alpha=0.15, lw=0)
            ax.plot(s.median().index, s.median(), color=colors[k], lw=2, marker="o", ms=4, label=f"k = {k}")
        d_axis(ax)
        y_axis(ax, ylabel)

    figure("plot3_zones_maps", "aniso_zone_N: median over map subsets of k maps, band from minimum to maximum",
           draw3)

    # plot 4: per-face vs Roney, one point per triplet and scenario
    triplets = comparison[comparison["sub"] != "pooled"]
    g4 = triplets[triplets["method"] == "roney"].merge(triplets[triplets["method"] == "per_face"],
                                                       on=["scenario", "N", "d", "sub"],
                                                       suffixes=("_roney", "_per_face"))
    scenarios = g4.drop_duplicates("scenario").sort_values("d")

    def draw4(ax, cols, ylabel):
        (col,) = cols.values()
        for (sc, n, d), marker in zip(scenarios[["scenario", "N", "d"]].itertuples(index=False), MARKERS):
            s = g4[g4["scenario"] == sc]
            ax.scatter(s[f"{col}_roney"], s[f"{col}_per_face"], s=22, color=COLORS_N[n], marker=marker,
                       label=f"N = {n} (d = {d:.0f} mm)" if d >= 10 else f"N = {n} (d = {d:.2g} mm)", zorder=3)
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
           "on the axes: min, centroid (mean over triplets) and max", draw4, after=separate_tick_labels)


if __name__ == "__main__":
    main()
