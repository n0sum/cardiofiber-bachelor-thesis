#!/usr/bin/env python3
"""
Analysis of the experiments, all on the painted mean speed (mu of
simulator/painted_mu.ply). For each scenario of SCENARIOS (the same list, in
the same order, as run_all.sh: GRID, then CORE) the 42 subsets of at least 3
maps with the per-face estimator (output/SCENARIO/per_face/sub_*/); for each
scenario of GRID also the 20 triplets with method 2 of Roney et al. 2019
(output/SCENARIO/roney/sub_i-j-k/). Each run has estimate.ply (format D) and
estimate.csv.
    GRID   painted at the LAT noise levels sigma = 0, 0.05, 0.1, 0.3 ms: the two
           main axes, number of maps k and noise sigma
    CORE   the core sample in anisotropy: mu of painted and constant
           r = CV_L / CV_T = 1.2, 2, 2.8, at sigma = CORE_NOISE = 0.05 ms; its
           single condition is k = CORE_K = 4; per-face only, without Roney

Usage, from the repository root:
    python analysis/analysis.py              # -> output/figures/

load() is the only place where the errors are computed, the same for the two
methods: estimate.ply against the ground truth data/reference/SCENARIO/truth.ply,
both in format D (mu, alpha, fiber f).
One error per estimable face (mu of the estimate finite), the log-Euclidean
distance between the estimated and the true velocity tensor (Arsigny et al.
2006). The tensor of a face is M = CV_L^2 f f^T + CV_T^2 t t^T, t the transverse
direction in the plane of the face; its matrix logarithm has the same
directions and the logarithms of the eigenvalues:
log M = 2 mu I + 2 alpha (f f^T - t t^T). Then
  error = ||log M_est - log M_true||_F
        = 2 sqrt 2 sqrt((mu_est - mu)^2 + alpha_est^2 + alpha^2 - 2 alpha_est alpha cos(2 theta))
with theta the angle between estimated and true fiber, cos(2 theta) = 2 c^2 - 1,
c = |f_est . f_true| (the fibers of format D are unit and tangent to the face).
Dimensionless; with exact fiber and anisotropy, a relative error eps on both
speeds gives error ~ 2 sqrt 2 eps (0.05 ~ 1.8% on the speeds). With alpha = 0 (isotropic tissue) the fiber does not enter the error; it
does not depend on the unit of the speeds (mu shifts by a constant); it is the
same for speeds or slownesses (log M^-1 = -log M).
AREA-WEIGHTED statistics: q50, q90 and q99 (the values below which 50, 90, 99%
of the area lies). cond = sigma_max / sigma_min of the system of a face
(per-face estimator only, from estimate.csv): only in per_face.csv, no figure
uses it.

per_face.csv, one row per scenario and subset of the per-face method: non-estimable
area, errors on the estimable faces of the subset and q50, q90 and q99
of cond (area-weighted, on the faces with computed cond).

comparison.csv, one row per scenario of GRID, triplet and method (per_face, roney): non-estimable
area of the method alone, paired area (estimable by both methods with that triplet)
and errors on the paired faces. For the Roney rows the non-estimable area is also
split in two: fit_area_% (k_used < 3, lost to the planar fit) and conic_area_%
(k_used = 3 and estimate NaN: the conic is not an ellipse or the system is nearly
singular); empty on the per-face rows. The rows with sub = pooled gather the
paired faces of the 20 triplets of the scenario; their areas are the means over
the triplets.

grid.csv, one row per (sigma, k) of the scenarios of GRID; core.csv, one row per
scenario of CORE at k = CORE_K (r in place of noise_ms). The cell of (scenario, k),
cell(): all the subsets with k maps (20 triplets, 15 quadruplets, 6 quintuplets,
1 sextuplet) concatenated, each face once per subset; subsets = their number;
error_mean = sum area E / sum area and error_q50/q90/q99 area-weighted, on the
estimable rows; non_estimable_area_% = area of the non-estimable rows / area of
all the rows (the mean over the subsets of the non-estimable area of each).

Figures:
  plot_grid.png        the cells of grid.csv as a 4 x 4 colored table: k = 3..6 on
                       the horizontal axis, sigma = 0, 0.05, 0.1, 0.3 ms on the
                       vertical axis as categories (equispaced, sigma = 0 at the
                       bottom); color error_mean, sequential, log scale; in each
                       cell error_mean, error_q99 and the non-estimable area
  plot_core.png        the rows of core.csv as a function of r: left error_mean
                       and error_q99, right the non-estimable area
  plot3_error_k_ratio.png
                       boxes of the log-Euclidean error E grouped by k (3 to 6), in
                       each group one box per scenario of CORE, one color per r, on
                       the estimable faces of the subsets with that k, log scale;
                       boxes from q25 to q75, line at q50, whiskers from q10 to q90,
                       area-weighted; in the title Spearman's rho of E with k and
                       with r over the same faces, unweighted
  plot4_roney_noise_{area,error}.png
                       per-face vs Roney, one point per triplet and scenario of
                       GRID (80 points, one color and marker per sigma). _area:
                       non-estimable area of each
                       method alone, one panel; _error: three panels q50, q90, q99
                       on the paired faces, each with its own scale. Roney on the
                       horizontal axis and per-face on the vertical axis, linear
                       from 0 with the same limits, y = x; on the axes only 0,
                       minimum, centroid (mean of all the points, cross) and maximum
                       of the cloud, with dotted projections and without grid
"""

import functools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, to_rgba
from matplotlib.ticker import FormatStrFormatter, LogLocator, NullFormatter
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]                # repository root
sys.path.insert(0, str(ROOT / "estimator"))
from roney import read_ply                                # the PLY reader

REF = ROOT / "data" / "reference"
OUT = ROOT / "output"
FIG = OUT / "figures"
GRID = [("painted", 0.0), ("painted_noise0.05", 0.05), ("painted_noise0.1", 0.1),
        ("painted_noise0.3", 0.3)]                        # (scenario, LAT noise sigma in ms), increasing sigma
CORE = [("painted_r1.2_noise0.05", 1.2), ("painted_r2_noise0.05", 2.0),
        ("painted_r2.8_noise0.05", 2.8)]                  # (scenario, constant r), increasing r
CORE_K = 4                                                # maps of the core sample
CORE_NOISE = 0.05                                         # sigma (ms) of the scenarios of CORE
SCENARIOS = [sc for sc, _ in GRID] + [sc for sc, _ in CORE]   # per_face.csv: the order of run_all.sh
ERRORS = ("error",)
QUANTILES = (50, 90, 99)
YLABEL = {"area": "non-estimable area [%]",
          "error": "log-Euclidean error ‖log M̂ − log M‖"}
METHODS = {"per_face": "per-face", "roney": "Roney"}
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
    """Fiber, mu and alpha per face, from the face element of a format D PLY."""
    f = np.stack([face[k] for k in ("fx", "fy", "fz")], axis=1)
    return f, face["mu"], face["alpha"]


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
    indexed by face_id with area, estimable, cond (per-face) or k_used (Roney) and the log-Euclidean
    error, NaN on the non-estimable faces. The only place with the formula of the error."""
    area, (f_true, mu_true, alpha_true) = truth(scenario)
    case = OUT / scenario / method / sub
    face = read_ply(case / "estimate.ply")["face"]
    f, mu, alpha = field(face)
    df = pd.read_csv(case / "estimate.csv").set_index("face_id")
    df["area"] = area
    df["estimable"] = np.isfinite(face["mu"])
    c = np.abs(np.sum(f * f_true, axis=1))                # |cos theta|, theta angle between the fibers
    # log-Euclidean distance ||log M_hat - log M||_F (Arsigny et al. 2006)
    #   = 2 sqrt 2 sqrt(dmu^2 + alpha_hat^2 + alpha^2 - 2 alpha_hat alpha cos 2theta)
    # np.maximum: by rounding the radicand can be -1e-17
    df["error"] = 2 * np.sqrt(2) * np.sqrt(np.maximum(0, (mu - mu_true) ** 2 + alpha ** 2 + alpha_true ** 2
                                                      - 2 * alpha * alpha_true * (2 * c ** 2 - 1)))
    return df


def non_estimable(df):
    """% of non-estimable area."""
    w = df["area"]
    return 100 * w[~df["estimable"]].sum() / w.sum()


def quantiles(df):
    """q50, q90 and q99 of the error on the estimable faces."""
    est = df[df["estimable"]]
    return {f"{e}_q{q}": wquantile(est[e], est["area"], q / 100) for e in ERRORS for q in QUANTILES}


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
    """per_face.csv (all the scenarios of SCENARIOS) and comparison.csv (only those of GRID, the only
    ones with Roney), as DataFrames."""
    per_face_table, comparison = [], []
    for sc in SCENARIOS:
        per_face = {sub: load(sc, "per_face", sub) for sub in subs(sc, "per_face")}
        for sub, df in per_face.items():
            per_face_table.append({"scenario": sc, "sub": sub, "k": sub.count("-") + 1,
                                   "non_estimable_area_%": non_estimable(df), **quantiles(df), **cond_quantiles(df)})
        if sc not in dict(GRID):
            continue
        roney = {sub: load(sc, "roney", sub) for sub in subs(sc, "roney")}

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


def figure(name, title, draw, after=None):
    """Two PNGs, FIG/NAME_{area,error}.png. draw(ax, cols, ylabel) draws
    in one panel the columns of cols, {label: column}: the area in one panel, the error in one panel
    per quantile. after(fig), if given, is called before saving."""
    for q, ylabel in YLABEL.items():
        if q == "area":
            panels = {"non-estimable area": {"non-estimable area": "non_estimable_area_%"}}
        else:
            panels = {f"q{p}": {f"q{p}": f"{q}_q{p}"} for p in QUANTILES}
        fig, axes = plt.subplots(1, len(panels), figsize=(max(4.6 * len(panels), 7.5), 4.6), squeeze=False,
                                 layout="constrained")
        for ax, (panel, cols) in zip(axes[0], panels.items()):
            ax.grid(alpha=0.3)
            draw(ax, cols, ylabel)
            ax.set_title(panel)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside lower center", ncol=min(len(labels), 6 if len(panels) > 1 else 3),
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


def boxes(ax, x, x0, x1, groups, col, color, label=None):
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


def cell(scenario, k):
    """A cell of grid.csv and core.csv: all the per-face subsets with k maps of the scenario
    concatenated, each face once per subset. Mean and quantiles of E area-weighted on the estimable
    rows; the non-estimable area over all the rows, i.e. the mean over the subsets (same total area)."""
    d = all_subsets(scenario)
    d = d[d["k"] == k]
    est = d[d["estimable"]]
    return {"k": k, "subsets": d["sub"].nunique(), "non_estimable_area_%": non_estimable(d),
            "error_mean": (est["area"] * est["error"]).sum() / est["area"].sum(), **quantiles(d)}


def plot_grid(grid):
    """plot_grid.png: the cells of grid.csv as a colored table, k on x, sigma on y as categories."""
    ks = sorted(grid["k"].unique())
    noises = [n for _, n in GRID]                         # increasing sigma: row 0 (bottom) is sigma = 0
    cells = grid.set_index(["noise_ms", "k"])
    mean = np.array([[cells.loc[(n, k), "error_mean"] for k in ks] for n in noises])
    cmap = plt.get_cmap("viridis")
    norm = LogNorm(vmin=mean.min(), vmax=mean.max())
    fig, ax = plt.subplots(figsize=(8.2, 6.2), layout="constrained")
    im = ax.imshow(mean, origin="lower", cmap=cmap, norm=norm, aspect="auto")
    for i, n in enumerate(noises):
        for j, k in enumerate(ks):
            c = cells.loc[(n, k)]
            r, g, b, _ = cmap(norm(c["error_mean"]))
            luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
            ax.text(j, i, f"mean {c['error_mean']:.3g}\nq99 {c['error_q99']:.3g}\n"
                          f"n.e. {c['non_estimable_area_%']:.2g}%",
                    ha="center", va="center", fontsize=10, linespacing=1.5,
                    color="black" if luminance > 0.45 else "white")
    ax.set_xticks(range(len(ks)), [f"{k}" for k in ks])
    ax.set_yticks(range(len(noises)), [f"{n:g}" for n in noises])
    ax.set_xticks(np.arange(len(ks) + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(noises) + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="white", lw=2)
    ax.tick_params(which="minor", length=0)
    ax.set_xlabel("maps used, k (subsets: " + ", ".join(
        f"{int(cells.loc[(noises[0], k), 'subsets'])}" for k in ks) + ")")
    ax.set_ylabel("LAT noise σ [ms]")
    bar = fig.colorbar(im, ax=ax, shrink=0.9)
    bar.set_label("mean log-Euclidean error (area-weighted)")
    ticks = [t for t in (0.05, 0.1, 0.2, 0.5, 1, 2, 5) if norm.vmin <= t <= norm.vmax]
    bar.set_ticks([norm.vmin, *ticks, norm.vmax], labels=[f"{t:.2g}" for t in (norm.vmin, *ticks, norm.vmax)])
    bar.ax.yaxis.set_minor_formatter(NullFormatter())
    fig.suptitle("painted, per-face: error by number of maps and LAT noise")
    fig.supxlabel("mean and q99 on the estimable faces of all the subsets with k maps, area-weighted; "
                  "n.e.: non-estimable area", fontsize="medium")
    fig.savefig(FIG / "plot_grid.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_core(core):
    """plot_core.png: the rows of core.csv as a function of r, error left, non-estimable area right."""
    rs = core["r"].to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained")
    ax = axes[0]
    ax.plot(rs, core["error_mean"], color="#0072B2", marker="o", ms=6, lw=1.5, label="mean")
    ax.plot(rs, core["error_q99"], color="#D55E00", marker="s", ms=6, lw=1.5, ls="--", label="q99")
    ax.legend(frameon=False)
    y_axis(ax, YLABEL["error"])
    ax.set_title("error on the estimable faces, area-weighted")
    ax = axes[1]
    ax.plot(rs, core["non_estimable_area_%"], color="#000000", marker="o", ms=6, lw=1.5)
    y_axis(ax, YLABEL["area"])
    ax.set_title("non-estimable area")
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.set_xticks(rs, [f"{r:g}" for r in rs])
        ax.set_xlabel("r = CV_L / CV_T")
    n = int(core["subsets"].iloc[0])
    fig.suptitle(f"painted mean speed, constant r, k = {CORE_K} ({n} subsets), LAT noise σ = {CORE_NOISE:g} ms")
    fig.savefig(FIG / "plot_core.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot3():
    """plot3_error_k_ratio.png: boxes of E grouped by k, one per scenario of CORE."""
    s = pd.concat([all_subsets(sc).assign(r=r) for sc, r in CORE])
    s = s[s["error"].notna()]
    ks = sorted(s["k"].unique())
    fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
    for j, ((sc, r), color) in enumerate(zip(CORE, OKABE_ITO)):
        x = np.array(ks) + 0.27 * (j - (len(CORE) - 1) / 2)
        boxes(ax, x, x - 0.11, x + 0.11, s[s["r"] == r].groupby("k"), "error", color=color, label=f"r = {r:g}")
    rho_k = spearmanr(s["error"], s["k"]).statistic
    rho_r = spearmanr(s["error"], s["r"]).statistic
    ax.grid(alpha=0.3)
    ax.set_xticks(ks)
    ax.set_xlabel("maps used, k")
    y_axis(ax, YLABEL["error"], log=True)
    fig.legend(*ax.get_legend_handles_labels(), loc="outside right center", frameon=False)
    boxes_note(fig)
    fig.suptitle(f"painted mean speed, constant r, LAT noise σ = {CORE_NOISE:g} ms, all map subsets: "
                 f"log-Euclidean error E by k and by r\n"
                 f"Spearman's ρ(E, k) = {rho_k:.2f}, ρ(E, r) = {rho_r:.2f} (estimable faces, unweighted)")
    fig.savefig(FIG / "plot3_error_k_ratio.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot4(comparison, name, labels):
    """FIG/NAME_{area,error}.png: per-face vs Roney, one point per triplet and scenario of labels,
    {scenario: legend label}, one color and marker per scenario in the order of labels."""
    triplets = comparison[(comparison["sub"] != "pooled") & comparison["scenario"].isin(labels)]
    g4 = triplets[triplets["method"] == "roney"].merge(triplets[triplets["method"] == "per_face"],
                                                       on=["scenario", "sub"],
                                                       suffixes=("_roney", "_per_face")).reset_index(drop=True)

    def draw4(ax, cols, ylabel):
        (col,) = cols.values()
        for (sc, label), color, marker in zip(labels.items(), OKABE_ITO, MARKERS):
            s = g4[g4["scenario"] == sc]
            ax.scatter(s[f"{col}_roney"], s[f"{col}_per_face"], s=22, color=color, marker=marker,
                       label=label, zorder=3)
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
            setter(vals, ["0"] + [f"{t:.0f}" if t >= 100 else f"{t:#.3g}".rstrip(".")   # 19.0, not 19;
                                  for t in vals[1:]])                              # 1172, not 1.17e+03
        ax.grid(False)
        ax.set_aspect("equal")
        ax.set_xlabel(f"{METHODS['roney']}:\n{ylabel}")                 # two lines: the labels are long
        ax.set_ylabel(f"{METHODS['per_face']}:\n{ylabel}")

    figure(name, f"map triplets: per-face vs Roney, paired faces (area: each method alone), {len(g4)} points\n"
           "on the axes: min, centroid (mean of all the points) and max", draw4, after=separate_tick_labels)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    per_face_table, comparison = tables()
    per_face_table.to_csv(FIG / "per_face.csv", index=False, float_format="%.6g")
    comparison.to_csv(FIG / "comparison.csv", index=False, float_format="%.6g")

    ks = sorted(per_face_table["k"].unique())
    grid = pd.DataFrame([{"scenario": sc, "noise_ms": n, **cell(sc, k)} for sc, n in GRID for k in ks])
    core = pd.DataFrame([{"scenario": sc, "r": r, **cell(sc, CORE_K)} for sc, r in CORE])
    cols = ["k", "subsets", "non_estimable_area_%", "error_mean", "error_q50", "error_q90", "error_q99"]
    grid[["scenario", "noise_ms", *cols]].to_csv(FIG / "grid.csv", index=False, float_format="%.6g")
    core[["scenario", "r", *cols]].to_csv(FIG / "core.csv", index=False, float_format="%.6g")

    plot_grid(grid)
    plot_core(core)
    plot3()
    plot4(comparison, "plot4_roney_noise", {sc: f"σ = {n:g} ms" for sc, n in GRID})


if __name__ == "__main__":
    main()
