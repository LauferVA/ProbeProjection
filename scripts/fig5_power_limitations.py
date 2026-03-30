#!/usr/bin/env python3
"""
Figure 5 (Supplementary): Statistical Power and Technology Detection Limitations

Panel A: Minimum detectable absolute disparity (%) as a function of haplotype
         count in the smaller comparison group, for three Bonferroni-correction
         tiers (FISH, NGS, CMA).

Panel B: Technology detection capability matrix showing which probe-failure
         modes each clinical technology can identify through coordinate
         projection alone.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,       # TrueType in PDF for editability
    "ps.fonttype": 42,
})

# Colorblind-safe palette (Wong 2011, Nature Methods)
CB_BLUE = "#0072B2"
CB_ORANGE = "#E69F00"
CB_RED = "#D55E00"
CB_GREEN = "#009E73"
CB_PURPLE = "#CC79A7"
CB_GREY = "#999999"

PANEL_LABEL_KW = dict(fontsize=13, fontweight="bold", va="top", ha="left")


# ---------------------------------------------------------------------------
# Panel A helpers
# ---------------------------------------------------------------------------
def min_detectable_delta(
    n: np.ndarray,
    p2: float,
    alpha: float,
    power: float = 0.80,
    tol: float = 1e-6,
    max_iter: int = 200,
) -> np.ndarray:
    """
    Compute the minimum detectable absolute disparity delta = |p1 - p2| for a
    two-sided two-proportion z-test via bisection.

    Parameters
    ----------
    n : array-like
        Sample sizes for the smaller group (equal allocation assumed).
    p2 : float
        Baseline proportion (reference population coverage).
    alpha : float
        Significance level (already Bonferroni-adjusted).
    power : float
        Desired statistical power.
    tol : float
        Convergence tolerance for bisection.
    max_iter : int
        Maximum bisection iterations.

    Returns
    -------
    deltas : np.ndarray
        Minimum detectable absolute disparity for each n.
    """
    z_alpha = norm.ppf(1 - alpha / 2)  # two-sided
    n = np.asarray(n, dtype=float)
    deltas = np.empty_like(n)

    for i, ni in enumerate(n):
        lo, hi = 0.0, p2  # delta cannot exceed p2 (would make p1 < 0)
        for _ in range(max_iter):
            mid = (lo + hi) / 2
            p1 = p2 - mid
            p_bar = (p1 + p2) / 2
            se_null = np.sqrt(2 * p_bar * (1 - p_bar) / ni)
            se_alt = np.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / ni)
            if se_null == 0 or se_alt == 0:
                lo = mid
                continue
            achieved_power = norm.cdf((mid / se_alt) - (z_alpha * se_null / se_alt))
            if achieved_power < power:
                lo = mid
            else:
                hi = mid
            if hi - lo < tol:
                break
        deltas[i] = (lo + hi) / 2
    return deltas * 100  # convert to percentage


def draw_panel_a(ax: plt.Axes) -> None:
    """Draw the minimum detectable effect size panel."""
    ns = np.arange(5, 121)
    p2 = 0.95  # reference coverage

    corrections = [
        ("FISH (4 probes)", 0.05 / 4, CB_BLUE, "-"),
        ("NGS (~500 probes)", 0.05 / 500, CB_ORANGE, "--"),
        ("CMA (~6.85M probes)", 0.05 / 6_850_000, CB_RED, "-."),
    ]

    for label, alpha_adj, color, ls in corrections:
        deltas = min_detectable_delta(ns, p2, alpha_adj)
        ax.plot(ns, deltas, color=color, ls=ls, lw=1.6, label=label)

    # Actual sample sizes
    samples = [
        ("AFR", 106),
        ("NFE", 70),
        ("EAS", 48),
        ("SAS", 38),
        ("AMR", 10),
    ]
    for pop_label, n_val in samples:
        ax.axvline(n_val, color=CB_GREY, ls=":", lw=0.8, zorder=0)
        ax.text(
            n_val,
            1.02,
            f"  {pop_label}\n  n={n_val}",
            transform=ax.get_xaxis_transform(),
            fontsize=7,
            va="bottom",
            ha="center",
            color="#444444",
        )

    ax.set_xlabel("Number of haplotypes in smaller group")
    ax.set_ylabel("Minimum detectable absolute disparity (%)\n(80% power, two-sided)")
    ax.set_xlim(5, 120)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", frameon=True, edgecolor="#cccccc", fancybox=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(-0.08, 1.08, "A", transform=ax.transAxes, **PANEL_LABEL_KW)


# ---------------------------------------------------------------------------
# Panel B helpers
# ---------------------------------------------------------------------------
TECHNOLOGIES = ["FISH", "CMA", "MLPA", "NGS", "OGM"]
FAILURE_MODES = ["ABSENT\n(large SV)", "DIVERGENT\n(seq. mismatch)",
                 "DUPLICATED\n(multi-map.)", "REARRANGED\n(translocation)"]

# 1 = full (detectable), 0.5 = half (partially), 0 = empty (not detectable)
MATRIX = np.array([
    [1.0, 0.0, 0.5, 1.0],   # FISH
    [0.5, 0.0, 0.0, 0.0],   # CMA
    [0.5, 0.0, 0.0, 0.0],   # MLPA
    [1.0, 0.0, 0.5, 0.5],   # NGS
    [1.0, 0.0, 0.5, 1.0],   # OGM
])

# Colors
FILL_COLOR = CB_BLUE
CIRCLE_EDGE = "#333333"


def draw_panel_b(ax: plt.Axes) -> None:
    """Draw the technology detection capability matrix using scatter markers.

    Uses display-coordinate marker sizes so circles remain round regardless
    of axis aspect ratio. Column positions are simple integers (0..3) and
    the subplot stretches to fill the available width, giving labels room.
    """
    from matplotlib.patches import Wedge

    n_rows = len(TECHNOLOGIES)
    n_cols = len(FAILURE_MODES)
    marker_sz = 420  # scatter marker area (points^2)

    # Light grid background
    for r in range(n_rows):
        bg_color = "#f7f7f7" if r % 2 == 0 else "white"
        ax.add_patch(plt.Rectangle((-0.5, r - 0.5), n_cols, 1,
                                   fc=bg_color, ec="none", zorder=0))

    # Plot filled, empty, half-filled glyphs
    for r in range(n_rows):
        for c in range(n_cols):
            val = MATRIX[r, c]
            if val == 1.0:
                ax.scatter(c, r, s=marker_sz, c=FILL_COLOR,
                           edgecolors=CIRCLE_EDGE, linewidths=0.8, zorder=3)
            elif val == 0.0:
                ax.scatter(c, r, s=marker_sz, c="white",
                           edgecolors=CIRCLE_EDGE, linewidths=0.8, zorder=3)
            else:
                # Half-filled: draw two wedges in data coords with a fixed
                # radius that we compute from the axes transform so they
                # look circular. Fall back to overlay approach: white circle
                # underneath + left-half wedge on top.
                ax.scatter(c, r, s=marker_sz, c="white",
                           edgecolors=CIRCLE_EDGE, linewidths=0.8, zorder=3)
                # Overlay a left-half-circle marker via custom path
                from matplotlib.path import Path as MPath
                theta = np.linspace(np.pi / 2, 3 * np.pi / 2, 60)
                lh_verts = np.column_stack([np.cos(theta), np.sin(theta)])
                lh_verts = np.vstack([lh_verts, lh_verts[0]])  # close
                lh_codes = [MPath.MOVETO] + [MPath.LINETO] * (len(theta) - 1) + [MPath.CLOSEPOLY]
                lh_path = MPath(lh_verts, lh_codes)
                ax.scatter(c, r, s=marker_sz, marker=lh_path, c=FILL_COLOR,
                           edgecolors=CIRCLE_EDGE, linewidths=0.8, zorder=4)

    # Axes formatting
    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(-0.6, n_rows - 0.4)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(FAILURE_MODES, fontsize=8, ha="center")
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(TECHNOLOGIES, fontsize=9, fontweight="bold")
    ax.invert_yaxis()
    # Do NOT set aspect="equal" -- let columns spread to fill the width

    # Remove spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="both", length=0)

    # Legend
    legend_elements = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=FILL_COLOR,
                   markeredgecolor=CIRCLE_EDGE, markersize=10, lw=0,
                   label="Detectable"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=CB_ORANGE,
                   markeredgecolor=CIRCLE_EDGE, markersize=10, lw=0,
                   label="Partially detectable"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="white",
                   markeredgecolor=CIRCLE_EDGE, markersize=10, lw=0,
                   label="Not detectable"),
    ]
    ax.legend(
        handles=legend_elements,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.18),
        ncol=3,
        frameon=True,
        edgecolor="#cccccc",
        fancybox=False,
        fontsize=8,
    )

    ax.text(-0.12, 1.15, "B", transform=ax.transAxes, **PANEL_LABEL_KW)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    fig = plt.figure(figsize=(7, 6.5))

    # Panel A occupies upper ~55% of figure, Panel B the lower ~45%
    gs = fig.add_gridspec(2, 1, height_ratios=[1.1, 1], hspace=0.45)

    ax_a = fig.add_subplot(gs[0])
    draw_panel_a(ax_a)

    ax_b = fig.add_subplot(gs[1])
    draw_panel_b(ax_b)

    # Output paths
    out_dir = Path(__file__).resolve().parent.parent / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    for suffix in (".png", ".pdf"):
        out_path = out_dir / f"fig5_power_limitations{suffix}"
        fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"Saved: {out_path}")

    plt.close(fig)


if __name__ == "__main__":
    main()
