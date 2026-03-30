#!/usr/bin/env python3
"""
Figure 2: Technology Landscape -- Catalog Scope, Projection Summary, and
Population Disparity Profile.

Panel A uses real catalog data (profile and probe counts).
Panels B and C use simulated data labelled as such.

Output
------
figures/fig2_technology_landscape.png  (300 DPI)
figures/fig2_technology_landscape.pdf
"""

from __future__ import annotations

import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,  # TrueType in PDF (editable text)
        "ps.fonttype": 42,
    }
)

# Colorblind-safe palette (Wong 2011, Nature Methods)
WONG = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "amber": "#F0E442",
    "sky": "#56B4E9",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "black": "#000000",
}

# Projection-summary bar colours
COLOR_UNIVERSAL = WONG["green"]
COLOR_NONSPECIFIC = WONG["orange"]
COLOR_POPULATION = WONG["vermillion"]

# Technology colours for Panel C dot plot (4 technologies)
TECH_COLORS = {
    "FISH": WONG["blue"],
    "HLA": WONG["vermillion"],
    "PGx": WONG["purple"],
    "NGS": WONG["sky"],
}

# ---------------------------------------------------------------------------
# Panel A -- Catalog Scope (REAL DATA)
# ---------------------------------------------------------------------------
TECHNOLOGIES = [
    "CMA",
    "OGM",
    "Long-Read",
    "Southern Blot",
    "MLPA",
    "ddPCR",
    "CE/FA",
    "RT-PCR",
    "Methylation",
    "PGx",
    "HLA",
    "NGS",
    "FISH",
]

PROFILES = [1, 2, 4, 6, 6, 6, 9, 11, 16, 23, 38, 184, 319]
PROBES = [
    6_850_000,
    1_000_000,
    50,
    30,
    300,
    30,
    100,
    50,
    5_000,
    200,
    300,
    500_000,
    1_000,
]

# ---------------------------------------------------------------------------
# Panel B -- Projection Results Summary (SIMULATED)
# ---------------------------------------------------------------------------
PROJECTION_TECHS = [
    "CMA",
    "OGM",
    "MLPA",
    "NGS",
    "Methylation",
    "CE/FA",
    "RT-PCR",
    "ddPCR",
    "Southern Blot",
    "Long-Read",
    "PGx",
    "HLA",
    "FISH",
]

# (universal%, non-specific%, population-specific%)
PROJECTION_DATA: dict[str, tuple[float, float, float]] = {
    "CMA": (97.0, 2.9, 0.1),
    "OGM": (96.0, 3.0, 1.0),
    "MLPA": (98.0, 1.5, 0.5),
    "NGS": (95.0, 3.0, 2.0),
    "Methylation": (97.0, 2.0, 1.0),
    "CE/FA": (96.0, 3.0, 1.0),
    "RT-PCR": (97.0, 2.0, 1.0),
    "ddPCR": (96.0, 3.0, 1.0),
    "Southern Blot": (95.0, 3.5, 1.5),
    "Long-Read": (96.0, 3.0, 1.0),
    "PGx": (70.0, 15.0, 15.0),
    "HLA": (50.0, 20.0, 30.0),
    "FISH": (85.0, 10.0, 5.0),
}

# ---------------------------------------------------------------------------
# Panel C -- Population Disparity Profile (SIMULATED)
# ---------------------------------------------------------------------------
SUPERPOPS = ["AMR", "SAS", "NFE", "EAS", "AFR"]

DISPARITY_DATA: dict[str, dict[str, float]] = {
    "FISH": {"AFR": 0.96, "EAS": 0.97, "NFE": 0.92, "SAS": 0.96, "AMR": 0.95},
    "HLA": {"AFR": 0.82, "EAS": 0.85, "NFE": 0.88, "SAS": 0.84, "AMR": 0.80},
    "PGx": {"AFR": 0.88, "EAS": 0.90, "NFE": 0.93, "SAS": 0.91, "AMR": 0.87},
}

# Simulated 95% CI half-widths; AMR gets wider error bars.
CI_HALF: dict[str, dict[str, float]] = {
    "FISH": {"AFR": 0.010, "EAS": 0.008, "NFE": 0.012, "SAS": 0.009, "AMR": 0.025},
    "HLA": {"AFR": 0.018, "EAS": 0.016, "NFE": 0.014, "SAS": 0.017, "AMR": 0.035},
    "PGx": {"AFR": 0.015, "EAS": 0.013, "NFE": 0.011, "SAS": 0.014, "AMR": 0.030},
}


# ===================================================================
# Drawing functions
# ===================================================================


def draw_panel_a(ax_profiles: plt.Axes, ax_probes: plt.Axes) -> None:
    """Panel A: dual bar chart -- profiles (linear) and probes (log)."""
    y = np.arange(len(TECHNOLOGIES))
    bar_h = 0.55

    # -- Left: profile counts (grows leftward via inverted x-axis) --
    ax_profiles.barh(
        y, PROFILES, height=bar_h, color=WONG["blue"], edgecolor="white", linewidth=0.3
    )
    ax_profiles.set_yticks(y)
    ax_profiles.set_yticklabels(TECHNOLOGIES, fontsize=7, fontweight="medium")
    ax_profiles.set_xlabel("Number of profiles")
    ax_profiles.set_xlim(0, max(PROFILES) * 1.25)
    ax_profiles.invert_xaxis()  # profiles grow leftward
    ax_profiles.yaxis.tick_right()
    ax_profiles.yaxis.set_label_position("right")

    # Place count labels inside bars (white on blue) for large bars,
    # or just outside (dark text) for small bars.
    for yi, val in zip(y, PROFILES):
        if val >= 30:
            # Inside the bar
            ax_profiles.text(
                val * 0.5,
                yi,
                str(val),
                va="center",
                ha="center",
                fontsize=6,
                color="white",
                fontweight="bold",
            )
        else:
            # Outside the bar tip (towards the gap). With inverted x-axis,
            # x=0 is on the right, so place label at a negative x offset.
            ax_profiles.text(
                -max(PROFILES) * 0.02,
                yi,
                str(val),
                va="center",
                ha="right",
                fontsize=6,
                color="#333333",
            )

    # -- Right: probe counts (log scale) --
    ax_probes.barh(
        y,
        PROBES,
        height=bar_h,
        color=WONG["orange"],
        edgecolor="white",
        linewidth=0.3,
    )
    ax_probes.set_xscale("log")
    ax_probes.set_xlabel("Total probes (log scale)")
    ax_probes.set_yticks(y)
    ax_probes.set_yticklabels([])  # labels shared via left axis
    ax_probes.set_xlim(10, max(PROBES) * 8)

    for yi, val in zip(y, PROBES):
        label = f"{val:,}"
        ax_probes.text(
            val * 1.4,
            yi,
            label,
            va="center",
            ha="left",
            fontsize=5.5,
            color="#444444",
        )

    # Synchronise y-limits
    for ax in (ax_profiles, ax_probes):
        ax.set_ylim(-0.6, len(TECHNOLOGIES) - 0.4)
        ax.tick_params(axis="y", length=0)

    ax_profiles.spines["left"].set_visible(False)
    ax_profiles.spines["top"].set_visible(False)
    ax_probes.spines["right"].set_visible(False)
    ax_probes.spines["top"].set_visible(False)


def draw_panel_b(ax: plt.Axes) -> None:
    """Panel B: stacked horizontal bar chart of projection results."""
    techs = PROJECTION_TECHS
    universal = [PROJECTION_DATA[t][0] for t in techs]
    nonspecific = [PROJECTION_DATA[t][1] for t in techs]
    population = [PROJECTION_DATA[t][2] for t in techs]

    y = np.arange(len(techs))
    bar_h = 0.6

    ax.barh(
        y,
        universal,
        height=bar_h,
        color=COLOR_UNIVERSAL,
        label="Universal coverage",
        edgecolor="white",
        linewidth=0.3,
    )
    ax.barh(
        y,
        nonspecific,
        left=universal,
        height=bar_h,
        color=COLOR_NONSPECIFIC,
        label="Non-specific failure",
        edgecolor="white",
        linewidth=0.3,
    )
    left2 = [u + n for u, n in zip(universal, nonspecific)]
    ax.barh(
        y,
        population,
        left=left2,
        height=bar_h,
        color=COLOR_POPULATION,
        label="Population-specific disparity",
        edgecolor="white",
        linewidth=0.3,
    )

    ax.set_yticks(y)
    ax.set_yticklabels(techs)
    ax.set_xlim(0, 105)
    ax.set_xlabel("Probes (%)")
    ax.xaxis.set_major_formatter(mticker.PercentFormatter(100, decimals=0))

    ax.legend(
        loc="lower right",
        frameon=True,
        framealpha=0.9,
        edgecolor="#cccccc",
        fontsize=6.5,
    )

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Watermark
    ax.text(
        0.5,
        -0.13,
        "Simulated -- based on expected findings",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=7,
        fontstyle="italic",
        color="#888888",
    )


def draw_panel_c(ax: plt.Axes) -> None:
    """Panel C: grouped dot plot of coverage by superpopulation."""
    tech_list = list(DISPARITY_DATA.keys())
    n_techs = len(tech_list)
    offsets = np.linspace(-0.2, 0.2, n_techs)

    for i, tech in enumerate(tech_list):
        means = [DISPARITY_DATA[tech][sp] for sp in SUPERPOPS]
        cis = [CI_HALF[tech][sp] for sp in SUPERPOPS]
        y_pos = np.arange(len(SUPERPOPS)) + offsets[i]

        # Marker style: AMR (index 0 in SUPERPOPS) gets hollow marker
        for j, (yp, m, ci, sp) in enumerate(zip(y_pos, means, cis, SUPERPOPS)):
            is_amr = sp == "AMR"
            ax.errorbar(
                m,
                yp,
                xerr=ci,
                fmt="o" if not is_amr else "o",
                markersize=5 if not is_amr else 5,
                markerfacecolor=TECH_COLORS[tech] if not is_amr else "white",
                markeredgecolor=TECH_COLORS[tech],
                markeredgewidth=1.2 if is_amr else 0.8,
                color=TECH_COLORS[tech],
                elinewidth=1.5 if is_amr else 0.8,
                capsize=2.5 if is_amr else 1.5,
                label=tech if j == 0 else None,
                zorder=3,
            )

    ax.set_yticks(np.arange(len(SUPERPOPS)))
    ax.set_yticklabels(SUPERPOPS)
    ax.set_xlim(0.74, 1.01)
    ax.set_xlabel("Mean probe coverage rate")
    ax.xaxis.set_major_formatter(mticker.PercentFormatter(1.0, decimals=0))

    ax.axvline(1.0, color="#cccccc", linewidth=0.5, linestyle="--", zorder=0)

    # Reference line at x=0.90 for a "clinical concern" threshold
    ax.axvline(
        0.90, color="#dddddd", linewidth=0.7, linestyle=":", zorder=0, label=None
    )
    ax.text(
        0.898,
        len(SUPERPOPS) - 0.5,
        "90% threshold",
        fontsize=5.5,
        color="#999999",
        ha="right",
        va="bottom",
        rotation=90,
    )

    handles, labels = ax.get_legend_handles_labels()
    # De-duplicate
    seen: dict[str, int] = {}
    unique_handles, unique_labels = [], []
    for h, l in zip(handles, labels):
        if l not in seen:
            seen[l] = 1
            unique_handles.append(h)
            unique_labels.append(l)

    legend = ax.legend(
        unique_handles,
        unique_labels,
        loc="lower left",
        frameon=True,
        framealpha=0.9,
        edgecolor="#cccccc",
        fontsize=6.5,
        title="Technology",
        title_fontsize=7,
    )

    # Add note about AMR marker meaning
    ax.annotate(
        "Hollow markers = AMR (wider CI due to smaller reference panel)",
        xy=(0.98, 0.97),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize=5.5,
        color="#888888",
        fontstyle="italic",
    )

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Watermark
    ax.text(
        0.5,
        -0.16,
        "Simulated -- based on expected findings",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=7,
        fontstyle="italic",
        color="#888888",
    )


# ===================================================================
# Main
# ===================================================================


def main() -> None:
    fig = plt.figure(figsize=(7.0, 11.0), constrained_layout=False)

    # Create a 3-row gridspec: Panel A gets two sub-columns
    gs = fig.add_gridspec(
        nrows=3,
        ncols=1,
        hspace=0.40,
        left=0.14,
        right=0.94,
        top=0.96,
        bottom=0.06,
        height_ratios=[1.0, 1.0, 0.85],
    )

    # Panel A: split into two columns (profiles | probes)
    gs_a = gs[0].subgridspec(1, 2, wspace=0.30, width_ratios=[1, 1.3])
    ax_a_profiles = fig.add_subplot(gs_a[0, 0])
    ax_a_probes = fig.add_subplot(gs_a[0, 1])
    draw_panel_a(ax_a_profiles, ax_a_probes)

    # Panel A title spanning both sub-axes
    fig.text(
        0.53,
        0.965,
        "A    Catalog scope: profiles and probes per technology",
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
    )

    # Panel B
    ax_b = fig.add_subplot(gs[1])
    draw_panel_b(ax_b)
    ax_b.set_title(
        "B    Projected probe coverage by technology",
        loc="left",
        fontweight="bold",
        pad=8,
    )

    # Panel C
    ax_c = fig.add_subplot(gs[2])
    draw_panel_c(ax_c)
    ax_c.set_title(
        "C    Population disparity in probe coverage",
        loc="left",
        fontweight="bold",
        pad=8,
    )

    # Save
    out_dir = pathlib.Path(__file__).resolve().parent.parent / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    fig.savefig(out_dir / "fig2_technology_landscape.png", dpi=300)
    fig.savefig(out_dir / "fig2_technology_landscape.pdf")
    print(f"Saved to {out_dir / 'fig2_technology_landscape.png'}")
    print(f"Saved to {out_dir / 'fig2_technology_landscape.pdf'}")
    plt.close(fig)


if __name__ == "__main__":
    main()
