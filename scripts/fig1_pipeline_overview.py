#!/usr/bin/env python3
"""
fig1_pipeline_overview.py

Generate Figure 1: ProbeProjection pipeline schematic.
Produces publication-quality PNG (300 DPI) and PDF versions.

Dependencies: matplotlib, numpy (standard scientific Python stack).
"""

import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np


# ---------------------------------------------------------------------------
# Colour palette -- muted, publication-friendly
# ---------------------------------------------------------------------------
COL_INPUT_BOX = "#D6E4F0"       # light steel blue
COL_INPUT_BORDER = "#7A9CC6"    # medium slate blue
COL_PROC_BOX = "#E8E8E8"        # light warm grey
COL_PROC_BORDER = "#999999"     # mid grey
COL_OUTPUT_BOX = "#D5E8D4"      # soft sage green
COL_OUTPUT_BORDER = "#82B366"   # muted green
COL_SECTION_LABEL = "#3B3B3B"   # near-black
COL_BODY_TEXT = "#4A4A4A"       # dark grey
COL_ARROW = "#666666"           # medium grey
COL_TECH_TAG = "#F0E6D3"        # warm parchment
COL_TECH_TAG_BORDER = "#C4A882" # warm tan


# ---------------------------------------------------------------------------
# Layout constants (figure coordinates, 0-1 normalised then scaled)
# ---------------------------------------------------------------------------
FIG_W, FIG_H = 7.0, 4.2  # inches

# Column centres (x)
X_INPUT = 0.14
X_PROC = 0.50
X_OUTPUT = 0.86

# Box dimensions
BOX_W_INPUT = 0.22
BOX_W_PROC = 0.22
BOX_W_OUTPUT = 0.20
BOX_H = 0.145
BOX_RAD = 0.015  # corner radius (in axes fraction)

# Y positions for three rows (top to bottom)
Y_TOP = 0.78
Y_MID = 0.50
Y_BOT = 0.22


# ---------------------------------------------------------------------------
# Helper: draw a rounded-rectangle box with title + subtitle
# ---------------------------------------------------------------------------
def draw_box(ax, cx, cy, w, h, title, subtitle, facecolor, edgecolor,
             title_size=7.5, sub_size=6.0, bold_title=True):
    """Draw a FancyBboxPatch centred at (cx, cy) with text inside."""
    rect = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle=f"round,pad={BOX_RAD}",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=1.0,
        transform=ax.transAxes,
        zorder=2,
    )
    ax.add_patch(rect)

    weight = "bold" if bold_title else "normal"
    ax.text(cx, cy + 0.025, title,
            ha="center", va="center", fontsize=title_size,
            fontweight=weight, color=COL_SECTION_LABEL,
            transform=ax.transAxes, zorder=3)
    ax.text(cx, cy - 0.020, subtitle,
            ha="center", va="center", fontsize=sub_size,
            color=COL_BODY_TEXT, transform=ax.transAxes, zorder=3,
            style="italic", wrap=False)
    return rect


# ---------------------------------------------------------------------------
# Helper: draw an arrow between two boxes
# ---------------------------------------------------------------------------
def draw_arrow(ax, x0, y0, x1, y1, color=COL_ARROW):
    """Straight arrow from (x0,y0) to (x1,y1) in axes coordinates."""
    arrow = FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle="->,head_width=4,head_length=4",
        color=color,
        linewidth=1.2,
        transform=ax.transAxes,
        zorder=1,
        connectionstyle="arc3,rad=0.0",
    )
    ax.add_patch(arrow)


def draw_curved_arrow(ax, x0, y0, x1, y1, rad=0.15, color=COL_ARROW):
    """Curved arrow between two points."""
    arrow = FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle="->,head_width=4,head_length=4",
        color=color,
        linewidth=1.2,
        transform=ax.transAxes,
        zorder=1,
        connectionstyle=f"arc3,rad={rad}",
    )
    ax.add_patch(arrow)


# ---------------------------------------------------------------------------
# Helper: section header
# ---------------------------------------------------------------------------
def section_header(ax, cx, y, label, color):
    """Draw a section title above a column."""
    ax.text(cx, y, label,
            ha="center", va="center", fontsize=9, fontweight="bold",
            color=color, transform=ax.transAxes, zorder=3)


# ---------------------------------------------------------------------------
# Helper: technology-count tag
# ---------------------------------------------------------------------------
def tech_tag(ax, cx, cy, text, w=0.22):
    """Small rounded tag below the input column."""
    tag_h = 0.055
    rect = FancyBboxPatch(
        (cx - w / 2, cy - tag_h / 2), w, tag_h,
        boxstyle=f"round,pad=0.008",
        facecolor=COL_TECH_TAG,
        edgecolor=COL_TECH_TAG_BORDER,
        linewidth=0.8,
        transform=ax.transAxes,
        zorder=2,
    )
    ax.add_patch(rect)
    ax.text(cx, cy, text,
            ha="center", va="center", fontsize=5.2,
            color=COL_BODY_TEXT, transform=ax.transAxes, zorder=3)


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------
def make_figure():
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # ---- Section headers ----
    section_header(ax, X_INPUT, 0.95, "INPUT", COL_INPUT_BORDER)
    section_header(ax, X_PROC, 0.95, "PROCESSING", COL_PROC_BORDER)
    section_header(ax, X_OUTPUT, 0.95, "OUTPUT", COL_OUTPUT_BORDER)

    # Thin vertical separator lines
    for xsep in [0.32, 0.70]:
        ax.plot([xsep, xsep], [0.05, 0.92], color="#D0D0D0", linewidth=0.6,
                linestyle="--", transform=ax.transAxes, zorder=0)

    # ---- INPUT boxes ----
    draw_box(ax, X_INPUT, Y_TOP, BOX_W_INPUT, BOX_H,
             "Profile Catalog",
             "1,628 profiles, 16 technologies\n625 eligible (GRCh38)",
             COL_INPUT_BOX, COL_INPUT_BORDER, sub_size=5.5)

    draw_box(ax, X_INPUT, Y_MID, BOX_W_INPUT, BOX_H,
             "HPRC Pangenome",
             "24 per-chr .og graphs\n~272 haplotype paths",
             COL_INPUT_BOX, COL_INPUT_BORDER, sub_size=5.5)

    draw_box(ax, X_INPUT, Y_BOT, BOX_W_INPUT, BOX_H,
             "Ancestry Registry",
             "136 samples\n5 superpopulations",
             COL_INPUT_BOX, COL_INPUT_BORDER, sub_size=5.5)

    # Superpopulation labels below the ancestry box
    ax.text(X_INPUT, Y_BOT - 0.09,
            "AFR  |  EAS  |  NFE  |  SAS  |  AMR",
            ha="center", va="center", fontsize=5.0,
            color=COL_INPUT_BORDER, fontweight="bold",
            transform=ax.transAxes, zorder=3,
            fontfamily="monospace")

    # ---- Technology scope tag below input column ----
    tech_tag(ax, X_INPUT, 0.065,
             "FISH (319)  NGS (184)  HLA (38)  PGx (23)  +9 others",
             w=0.24)

    # ---- PROCESSING boxes ----
    draw_box(ax, X_PROC, Y_TOP, BOX_W_PROC, BOX_H,
             "Probe Extraction",
             "Per-profile coordinate\nextraction (chr, start, end)",
             COL_PROC_BOX, COL_PROC_BORDER, sub_size=5.5)

    draw_box(ax, X_PROC, Y_MID, BOX_W_PROC, BOX_H,
             "Coordinate Projection",
             "odgi position: project each\nprobe to each haplotype path",
             COL_PROC_BOX, COL_PROC_BORDER, sub_size=5.5)

    draw_box(ax, X_PROC, Y_BOT, BOX_W_PROC, BOX_H,
             "Population Stratification",
             "Chi-squared / Fisher exact\nBonferroni correction",
             COL_PROC_BOX, COL_PROC_BORDER, sub_size=5.5)

    # ---- OUTPUT boxes ----
    draw_box(ax, X_OUTPUT, Y_TOP, BOX_W_OUTPUT, BOX_H,
             "Per-Profile",
             "coverage.json\ndiversity.json",
             COL_OUTPUT_BOX, COL_OUTPUT_BORDER, sub_size=5.5)

    draw_box(ax, X_OUTPUT, Y_MID, BOX_W_OUTPUT, BOX_H,
             "Cross-Profile",
             "summary.json\nBED files",
             COL_OUTPUT_BOX, COL_OUTPUT_BORDER, sub_size=5.5)

    draw_box(ax, X_OUTPUT, Y_BOT, BOX_W_OUTPUT, BOX_H,
             "Distribution",
             "UCSC Track Hub",
             COL_OUTPUT_BOX, COL_OUTPUT_BORDER, sub_size=5.5)

    # ---- Arrows: INPUT -> PROCESSING ----
    x_in_right = X_INPUT + BOX_W_INPUT / 2 + 0.01
    x_proc_left = X_PROC - BOX_W_PROC / 2 - 0.01

    # Profile Catalog -> Probe Extraction
    draw_arrow(ax, x_in_right, Y_TOP, x_proc_left, Y_TOP)
    # HPRC Pangenome -> Coordinate Projection
    draw_arrow(ax, x_in_right, Y_MID, x_proc_left, Y_MID)
    # Ancestry Registry -> Population Stratification
    draw_arrow(ax, x_in_right, Y_BOT, x_proc_left, Y_BOT)

    # ---- Arrows: PROCESSING vertical chain ----
    draw_arrow(ax, X_PROC, Y_TOP - BOX_H / 2 - 0.01,
               X_PROC, Y_MID + BOX_H / 2 + 0.01)
    draw_arrow(ax, X_PROC, Y_MID - BOX_H / 2 - 0.01,
               X_PROC, Y_BOT + BOX_H / 2 + 0.01)

    # ---- Arrows: PROCESSING -> OUTPUT ----
    x_proc_right = X_PROC + BOX_W_PROC / 2 + 0.01
    x_out_left = X_OUTPUT - BOX_W_OUTPUT / 2 - 0.01

    # Probe Extraction -> Per-Profile (coverage comes from extraction too)
    draw_arrow(ax, x_proc_right, Y_TOP, x_out_left, Y_TOP)
    # Coordinate Projection -> Cross-Profile
    draw_arrow(ax, x_proc_right, Y_MID, x_out_left, Y_MID)
    # Population Stratification -> Distribution
    draw_arrow(ax, x_proc_right, Y_BOT, x_out_left, Y_BOT)

    # ---- Diagonal feed-forward arrows (inputs feed downstream steps) ----
    # Profile Catalog also feeds Coordinate Projection (probes projected)
    draw_curved_arrow(ax, x_in_right, Y_TOP - 0.03,
                      x_proc_left, Y_MID + 0.03, rad=0.12)
    # HPRC also feeds Population Stratification (haplotype context)
    draw_curved_arrow(ax, x_in_right, Y_MID - 0.03,
                      x_proc_left, Y_BOT + 0.03, rad=0.12)

    # ---- Subtle step numbers on processing boxes ----
    for i, y in enumerate([Y_TOP, Y_MID, Y_BOT], start=1):
        ax.text(X_PROC - BOX_W_PROC / 2 + 0.015, y + BOX_H / 2 - 0.018,
                str(i), ha="center", va="center", fontsize=7,
                fontweight="bold", color="#FFFFFF",
                transform=ax.transAxes, zorder=4,
                bbox=dict(boxstyle="circle,pad=0.15", facecolor=COL_PROC_BORDER,
                          edgecolor="none"))

    # ---- Figure title ----
    ax.text(0.50, 0.995,
            "ProbeProjection Pipeline Overview",
            ha="center", va="top", fontsize=10.5, fontweight="bold",
            color=COL_SECTION_LABEL, transform=ax.transAxes)

    fig.subplots_adjust(left=0.01, right=0.99, top=0.97, bottom=0.01)
    return fig


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
def main():
    out_dir = pathlib.Path(__file__).resolve().parent.parent / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    fig = make_figure()

    png_path = out_dir / "fig1_pipeline_overview.png"
    pdf_path = out_dir / "fig1_pipeline_overview.pdf"

    fig.savefig(str(png_path), dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    fig.savefig(str(pdf_path), bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)

    print(f"Saved PNG: {png_path}")
    print(f"Saved PDF: {pdf_path}")


if __name__ == "__main__":
    main()
