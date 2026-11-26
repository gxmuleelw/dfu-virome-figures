"""Draw Figure 1b–f from catalogue summaries, detection frequencies and cohort counts."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from main_1 import bars


BLUE = "#0072B2"
ORANGE = "#E69F00"
INK = "#20343D"
GREY = "#637A85"

STYLE = {
    "font.family": "Arial",
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "xtick.labelsize": 6.7,
    "ytick.labelsize": 6.7,
    "legend.fontsize": 6.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.55,
    "xtick.major.width": 0.55,
    "ytick.major.width": 0.55,
    "xtick.major.size": 2,
    "ytick.major.size": 2,
    "pdf.fonttype": 42,
    "axes.unicode_minus": False,
}


def heading(fig, letter, title, x, y):
    fig.text(x, y, letter, fontsize=10, weight="bold", color=INK)
    fig.text(x + 0.032, y, title, fontsize=8, weight="bold", color=INK)


def draw_panels(cfg):
    frequency = cfg["frequency"]
    cohorts = cfg["cohorts"]
    n = len(frequency)
    if (len(cfg["quality"]), len(cfg["clusters"]), len(cfg["ranks"])) != (2, 2, 4):
        raise ValueError("Expected two quality counts, two cluster counts and four rank counts")
    if not n or len(cohorts) != 3 or any(len(pair) != 2 for pair in cohorts):
        raise ValueError("Expected detection frequencies and three pairs of cohort counts")

    # Dimensions in inches of the quantitative portion of the manuscript figure.
    fig = plt.figure(figsize=(183 / 25.4, 220 / 25.4 - 386 / 72), facecolor="white")

    # b–d. Counts of viral operational taxonomic units by catalogue category.
    catalogue_panels = [
        ("b", "CheckV quality", 0.012, (0.075, 0.64, 0.225, 0.25),
         ["Complete", "High-quality"], cfg["quality"], ["#A9C3CE", GREY]),
        ("c", "Viral clusters", 0.342, (0.405, 0.64, 0.225, 0.25),
         ["VC members", "Residual"], cfg["clusters"], ["#8CACB8", "#CBD4D9"]),
        ("d", "Taxonomic resolution", 0.672, (0.735, 0.64, 0.245, 0.25),
         ["Class", "Order", "Family", "Uncl."], cfg["ranks"], ["#9CB4BF"] * 3 + ["#D1D8DD"]),
    ]
    for letter, title, title_x, bounds, categories, counts, colors in catalogue_panels:
        axis = fig.add_axes(bounds)
        heading(fig, letter, title, title_x, 0.956)
        bars(axis, categories, counts, colors, (0, 1000), [0, 400, 800], "vOTUs", 30)
        if letter == "c":
            axis.text(0.98, 0.95, cfg["cluster_label"], transform=axis.transAxes,
                      ha="right", va="top", fontsize=6.7)

    # e. Entry k contains the number of vOTUs detected in k + 1 participants.
    heading(fig, "e", "Detection frequency", 0.012, 0.49)
    axis = fig.add_axes((0.075, 0.17, 0.41, 0.24))
    axis.bar(range(1, n + 1), frequency, width=0.82, color=GREY, linewidth=0)
    axis.set_xlim(0, n + 1)
    axis.set_ylim(0, 420)
    axis.set_xticks(cfg["detection_xticks"])
    axis.set_yticks([0, 200, 400])
    axis.set_xlabel(f"Participants with detection (of {n})", labelpad=2)
    axis.set_ylabel("vOTUs", labelpad=2)
    axis.axvline(n / 2, color="#A3B6BF", ls="--", lw=0.6)
    axis.grid(axis="y", color="#E5ECEF", lw=0.45)
    axis.set_axisbelow(True)
    axis.text(0.98, 0.95, cfg["frequency_note"], transform=axis.transAxes,
              ha="right", va="top", fontsize=7)

    # f. The two counts in each analysis population follow the Wagner group order.
    heading(fig, "f", "Analysis populations", 0.55, 0.49)
    axis = fig.add_axes((0.68, 0.17, 0.29, 0.24))
    for index, (grades_1_3, grade_4) in enumerate(cohorts):
        total = grades_1_3 + grade_4
        axis.barh(index, grades_1_3, color=BLUE, height=0.55)
        axis.barh(index, grade_4, left=grades_1_3, color=ORANGE, height=0.55)
        axis.text(grades_1_3 / 2, index, str(grades_1_3), ha="center", va="center",
                  color="white", fontsize=7.3)
        axis.text(grades_1_3 + grade_4 / 2, index, str(grade_4), ha="center", va="center",
                  color=INK, fontsize=7.3)
        axis.text(total + 0.7, index, str(total), va="center", fontsize=7)
    axis.set_yticks([0, 1, 2], ["Detection", "Composition", "ANCOM-BC2"])
    axis.invert_yaxis()
    axis.set_xlim(0, 49)
    axis.set_xticks([0, 20, 40])
    axis.set_xlabel("Participants", labelpad=2)
    handles = [
        Line2D([], [], marker="o", color=BLUE, ls="", label="Wagner 1–3"),
        Line2D([], [], marker="s", color=ORANGE, ls="", label="Wagner 4"),
    ]
    fig.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.985, 0.005),
               frameon=False, ncol=2, handletextpad=0.25, columnspacing=0.8)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, metavar="INPUT_JSON")
    parser.add_argument("output", type=Path, metavar="OUTPUT_DIR")
    args = parser.parse_args()
    cfg = json.loads(args.inputs.read_text(encoding="utf-8-sig"))
    args.output.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(STYLE):
        figure = draw_panels(cfg)
        figure.savefig(args.output / "Figure1b-f.pdf")
        plt.close(figure)


if __name__ == "__main__":
    main()
