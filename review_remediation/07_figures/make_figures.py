#!/usr/bin/env python3
"""Regenerate the four adjudicated PLOS ONE main figures from audited tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path as MplPath

from figure_style import (
    BLUE, ORANGE, GREEN, MAGENTA, GREY, LIGHT, DARK, CONTRACT, STYLE,
    panel_label, panel_title, reference_line, register_box_text,
    register_flow_box, register_legend_data, register_orthogonal_arrow,
    register_phase_label, save_main_figure,
)


PROJECT = Path(__file__).resolve().parents[2]
ROOT = PROJECT / "review_remediation"
OUT = ROOT / "07_figures"
OUT.mkdir(parents=True, exist_ok=True)


def read_counts() -> dict:
    path = ROOT / "03_systematic_search/search_counts.json"
    counts = json.loads(path.read_text(encoding="utf-8"))
    required = {
        "records_identified_total", "records_identified_GEO", "records_identified_PubMed",
        "records_identified_BioProject", "records_identified_SRA",
        "records_identified_Europe_PMC", "records_identified_BioStudies",
        "duplicates_or_component_records_removed", "records_screened",
        "records_excluded_before_report_assessment", "reports_sought_for_retrieval",
        "reports_not_retrieved", "reports_assessed_for_eligibility", "reports_excluded",
        "reports_excluded_no_P17_contrast", "reports_excluded_insufficient_replication",
        "reports_excluded_wrong_tissue_or_modality",
        "reports_excluded_no_recoverable_expression_or_contrast",
        "reports_excluded_other", "eligible_GEO_cohorts", "eligible_non_GEO_cohorts",
        "unique_eligible_cohorts_included", "unique_GEO_series",
        "unique_GEO_series_included", "unique_GEO_series_excluded",
    }
    missing = required - set(counts)
    if missing:
        raise RuntimeError(f"search_counts.json missing fields: {sorted(missing)}")
    if counts["unique_GEO_series"] != counts["unique_GEO_series_included"] + counts["unique_GEO_series_excluded"]:
        raise RuntimeError("GEO series arithmetic is inconsistent")
    assertions = {
        "1029 - 774 = 255": counts["records_identified_total"] - counts["duplicates_or_component_records_removed"] == counts["records_screened"],
        "255 - 199 = 56": counts["records_screened"] - counts["records_excluded_before_report_assessment"] == counts["reports_sought_for_retrieval"],
        "56 - 0 - 40 = 16": counts["reports_sought_for_retrieval"] - counts["reports_not_retrieved"] - counts["reports_excluded"] == counts["unique_eligible_cohorts_included"],
        "14 + 2 = 16": counts["eligible_GEO_cohorts"] + counts["eligible_non_GEO_cohorts"] == counts["unique_eligible_cohorts_included"],
    }
    failed = [expression for expression, passed in assertions.items() if not passed]
    if failed:
        raise RuntimeError(f"PRISMA count arithmetic failed: {failed}")
    return counts


def add_flow_box(ax, fig, *, x: float, y: float, width: float, height: float,
                 lines: list[str], name: str, column: str, fill: str = "white"):
    patch = FancyBboxPatch(
        (x, y), width, height, boxstyle="round,pad=0,rounding_size=0.004",
        facecolor=fill, edgecolor="#333333", linewidth=0.9,
    )
    ax.add_patch(patch)
    text = ax.text(
        x + 0.018, y + height / 2, "\n".join(lines),
        ha="left", va="center", multialignment="left",
        fontsize=8.5, linespacing=1.18,
    )
    register_box_text(fig, patch, text, name)
    register_flow_box(fig, patch, name, column)
    return {"x": x, "y": y, "w": width, "h": height, "patch": patch, "text": text}


def orthogonal_arrow(ax, fig, source, target, name: str, vertices=None):
    if vertices is None:
        vertices = [
            (source["x"] + source["w"] / 2, source["y"]),
            (target["x"] + target["w"] / 2, target["y"] + target["h"]),
        ]
    path = MplPath(vertices, [MplPath.MOVETO] + [MplPath.LINETO] * (len(vertices) - 1))
    arrow_patch = FancyArrowPatch(
        path=path, arrowstyle="-|>", color="#555555", linewidth=0.9,
        mutation_scale=9, shrinkA=0, shrinkB=0, capstyle="butt", joinstyle="miter",
    )
    ax.add_patch(arrow_patch)
    register_orthogonal_arrow(fig, ax, vertices, source["patch"], target["patch"], name)


def add_phase_label(ax, fig, label: str, y: float):
    text = ax.text(
        .070, y, label, ha="center", va="center", fontsize=9.2,
        fontweight="bold", color="#777777", rotation=90,
    )
    register_phase_label(fig, text, label)


def build_fig1() -> plt.Figure:
    c = read_counts()
    fig, ax = plt.subplots(figsize=(CONTRACT["output"]["width_in"], CONTRACT["figures"]["Fig1"]["height_in"]))
    fig.subplots_adjust(left=.03, right=.985, top=.985, bottom=.03)
    ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    add_phase_label(ax, fig, "IDENTIFICATION", .865)
    add_phase_label(ax, fig, "SCREENING", .555)
    add_phase_label(ax, fig, "ELIGIBILITY", .345)
    add_phase_label(ax, fig, "INCLUDED", .108)

    identified = add_flow_box(ax, fig, x=.185, y=.735, width=.47, height=.235,
                              name="records_identified", column="main", lines=[
        "Records identified from databases and",
        "repositories",
        f"Total records (n = {c['records_identified_total']})",
        "",
        f"GEO objects (n = {c['records_identified_GEO']})",
        f"PubMed reports (n = {c['records_identified_PubMed']})",
        f"BioProject/SRA records (n = {c['records_identified_BioProject'] + c['records_identified_SRA']})",
        f"Europe PMC records (n = {c['records_identified_Europe_PMC']})",
        f"BioStudies records (n = {c['records_identified_BioStudies']})",
    ])
    removed = add_flow_box(ax, fig, x=.678, y=.795, width=.302, height=.145,
                           name="removed_before_screening", column="exclusion", lines=[
        "Records removed before",
        "screening",
        "",
        "Duplicate/component",
        f"records removed (n = {c['duplicates_or_component_records_removed']})",
    ])
    screened = add_flow_box(ax, fig, x=.185, y=.585, width=.47, height=.085,
                            name="screened", column="main", lines=[
        "Records screened",
        f"n = {c['records_screened']}",
    ])
    screen_excl = add_flow_box(ax, fig, x=.678, y=.585, width=.302, height=.085,
                               name="screen_excluded", column="exclusion", lines=[
        "Records excluded",
        f"n = {c['records_excluded_before_report_assessment']}",
    ])
    sought = add_flow_box(ax, fig, x=.185, y=.445, width=.47, height=.085,
                          name="reports_sought", column="main", lines=[
        "Reports sought for retrieval",
        f"n = {c['reports_sought_for_retrieval']}",
    ])
    not_retrieved = add_flow_box(ax, fig, x=.678, y=.445, width=.302, height=.085,
                                 name="reports_not_retrieved", column="exclusion", lines=[
        "Reports not retrieved",
        f"n = {c['reports_not_retrieved']}",
    ])
    assessed = add_flow_box(ax, fig, x=.185, y=.305, width=.47, height=.085,
                            name="eligibility", column="main", lines=[
        "Reports assessed for eligibility",
        f"n = {c['reports_assessed_for_eligibility']}",
    ])
    excluded = add_flow_box(ax, fig, x=.678, y=.155, width=.302, height=.235,
                            name="eligibility_excluded", column="exclusion", lines=[
        f"Reports excluded (n = {c['reports_excluded']})",
        "",
        f"Non-P17 contrast (n = {c['reports_excluded_no_P17_contrast']})",
        f"Insufficient replication (n = {c['reports_excluded_insufficient_replication']})",
        f"Wrong tissue/modality (n = {c['reports_excluded_wrong_tissue_or_modality']})",
        "No reconstructable",
        f"expression/contrast (n = {c['reports_excluded_no_recoverable_expression_or_contrast']})",
        f"Other eligibility reasons (n = {c['reports_excluded_other']})",
    ])
    included = add_flow_box(ax, fig, x=.185, y=.030, width=.47, height=.130,
                            name="included", column="main", fill="#EEF5F1", lines=[
        "Cohorts included in meta-analysis",
        f"n = {c['unique_eligible_cohorts_included']}",
        "",
        f"GEO cohorts (n = {c['eligible_GEO_cohorts']})",
        f"Non-GEO cohorts (n = {c['eligible_non_GEO_cohorts']})",
    ])

    orthogonal_arrow(ax, fig, identified, removed, "identified_to_removed", [
        (identified["x"] + identified["w"], .865), (removed["x"], .865)])
    orthogonal_arrow(ax, fig, identified, screened, "identified_to_screened")
    orthogonal_arrow(ax, fig, screened, screen_excl, "screened_to_excluded", [
        (screened["x"] + screened["w"], screened["y"] + screened["h"] / 2),
        (screen_excl["x"], screen_excl["y"] + screen_excl["h"] / 2)])
    orthogonal_arrow(ax, fig, screened, sought, "screened_to_sought")
    orthogonal_arrow(ax, fig, sought, not_retrieved, "sought_to_not_retrieved", [
        (sought["x"] + sought["w"], sought["y"] + sought["h"] / 2),
        (not_retrieved["x"], not_retrieved["y"] + not_retrieved["h"] / 2)])
    orthogonal_arrow(ax, fig, sought, assessed, "sought_to_assessed")
    orthogonal_arrow(ax, fig, assessed, excluded, "assessed_to_excluded", [
        (assessed["x"] + assessed["w"], assessed["y"] + assessed["h"] / 2),
        (excluded["x"], assessed["y"] + assessed["h"] / 2)])
    orthogonal_arrow(ax, fig, assessed, included, "assessed_to_included", [
        (assessed["x"] + assessed["w"] / 2, assessed["y"]),
        (assessed["x"] + assessed["w"] / 2, .235),
        (included["x"] + included["w"] / 2, .235),
        (included["x"] + included["w"] / 2, included["y"] + included["h"]),
    ])
    return fig


eff = pd.read_csv(ROOT / "02_mouse_unit_audit/MOUSE_FINAL_PRIMARY_STUDY_EFFECTS.tsv", sep="\t")
comp = pd.read_csv(ROOT / "02_mouse_compartment_audit/MOUSE_COMPARTMENT_REGISTRY.tsv", sep="\t")
eff = eff.merge(comp[["accession", "compartment_category"]], left_on="study", right_on="accession")
meta = pd.read_csv(ROOT / "02_mouse_unit_audit/MOUSE_FINAL_PRIMARY_META_RESULTS.tsv", sep="\t").iloc[0]
strict = pd.read_csv(ROOT / "02_mouse_unit_audit/MOUSE_STRICT_UNIT_CONFIRMED_META_RESULTS.tsv", sep="\t").iloc[0]
loo = pd.read_csv(ROOT / "02_mouse_unit_audit/MOUSE_FINAL_PRIMARY_LOO.tsv", sep="\t")
sub = pd.read_csv(ROOT / "02_mouse_compartment_audit/MOUSE_COMPARTMENT_META_RESULTS.tsv", sep="\t")


def draw_fig2a(ax):
    e = eff.sort_values("hedges_g").reset_index(drop=True)
    markers = {
        "WHOLE_RETINA_OR_LYSATE": "o",
        "ENRICHED_OR_ISOLATED_CELL_COMPARTMENT": "s",
        "OTHER_RETINAL_COMPARTMENT": "D",
    }
    for category, marker in markers.items():
        z = e[e.compartment_category == category]
        ax.errorbar(
            z.hedges_g, z.index,
            xerr=[z.hedges_g - z.ci95_low, z.ci95_high - z.hedges_g],
            fmt=marker, color=DARK, ecolor=DARK, capsize=2,
            ms=STYLE["marker_size_pt"], lw=STYLE["line_width_pt"], zorder=3,
        )
    pooled_y = -1.25
    ax.errorbar(
        meta.pooled_hedges_g, pooled_y,
        xerr=[[meta.pooled_hedges_g - meta.ci95_low], [meta.ci95_high - meta.pooled_hedges_g]],
        fmt="D", color=ORANGE, ecolor=ORANGE, capsize=3, ms=5.5, zorder=4,
    )
    ax.plot(
        [meta.prediction_interval_low, meta.prediction_interval_high],
        [pooled_y - .35, pooled_y - .35], color=ORANGE, lw=3.0, alpha=.55, zorder=2,
    )
    reference_line(ax)
    ax.set_yticks([*range(len(e)), pooled_y], [*e.study, "REML/HK pooled"])
    ax.set_ylim(pooled_y - .75, len(e) - .2)
    ax.set_xlim(-2.0, 6.2)
    ax.set_xticks([-2, 0, 2, 4, 6])
    ax.set_xlabel("Hedges g (OIR versus normoxia)")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=DARK, markeredgecolor=DARK,
               label="Whole retina/lysate", markersize=5),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=DARK, markeredgecolor=DARK,
               label="Enriched/isolated cells", markersize=5),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=DARK, markeredgecolor=DARK,
               label="Other retinal compartment", markersize=5),
    ]
    legend = ax.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(.5, -.13), ncol=3,
        frameon=True, facecolor="white", edgecolor=LIGHT, columnspacing=1.4, handletextpad=.5,
    )
    return legend


def draw_fig2b(ax):
    z = loo.sort_values("pooled_hedges_g").reset_index(drop=True)
    yy = np.arange(len(z))
    ax.errorbar(
        z.pooled_hedges_g, yy,
        xerr=[z.pooled_hedges_g - z.ci95_low, z.ci95_high - z.pooled_hedges_g],
        fmt="o", color=DARK, ecolor=DARK, capsize=2, ms=3.5, lw=.8,
    )
    reference_line(ax, float(meta.pooled_hedges_g), linestyle="--", linewidth=1.0).set_color(ORANGE)
    reference_line(ax)
    ax.set_yticks(yy, z.omitted_accession)
    ax.set_xlabel("Pooled Hedges g after omission")


def draw_fig2c(ax):
    order = ["BROAD_ALL_RETINAL_CONTEXTS", "WHOLE_RETINA_OR_LYSATE", "ENRICHED_OR_ISOLATED_CELL_COMPARTMENT"]
    z = sub.set_index("analysis_set").loc[order].reset_index()
    yy = np.arange(3)[::-1]
    labels = ["Broad retinal\ncontexts", "Whole retina /\nlysate", "Enriched /\nisolated cells"]
    for i, row in z.iterrows():
        ax.errorbar(
            row.pooled_hedges_g, yy[i],
            xerr=[[row.pooled_hedges_g - row.ci95_low], [row.ci95_high - row.pooled_hedges_g]],
            fmt="D", color=DARK, ecolor=DARK, capsize=3, ms=4.8,
        )
        ax.text(row.ci95_high + .08, yy[i], f"k={int(row.k)}", va="center", fontsize=STYLE["annotation_pt"])
    reference_line(ax)
    ax.set_yticks(yy, labels)
    ax.set_xlim(-1.45, 2.75)
    ax.set_xlabel("Pooled Hedges g (95% CI)")


def draw_fig2d(ax):
    ax.errorbar(
        strict.pooled_hedges_g, 0,
        xerr=[[strict.pooled_hedges_g - strict.ci95_low], [strict.ci95_high - strict.pooled_hedges_g]],
        fmt="D", color=DARK, ecolor=DARK, capsize=3, ms=4.8,
    )
    reference_line(ax)
    ax.set_yticks([0], ["GSE234447 +\nGSE315511"])
    ax.set_ylim(-.7, .7)
    ax.set_xlim(min(-27, strict.ci95_low - 1), max(30, strict.ci95_high + 1))
    ax.set_xlabel("Pooled Hedges g (95% CI)")
    ax.text(.97, .84, "k=2", transform=ax.transAxes, ha="right", va="top", fontsize=STYLE["annotation_pt"])


def build_fig2() -> plt.Figure:
    fig = plt.figure(figsize=(CONTRACT["output"]["width_in"], CONTRACT["figures"]["Fig2"]["height_in"]))
    gs = fig.add_gridspec(
        2, 3, height_ratios=[1.78, 1.0], width_ratios=[1.22, 1.05, .92],
        hspace=.60, wspace=.72, left=.165, right=.945, top=.935, bottom=.105,
    )
    ax_a = fig.add_subplot(gs[0, :]); panel_label(ax_a, "A"); panel_title(ax_a, "A", "Cohort-level effects and pooled estimate")
    draw_fig2a(ax_a)
    ax_b = fig.add_subplot(gs[1, 0]); panel_label(ax_b, "B"); panel_title(ax_b, "B", "Leave-one-cohort-out estimates")
    draw_fig2b(ax_b)
    ax_c = fig.add_subplot(gs[1, 1]); panel_label(ax_c, "C"); panel_title(ax_c, "C", "Compartment sensitivity")
    draw_fig2c(ax_c)
    ax_d = fig.add_subplot(gs[1, 2]); panel_label(ax_d, "D"); panel_title(ax_d, "D", "Strict-unit sensitivity")
    draw_fig2d(ax_d)
    return fig


h = pd.read_csv(ROOT / "04_human_rebuild/human_primary_model_input.tsv", sep="\t")
h = h[h.analysis_universe == "HIGHER_CELL_COUNT_7"].copy()
mods = pd.read_csv(ROOT / "04_human_rebuild/human_primary_model_full_results.tsv", sep="\t")
mods = mods[(mods.analysis_universe == "HIGHER_CELL_COUNT_7") & (mods.term == "ADORA2A_z")]
lodo = pd.read_csv(ROOT / "04_human_rebuild/human_LODO.tsv", sep="\t")
lodo_m0 = lodo[(lodo.analysis_universe == "HIGHER_CELL_COUNT_7") & (lodo.model == "M0")]


def draw_fig3a(ax, fig=None):
    for dataset, color, marker in [("GSE165784", BLUE, "o"), ("GSE245561", ORANGE, "s")]:
        z = h[h.dataset == dataset]
        ax.scatter(z.ADORA2A_z, z.score_z, label=dataset, color=color, marker=marker,
                   s=38, edgecolor="white", linewidth=.45, zorder=3)
    ax.set_xlabel("ADORA2A abundance (z)")
    ax.set_ylabel("Six-gene donor score (sensitivity)")
    legend = ax.legend(
        loc="lower right", ncol=1, frameon=True, facecolor="white", edgecolor=LIGHT,
        framealpha=.96, borderpad=.5, labelspacing=.4,
    )
    if fig is not None:
        register_legend_data(fig, legend, ax, h.ADORA2A_z, h.score_z, "Fig3A donors", radius_pt=5)


def draw_fig3b(ax):
    z = mods.set_index("model").loc[["M0", "M1", "M2"]].reset_index()
    yy = np.arange(3)[::-1]
    ax.errorbar(
        z.estimate, yy, xerr=[z.estimate - z.ci95_low, z.ci95_high - z.estimate],
        fmt="D", color=DARK, ecolor=DARK, capsize=3, ms=4.8,
    )
    reference_line(ax)
    ax.set_yticks(yy, ["M0", "M1", "M2"])
    ax.set_xlabel("ADORA2A coefficient (95% CI)")


def draw_fig3c(ax):
    z = lodo_m0.sort_values("estimate")
    yy = np.arange(len(z))
    ax.errorbar(
        z.estimate, yy, xerr=[z.estimate - z.ci95_low, z.ci95_high - z.estimate],
        fmt="o", color=DARK, ecolor=DARK, capsize=2, ms=4,
    )
    reference_line(ax)
    ax.set_yticks(yy, z.left_out_GSM)
    ax.set_xlabel("M0 coefficient after donor omission")


def draw_fig3d(ax):
    ax.set_axis_off()
    z = mods.set_index("model").loc[["M0", "M1", "M2"]].reset_index()
    diag = [[r.model, int(r.residual_df), f"{r.condition_number:.1f}", f"{r.vif:.1f}"] for _, r in z.iterrows()]
    table = ax.table(
        cellText=diag, colLabels=["Model", "Residual\ndf", "Condition\nno.", "VIF"],
        bbox=[.02, .28, .96, .58], cellLoc="center", colColours=["#DDEBF7"] * 4,
    )
    table.auto_set_font_size(False); table.set_fontsize(8)
    for cell in table.get_celld().values():
        cell.set_text_props(va="center")
        cell.set_edgecolor(LIGHT)
        cell.set_linewidth(.6)
    ax.text(.5, .14, "M2: technical sensitivity", ha="center", va="center",
            transform=ax.transAxes, color=GREY, fontsize=STYLE["annotation_pt"])


def build_fig3() -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(CONTRACT["output"]["width_in"], CONTRACT["figures"]["Fig3"]["height_in"]))
    fig.subplots_adjust(left=.13, right=.945, top=.93, bottom=.105, hspace=.56, wspace=.52)
    titles = {
        "A": "Higher-cell-count PDR donors\n(≥20 endothelial cells)",
        "B": "Sensitivity-model estimates",
        "C": "M0 leave-one-donor-out",
        "D": "Identification diagnostics",
    }
    for ax, label in zip(axes.ravel(), "ABCD"):
        panel_label(ax, label); panel_title(ax, label, titles[label], pad=9)
    draw_fig3a(axes[0, 0], fig)
    draw_fig3b(axes[0, 1])
    draw_fig3c(axes[1, 0])
    draw_fig3d(axes[1, 1])
    return fig


def draw_fig4a(ax):
    z = pd.read_csv(ROOT / "05_contextual_human/Context_GSE160306.tsv", sep="\t")
    yy = np.arange(len(z))
    ax.errorbar(z.estimate, yy, xerr=[z.estimate - z.ci95_low, z.ci95_high - z.estimate],
                fmt="o", color=DARK, ecolor=DARK, capsize=3, ms=4)
    reference_line(ax)
    ax.set_yticks(yy, ["Stage", "Stage × peripheral\nregion"])
    ax.set_xlabel("Coefficient (95% CI)")


def draw_fig4b(ax):
    z = pd.read_csv(ROOT / "05_contextual_human/Context_GSE60436_samples.tsv", sep="\t")
    order = [x for x in ["normal_retina", "inactive_FVM", "active_FVM"] if x in set(z.group)]
    for i, group in enumerate(order):
        values = z.loc[z.group == group, "ADORA2A"].dropna()
        jitter = np.linspace(-.06, .06, len(values)) if len(values) > 1 else np.zeros(len(values))
        ax.scatter(np.full(len(values), i) + jitter, values, color=[GREY, BLUE, ORANGE][i],
                   s=28, edgecolor="white", linewidth=.4, zorder=3)
    ax.set_xticks(range(len(order)), [x.replace("_", "\n") for x in order])
    ax.set_ylabel("Deposited expression value")
    ax.margins(y=.16)


def draw_fig4c(ax):
    z = pd.read_csv(ROOT / "05_contextual_human/Context_GSE276892.tsv", sep="\t")
    groups = list(dict.fromkeys(z.group))
    for i, group in enumerate(groups):
        values = z.loc[z.group == group, "ADORA2A"].dropna()
        jitter = np.linspace(-.08, .08, len(values)) if len(values) > 1 else np.zeros(len(values))
        ax.scatter(np.full(len(values), i) + jitter, values, color=[BLUE, ORANGE, GREEN, MAGENTA][i % 4],
                   s=27, edgecolor="white", linewidth=.4, zorder=3)
        if len(values):
            ax.plot([i - .15, i + .15], [np.median(values)] * 2, color=DARK, zorder=4)
    ax.set_xticks(range(len(groups)), groups)
    ax.set_ylabel("ADORA2A expression")


def draw_fig4d(ax):
    z = pd.read_csv(ROOT / "05_contextual_human/Context_CD31_FVM_effects.tsv", sep="\t")
    yy = np.arange(len(z))
    ax.errorbar(z.effect, yy, xerr=[z.effect - z.lo, z.hi - z.effect],
                fmt="o", color=DARK, ecolor=DARK, capsize=3, ms=4)
    reference_line(ax)
    labels = [
        "GSE94019\nPDR FVM vs control",
        "GSE307925\ntreated vs untreated",
        "GSE179568\nPDR vs MP",
        "GSE179568\nPDR vs MH",
    ]
    ax.set_yticks(yy, labels)
    ax.set_xlabel("Dataset-specific coefficient (95% CI)")


def draw_fig4e(ax):
    z = pd.read_csv(ROOT / "05_contextual_human/Context_2026_MNV.tsv", sep="\t")
    z = z[z.file == "media-2.xlsx"].drop_duplicates("cell_type").sort_values("logFC_pseudobulk")
    yy = np.arange(len(z))
    ax.scatter(z.logFC_pseudobulk, yy, color=ORANGE, s=28, edgecolor="white", linewidth=.4, zorder=3)
    reference_line(ax)
    ax.set_yticks(yy, z.cell_type)
    ax.set_xlabel("MNV vs control pseudobulk logFC")


def draw_fig4f(ax):
    z = pd.read_csv(ROOT / "05_contextual_human/Context_GSE234047.tsv", sep="\t")
    wide = z.pivot(index="section", columns="stratum", values="mean_count")
    if {"NON_LESIONAL", "LESION_OVERLAP"}.issubset(wide.columns):
        for _, row in wide.iterrows():
            ax.plot([0, 1], [row.NON_LESIONAL, row.LESION_OVERLAP], color=GREY, marker="o", ms=4.5)
    ax.set_xticks([0, 1], ["Non-lesional", "Lesion overlap"])
    ax.set_xlim(-.15, 1.15)
    ax.set_ylabel("Mean ADORA2A count per spot")


def build_fig4() -> plt.Figure:
    fig, axes = plt.subplots(3, 2, figsize=(CONTRACT["output"]["width_in"], CONTRACT["figures"]["Fig4"]["height_in"]))
    fig.subplots_adjust(left=.185, right=.97, top=.945, bottom=.075, hspace=.72, wspace=.62)
    titles = [
        "GSE160306 retina", "GSE60436 FVM", "GSE276892 vitreous cells",
        "CD31/FVM contexts", "2026 MNV preprint v1", "GSE234047 spatial",
    ]
    drawers = [draw_fig4a, draw_fig4b, draw_fig4c, draw_fig4d, draw_fig4e, draw_fig4f]
    for ax, label, title, drawer in zip(axes.ravel(), "ABCDEF", titles, drawers):
        panel_label(ax, label); panel_title(ax, label, title, pad=10); drawer(ax)
    return fig


def export_panel_previews() -> None:
    panel_dir = OUT / "panel_previews"
    panel_dir.mkdir(exist_ok=True)
    specifications = [
        ("Fig2A", draw_fig2a, "Cohort-level effects and pooled estimate", (5.8, 4.8)),
        ("Fig2B", draw_fig2b, "Leave-one-cohort-out estimates", (4.4, 4.8)),
        ("Fig2C", draw_fig2c, "Compartment sensitivity", (4.4, 3.4)),
        ("Fig2D", draw_fig2d, "Strict-unit sensitivity", (4.4, 3.0)),
        ("Fig3A", lambda ax: draw_fig3a(ax), "Higher-cell-count PDR donors", (4.4, 3.5)),
        ("Fig3B", draw_fig3b, "Sensitivity-model estimates", (4.4, 3.5)),
        ("Fig3C", draw_fig3c, "M0 leave-one-donor-out", (4.4, 3.8)),
        ("Fig3D", draw_fig3d, "Identification diagnostics", (4.4, 3.4)),
        ("Fig4A", draw_fig4a, "GSE160306 retina", (4.4, 3.2)),
        ("Fig4B", draw_fig4b, "GSE60436 FVM", (4.4, 3.2)),
        ("Fig4C", draw_fig4c, "GSE276892 vitreous cells", (4.4, 3.2)),
        ("Fig4D", draw_fig4d, "CD31/FVM contexts", (4.8, 3.5)),
        ("Fig4E", draw_fig4e, "2026 MNV preprint v1", (4.4, 3.2)),
        ("Fig4F", draw_fig4f, "GSE234047 spatial", (4.4, 3.2)),
    ]
    for stem, drawer, title, size in specifications:
        fig, ax = plt.subplots(figsize=size)
        drawer(ax); ax.set_title(title, loc="left", fontsize=9, fontweight="semibold", pad=8)
        fig.subplots_adjust(left=.22, right=.96, top=.88, bottom=.18)
        fig.savefig(panel_dir / f"{stem}.png", dpi=180, facecolor="white")
        plt.close(fig)


parser = argparse.ArgumentParser()
parser.add_argument("--fig1-only", action="store_true", help="Regenerate only the PRISMA flow figure and preserve Fig 2–4 files.")
args = parser.parse_args()

if args.fig1_only:
    fig1_qa = save_main_figure(build_fig1(), "Fig1", OUT)
    existing = OUT / "FIGURE_LAYOUT_QA_MAIN.tsv"
    if existing.exists():
        prior = pd.read_csv(existing, sep="\t")
        main_qa = pd.concat([fig1_qa, prior.loc[prior.figure != "Fig1"]], ignore_index=True)
    else:
        main_qa = fig1_qa
    main_qa.to_csv(existing, sep="\t", index=False)
    combined = OUT / "FIGURE_LAYOUT_QA.tsv"
    if combined.exists():
        prior_combined = pd.read_csv(combined, sep="\t")
        pd.concat([fig1_qa, prior_combined.loc[prior_combined.figure != "Fig1"]], ignore_index=True).to_csv(
            combined, sep="\t", index=False
        )
    print(f"Fig1 written to {OUT}; Fig2–4 files preserved")
else:
    qa_frames = [
        save_main_figure(build_fig1(), "Fig1", OUT),
        save_main_figure(build_fig2(), "Fig2", OUT),
        save_main_figure(build_fig3(), "Fig3", OUT),
        save_main_figure(build_fig4(), "Fig4", OUT),
    ]
    export_panel_previews()
    pd.concat(qa_frames, ignore_index=True).to_csv(OUT / "FIGURE_LAYOUT_QA_MAIN.tsv", sep="\t", index=False)
    print(f"Four main figures written to {OUT}")
