#!/usr/bin/env python3
"""Regenerate the four adjudicated PLOS ONE main figures from audited tables."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

from figure_style import (
    BLUE, ORANGE, GREEN, MAGENTA, GREY, LIGHT, DARK, CONTRACT, STYLE,
    panel_label, panel_title, reference_line, register_box_text,
    register_legend_data, save_main_figure,
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
    return counts


def add_flow_box(ax, fig, *, x: float, top: float, width: float, lines: list[str], edge: str, name: str):
    line_height = 0.027
    height = 0.042 + line_height * len(lines)
    y = top - height
    patch = FancyBboxPatch(
        (x, y), width, height, boxstyle="round,pad=0.010",
        facecolor="white", edgecolor=edge, linewidth=1.0,
    )
    ax.add_patch(patch)
    text = ax.text(
        x + 0.020, y + height / 2, "\n".join(lines),
        ha="left", va="center", multialignment="left",
        fontsize=STYLE["annotation_pt"], linespacing=1.23,
    )
    register_box_text(fig, patch, text, name)
    return {"x": x, "y": y, "w": width, "h": height, "patch": patch, "text": text}


def arrow(ax, start, end):
    ax.annotate(
        "", xy=end, xytext=start,
        arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.0,
                        shrinkA=0, shrinkB=0, mutation_scale=10),
        annotation_clip=False,
    )


def build_fig1() -> plt.Figure:
    c = read_counts()
    fig, ax = plt.subplots(figsize=(CONTRACT["output"]["width_in"], CONTRACT["figures"]["Fig1"]["height_in"]))
    ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    left = add_flow_box(ax, fig, x=.045, top=.965, width=.42, edge=BLUE, name="database_identification", lines=[
        "Databases and registers",
        f"GEO objects (n = {c['records_identified_GEO']})",
        f"PubMed reports (n = {c['records_identified_PubMed']})",
    ])
    right = add_flow_box(ax, fig, x=.535, top=.965, width=.42, edge=BLUE, name="other_identification", lines=[
        "Other sources",
        f"BioProject/SRA records (n = {c['records_identified_BioProject'] + c['records_identified_SRA']})",
        f"Europe PMC reports (n = {c['records_identified_Europe_PMC']})",
        f"BioStudies records (n = {c['records_identified_BioStudies']})",
    ])
    norm = add_flow_box(ax, fig, x=.17, top=.745, width=.66, edge=GREY, name="normalization", lines=[
        "Record normalization and deduplication",
        f"Duplicate/component records removed (n = {c['duplicates_or_component_records_removed']})",
        f"GEO objects consolidated to unique series (n = {c['unique_GEO_series']})",
        "Cross-source duplicates resolved",
    ])
    screened = add_flow_box(ax, fig, x=.075, top=.555, width=.46, edge=DARK, name="screened", lines=[
        "Records screened",
        f"n = {c['records_screened']}",
    ])
    screen_excl = add_flow_box(ax, fig, x=.63, top=.555, width=.325, edge=GREY, name="screen_excluded", lines=[
        "Records excluded",
        f"n = {c['records_excluded_before_report_assessment']}",
    ])
    assessed = add_flow_box(ax, fig, x=.075, top=.385, width=.46, edge=DARK, name="eligibility", lines=[
        f"Reports sought and retrieved (n = {c['reports_sought_for_retrieval']})",
        f"Reports assessed for eligibility (n = {c['reports_assessed_for_eligibility']})",
        f"Reports not retrieved (n = {c['reports_not_retrieved']})",
    ])
    excluded = add_flow_box(ax, fig, x=.585, top=.405, width=.37, edge=GREY, name="eligibility_excluded", lines=[
        f"Reports/datasets excluded (n = {c['reports_excluded']})",
        f"Non-P17 contrast (n = {c['reports_excluded_no_P17_contrast']})",
        f"Insufficient replication (n = {c['reports_excluded_insufficient_replication']})",
        f"Wrong tissue/modality (n = {c['reports_excluded_wrong_tissue_or_modality']})",
        f"No reconstructable expression /\ncontrast (n = {c['reports_excluded_no_recoverable_expression_or_contrast']})",
        f"Other eligibility reasons (n = {c['reports_excluded_other']})",
    ])
    included = add_flow_box(ax, fig, x=.20, top=.165, width=.60, edge=GREEN, name="included", lines=[
        "Unique eligible mouse cohorts",
        f"GEO cohorts (n = {c['eligible_GEO_cohorts']})",
        f"Non-GEO cohorts (n = {c['eligible_non_GEO_cohorts']})",
        f"Total (k = {c['unique_eligible_cohorts_included']})",
    ])

    arrow(ax, (left["x"] + left["w"] / 2, left["y"]), (norm["x"] + norm["w"] * .36, norm["y"] + norm["h"]))
    arrow(ax, (right["x"] + right["w"] / 2, right["y"]), (norm["x"] + norm["w"] * .64, norm["y"] + norm["h"]))
    arrow(ax, (norm["x"] + norm["w"] / 2, norm["y"]), (screened["x"] + screened["w"] / 2, screened["y"] + screened["h"]))
    arrow(ax, (screened["x"] + screened["w"], screened["y"] + screened["h"] / 2), (screen_excl["x"], screen_excl["y"] + screen_excl["h"] / 2))
    arrow(ax, (screened["x"] + screened["w"] / 2, screened["y"]), (assessed["x"] + assessed["w"] / 2, assessed["y"] + assessed["h"]))
    arrow(ax, (assessed["x"] + assessed["w"], assessed["y"] + assessed["h"] / 2), (excluded["x"], excluded["y"] + excluded["h"] / 2))
    arrow(ax, (assessed["x"] + assessed["w"] / 2, assessed["y"]), (included["x"] + included["w"] / 2, included["y"] + included["h"]))
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


qa_frames = []
qa_frames.append(save_main_figure(build_fig1(), "Fig1", OUT))
qa_frames.append(save_main_figure(build_fig2(), "Fig2", OUT))
qa_frames.append(save_main_figure(build_fig3(), "Fig3", OUT))
qa_frames.append(save_main_figure(build_fig4(), "Fig4", OUT))
export_panel_previews()
pd.concat(qa_frames, ignore_index=True).to_csv(OUT / "FIGURE_LAYOUT_QA_MAIN.tsv", sep="\t", index=False)
print(f"Four main figures written to {OUT}")
