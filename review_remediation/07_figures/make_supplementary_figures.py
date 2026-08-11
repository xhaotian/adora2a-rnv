#!/usr/bin/env python3
"""Generate the five non-redundant supplementary figures with measured QA."""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

from figure_style import (
    BLUE, ORANGE, GREEN, GREY, LIGHT, DARK, STYLE,
    panel_label, panel_title, reference_line, save_supplement_page,
)


PROJECT = Path(__file__).resolve().parents[2]
ROOT = PROJECT / "review_remediation"
OUT = ROOT / "07_figures"
OUT.mkdir(parents=True, exist_ok=True)


qa_frames = []
with PdfPages(OUT / "Supplementary_Figures.pdf") as pdf:
    # New S1 (former S2): selected sample-level examples.
    data = pd.read_csv(ROOT / "02_mouse_unit_audit/MOUSE_FINAL_SAMPLE_EXPRESSION.tsv", sep="\t")
    required = ["GSE234447", "GSE241239", "GSE315511"]
    ranked = list(data.groupby("study").size().sort_values(ascending=False).index)
    selected = required + [study for study in ranked if study not in required][:5]
    data = data[data.study.isin(selected)].copy()
    data["display_group"] = np.where(
        data.group.str.lower().str.startswith("norm"), "Normoxia",
        np.where(data.group.str.upper().str.startswith("OIR"), "OIR", data.group),
    )
    data["display_z"] = data.groupby("study").Adora2a.transform(
        lambda values: (values - values.mean()) / (values.std(ddof=0) or 1)
    )
    fig, axes = plt.subplots(2, 4, figsize=(7.3, 5.65))
    fig.subplots_adjust(left=.09, right=.985, top=.92, bottom=.13, hspace=.52, wspace=.42)
    for ax, label, study in zip(axes.ravel(), "ABCDEFGH", selected):
        panel_label(ax, label); panel_title(ax, label, study, pad=8)
        subset = data[data.study == study]
        for index, group in enumerate(["Normoxia", "OIR"]):
            values = subset.loc[subset.display_group == group, "display_z"].dropna().to_numpy()
            jitter = np.linspace(-.10, .10, len(values)) if len(values) > 1 else np.zeros(len(values))
            ax.scatter(
                np.full(len(values), index) + jitter, values,
                color=[BLUE, ORANGE][index], s=26, alpha=.86,
                edgecolor="white", linewidth=.4, zorder=4,
            )
            if len(values):
                ax.plot([index - .16, index + .16], [np.median(values)] * 2,
                        color=DARK, linewidth=1.2, zorder=3)
        ax.axhline(0, color=LIGHT, linewidth=.7, zorder=0)
        ax.set_xticks([0, 1], ["Normoxia", "OIR"])
    fig.supxlabel("Deposited retinal sample/library units", fontsize=STYLE["axis_label_pt"])
    fig.supylabel("Within-study standardized ADORA2A expression", fontsize=STYLE["axis_label_pt"])
    qa_frames.append(save_supplement_page(fig, pdf, "S1", OUT))

    # New S2 (former S3): higher-cell-count versus expanded donor sensitivity.
    models = pd.read_csv(ROOT / "04_human_rebuild/human_primary_model_full_results.tsv", sep="\t")
    models = models[models.term == "ADORA2A_z"].copy()
    fig, ax = plt.subplots(figsize=(7.3, 4.35))
    markers = {"HIGHER_CELL_COUNT_7": "o", "EXPANDED_DONOR_9": "s"}
    colors = {"HIGHER_CELL_COUNT_7": BLUE, "EXPANDED_DONOR_9": ORANGE}
    for offset, universe in [(-.10, "HIGHER_CELL_COUNT_7"), (.10, "EXPANDED_DONOR_9")]:
        z = models[models.analysis_universe == universe].set_index("model").loc[["M0", "M1", "M2"]]
        yy = np.arange(3)[::-1] + offset
        ax.errorbar(
            z.estimate, yy, xerr=[z.estimate - z.ci95_low, z.ci95_high - z.estimate],
            fmt=markers[universe], color=colors[universe], ecolor=colors[universe],
            capsize=3, ms=5, linewidth=.9,
        )
    reference_line(ax)
    ax.set_yticks(np.arange(3)[::-1], ["M0", "M1", "M2"])
    ax.set_xlim(-10.5, 8.5)
    ax.set_xticks([-10, -5, 0, 5])
    ax.set_xlabel("ADORA2A coefficient (95% CI)")
    handles = [
        Line2D([0], [0], color=BLUE, marker="o", label="Higher-cell-count stratum (≥20 cells)"),
        Line2D([0], [0], color=ORANGE, marker="s", label="Expanded stratum (including <20 cells)"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, 1.16), ncol=2,
              frameon=True, facecolor="white", edgecolor=LIGHT)
    fig.subplots_adjust(left=.17, right=.97, top=.78, bottom=.20)
    qa_frames.append(save_supplement_page(fig, pdf, "S2", OUT))

    # New S3 (former S4): full LODO with public GSM identifiers.
    lodo = pd.read_csv(ROOT / "04_human_rebuild/human_LODO.tsv", sep="\t")
    lodo = lodo[lodo.analysis_universe == "HIGHER_CELL_COUNT_7"]
    donors = list(dict.fromkeys(lodo.left_out_GSM))
    positions = {donor: index for index, donor in enumerate(donors)}
    fig, axes = plt.subplots(1, 3, figsize=(7.3, 5.0), sharey=True)
    fig.subplots_adjust(left=.16, right=.98, top=.88, bottom=.16, wspace=.36)
    for ax, label, model in zip(axes, "ABC", ["M0", "M1", "M2"]):
        panel_label(ax, label); panel_title(ax, label, model, pad=8)
        z = lodo[lodo.model == model].copy()
        yy = z.left_out_GSM.map(positions)
        ax.errorbar(
            z.estimate, yy, xerr=[z.estimate - z.ci95_low, z.ci95_high - z.estimate],
            fmt="o", color=DARK, ecolor=DARK, capsize=2, ms=4,
        )
        reference_line(ax)
        limits = {"M0": (-1.5, 2.0), "M1": (-15, 15), "M2": (-25, 22)}
        ticks = {"M0": [-1, 0, 1, 2], "M1": [-15, 0, 15], "M2": [-20, 0, 20]}
        ax.set_xlim(*limits[model])
        ax.set_xticks(ticks[model])
        ax.set_xlabel("Coefficient (95% CI)")
    axes[0].set_yticks(range(len(donors)), donors)
    qa_frames.append(save_supplement_page(fig, pdf, "S3", OUT))

    # New S4 (former S5): all disclosed signatures with reader-facing labels.
    signatures = pd.read_csv(ROOT / "04_human_rebuild/human_signature_sensitivity_results.tsv", sep="\t")
    signatures = signatures[signatures.analysis_universe == "HIGHER_CELL_COUNT_7"].copy()
    label_map = {
        "HISTORICAL_SIX_GENE": "Historical six-gene score",
        "EXTERNAL_VEGF_RESPONSE": "External VEGF response",
        "HALLMARK_ANGIOGENESIS": "Hallmark angiogenesis",
        "REACTOME_VEGFA_VEGFR2": "Reactome VEGFA–VEGFR2",
    }
    signatures["reader_label"] = signatures.signature.map(label_map) + " · " + signatures.model
    signatures["signature_order"] = signatures.signature.map({key: i for i, key in enumerate(label_map)})
    signatures["model_order"] = signatures.model.map({"M0": 0, "M1": 1, "M2": 2})
    signatures = signatures.sort_values(["signature_order", "model_order"], ascending=[False, False])
    fig, ax = plt.subplots(figsize=(7.3, 6.3))
    yy = np.arange(len(signatures))
    ax.errorbar(
        signatures.estimate, yy,
        xerr=[signatures.estimate - signatures.ci95_low, signatures.ci95_high - signatures.estimate],
        fmt="o", color=DARK, ecolor=DARK, capsize=2, ms=4,
    )
    reference_line(ax)
    ax.set_yticks(yy, signatures.reader_label)
    ax.set_xlim(-10.5, 8.5)
    ax.set_xticks([-10, -5, 0, 5])
    ax.set_xlabel("ADORA2A coefficient (95% CI)")
    fig.subplots_adjust(left=.34, right=.97, top=.94, bottom=.12)
    qa_frames.append(save_supplement_page(fig, pdf, "S4", OUT))

    # New S5 (former S6): exploratory small-study-effect diagnostic.
    effects = pd.read_csv(ROOT / "02_mouse_unit_audit/MOUSE_FINAL_PRIMARY_STUDY_EFFECTS.tsv", sep="\t")
    test = pd.read_csv(ROOT / "02_mouse_compartment_audit/FUNNEL_ASYMMETRY_DIAGNOSTIC.tsv", sep="\t").iloc[0]
    effects["se"] = np.sqrt(effects.sampling_variance)
    fig, ax = plt.subplots(figsize=(7.3, 5.0))
    ax.scatter(effects.hedges_g, effects.se, color=DARK, s=32, edgecolor="white", linewidth=.4, zorder=3)
    ax.axvline(effects.hedges_g.mean(), color=GREY, linestyle="--", linewidth=.9)
    ax.invert_yaxis()
    ax.set_xlim(-.75, 4.5)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_ylim(1.35, .35)
    ax.set_yticks([.4, .6, .8, 1.0, 1.2])
    ax.set_xlabel("Hedges g")
    ax.set_ylabel("Standard error")
    ax.text(
        .98, .97,
        f"Pustejovsky–Rodgers\nt = {test.test_statistic_t:.3f}; P = {test.p:.3f}; k = {int(test.k)}",
        transform=ax.transAxes, ha="right", va="top", fontsize=STYLE["annotation_pt"],
        color=GREY, bbox=dict(facecolor="white", edgecolor=LIGHT, boxstyle="round,pad=.3"),
    )
    fig.subplots_adjust(left=.15, right=.96, top=.94, bottom=.16)
    qa_frames.append(save_supplement_page(fig, pdf, "S5", OUT))


main_qa = pd.read_csv(OUT / "FIGURE_LAYOUT_QA_MAIN.tsv", sep="\t")
all_qa = pd.concat([main_qa, *qa_frames], ignore_index=True)
all_qa.to_csv(OUT / "FIGURE_LAYOUT_QA.tsv", sep="\t", index=False)
if all_qa.status.eq("FAIL").any():
    raise RuntimeError("FIGURE_LAYOUT_QA contains failures")
print(OUT / "Supplementary_Figures.pdf")
