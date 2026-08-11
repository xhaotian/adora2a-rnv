#!/usr/bin/env python3
"""Generate five PLOS ONE production figures from audited tabular outputs."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle
from PIL import Image

PROJECT = Path(__file__).resolve().parents[2]
ROOT = PROJECT / "review_remediation"
OUT = ROOT / "07_figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, GREEN, MAGENTA = "#0072B2", "#D55E00", "#009E73", "#CC79A7"
GREY, LIGHT, DARK = "#6F6F6F", "#D4D4D4", "#222222"
plt.rcParams.update({
    "font.family": "Arimo", "font.size": 8,
    "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 8,
    "ytick.labelsize": 8, "legend.fontsize": 8, "axes.spines.top": False,
    "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42,
})


def panel(ax, label):
    ax.text(0.0, 1.04, label, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", ha="left")


def save(fig, stem):
    """Preserve physical canvas; write 600-dpi flattened RGB LZW TIFF."""
    fig.savefig(OUT / f"{stem}.pdf", facecolor="white")
    png = OUT / f".{stem}.png"
    fig.savefig(png, dpi=600, facecolor="white")
    with Image.open(png) as im:
        im.convert("RGB").save(OUT / f"{stem}.tif", compression="tiff_lzw", dpi=(600, 600))
    png.unlink()
    plt.close(fig)


# Fig 1: PRISMA 2020 flow. Counts come only from the all-source audit.
counts = pd.read_csv(ROOT / "03_systematic_search/PRISMA_FLOW_SOURCE.tsv", sep="\t").set_index("item")["count"]
def c(key): return int(counts.loc[key])

fig, ax = plt.subplots(figsize=(7.3, 8.45))
ax.set_axis_off()

def box(x, y, w, h, text, edge=DARK, align="left"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008",
                                facecolor="white", edgecolor=edge, linewidth=1.1))
    ax.text(x + (0.02 if align == "left" else w/2), y + h/2, text,
            ha=align, va="center", fontsize=8, linespacing=1.28)

box(.03, .77, .45, .20,
    "Records identified from databases\n"
    f"GEO objects (n = {c('records_identified_GEO')})\n"
    f"PubMed reports (n = {c('records_identified_PubMed')})\n"
    f"BioProject records (n = {c('records_identified_BioProject')})\n"
    f"SRA experiments (n = {c('records_identified_SRA')})\n"
    f"Europe PMC reports (n = {c('records_identified_Europe_PMC')})\n"
    f"BioStudies records (n = {c('records_identified_BioStudies')})", BLUE)
box(.57, .79, .40, .15,
    "Records removed before screening\n"
    f"Duplicate or component records (n = {c('duplicates_or_component_records_removed')})",
    GREY)
box(.03, .61, .45, .09, f"Records screened\n(n = {c('records_screened')})", DARK)
box(.57, .61, .40, .09,
    f"Records excluded\n(n = {c('records_excluded_before_report_assessment')})", GREY)
box(.03, .46, .45, .09,
    f"Reports sought for retrieval\n(n = {c('reports_sought_for_retrieval')})", DARK)
box(.57, .46, .40, .09,
    f"Reports not retrieved\n(n = {c('reports_not_retrieved')})", GREY)
box(.03, .31, .45, .09,
    f"Reports assessed for eligibility\n(n = {c('reports_assessed_for_eligibility')})", DARK)
box(.57, .25, .40, .21,
    "Reports excluded (n = 40)\n"
    f"No P17 contrast (n = {c('reports_excluded_no_P17_contrast')})\n"
    f"Insufficient replication (n = {c('reports_excluded_insufficient_replication')})\n"
    f"Wrong tissue or modality (n = {c('reports_excluded_wrong_tissue_or_modality')})\n"
    f"No recoverable expression/contrast (n = {c('reports_excluded_no_recoverable_expression_or_contrast')})\n"
    f"Other eligibility reason (n = {c('reports_excluded_other')})", GREY)
box(.03, .08, .45, .13,
    "Unique eligible cohorts included\n"
    f"GEO cohorts (n = {c('eligible_GEO_cohorts')})\n"
    f"Non-GEO cohorts (n = {c('eligible_non_GEO_cohorts')})\n"
    f"Total (k = {c('unique_eligible_cohorts_included')})", GREEN)
for y1, y2 in [(.77, .70), (.61, .55), (.46, .40), (.31, .21)]:
    ax.annotate("", xy=(.255, y2), xytext=(.255, y1),
                arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1))
for y in [.655, .505, .355]:
    ax.annotate("", xy=(.57, y), xytext=(.48, y),
                arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1))
ax.text(.03, .015, "Search date: 11 August 2026; records and decisions are listed in Search_All_Sources.tsv.",
        fontsize=8, color=GREY)
save(fig, "Fig1")


# Fig 2: mouse synthesis with compartment encoding.
eff = pd.read_csv(ROOT / "02_mouse_unit_audit/MOUSE_FINAL_PRIMARY_STUDY_EFFECTS.tsv", sep="\t")
comp = pd.read_csv(ROOT / "02_mouse_compartment_audit/MOUSE_COMPARTMENT_REGISTRY.tsv", sep="\t")
eff = eff.merge(comp[["accession", "compartment_category"]], left_on="study", right_on="accession")
meta = pd.read_csv(ROOT / "02_mouse_unit_audit/MOUSE_FINAL_PRIMARY_META_RESULTS.tsv", sep="\t").iloc[0]
loo = pd.read_csv(ROOT / "02_mouse_unit_audit/MOUSE_FINAL_PRIMARY_LOO.tsv", sep="\t")
sub = pd.read_csv(ROOT / "02_mouse_compartment_audit/MOUSE_COMPARTMENT_META_RESULTS.tsv", sep="\t")
fig = plt.figure(figsize=(7.3, 8.55))
gs = fig.add_gridspec(2, 2, height_ratios=[1.6, 1], hspace=.42, wspace=.42,
                      left=.20, right=.97, top=.97, bottom=.13)
ax = fig.add_subplot(gs[0, :]); panel(ax, "A")
e = eff.sort_values("hedges_g").reset_index(drop=True); y = np.arange(len(e))
styles = {
    "WHOLE_RETINA_OR_LYSATE": ("o", BLUE),
    "ENRICHED_OR_ISOLATED_CELL_COMPARTMENT": ("s", ORANGE),
    "OTHER_RETINAL_COMPARTMENT": ("^", GREEN),
}
for category, (marker, color) in styles.items():
    z = e[e.compartment_category == category]
    ax.errorbar(z.hedges_g, z.index, xerr=[z.hedges_g-z.ci95_low, z.ci95_high-z.hedges_g],
                fmt=marker, color=color, ecolor=color, capsize=2, ms=4.5, lw=.9)
py = -1.25
ax.errorbar(meta.pooled_hedges_g, py,
            xerr=[[meta.pooled_hedges_g-meta.ci95_low], [meta.ci95_high-meta.pooled_hedges_g]],
            fmt="D", color=DARK, capsize=3, ms=5)
ax.plot([meta.prediction_interval_low, meta.prediction_interval_high], [py-.35, py-.35],
        color=DARK, lw=3, alpha=.42)
ax.axvline(0, color=GREY, lw=.8)
ax.set_yticks(list(y)+[py], list(e.study)+["Pooled"])
ax.set_ylim(py-.72, len(e)-.2); ax.set_xlabel("Hedges g (OIR versus normoxia)")
ax.set_title(f"Average direction across retinal transcriptomic contexts: k={int(meta.k)}, "
             f"{int(meta.positive_studies)}/{int(meta.k)} positive; I²={meta.I2:.1f}%")
handles = [Line2D([0], [0], marker=m, color="none", markerfacecolor=col,
                  markeredgecolor=col, label=lab, markersize=5)
           for lab, (m, col) in [
               ("Whole retina/lysate", styles["WHOLE_RETINA_OR_LYSATE"]),
               ("Enriched/isolated cells", styles["ENRICHED_OR_ISOLATED_CELL_COMPARTMENT"]),
               ("Other retinal compartment", styles["OTHER_RETINAL_COMPARTMENT"]),
           ]]
ax.legend(handles=handles, frameon=False, loc="lower right")

ax = fig.add_subplot(gs[1, 0]); panel(ax, "B")
l = loo.sort_values("pooled_hedges_g").reset_index(drop=True); yy = np.arange(len(l))
ax.errorbar(l.pooled_hedges_g, yy,
            xerr=[l.pooled_hedges_g-l.ci95_low, l.ci95_high-l.pooled_hedges_g],
            fmt="o", color=BLUE, ecolor=BLUE, capsize=2, ms=3.5, lw=.8)
ax.axvline(meta.pooled_hedges_g, color=ORANGE, ls="--", lw=1)
ax.axvline(0, color=GREY, lw=.8); ax.set_yticks(yy, l.omitted_accession)
ax.set_xlabel("Pooled Hedges g after omission"); ax.set_title("Leave-one-accession-out")

ax = fig.add_subplot(gs[1, 1]); panel(ax, "C")
z = sub[sub.analysis_set.isin(["BROAD_ALL_RETINAL_CONTEXTS", "WHOLE_RETINA_OR_LYSATE",
                               "ENRICHED_OR_ISOLATED_CELL_COMPARTMENT"])].copy()
order = ["BROAD_ALL_RETINAL_CONTEXTS", "WHOLE_RETINA_OR_LYSATE",
         "ENRICHED_OR_ISOLATED_CELL_COMPARTMENT"]
z = z.set_index("analysis_set").loc[order].reset_index(); yy = np.arange(3)[::-1]
labels = ["All contexts (k=16)", "Whole retina/lysate (k=13)", "Enriched/isolated (k=2)"]
for i, r in z.iterrows():
    ax.errorbar(r.pooled_hedges_g, yy[i],
                xerr=[[r.pooled_hedges_g-r.ci95_low], [r.ci95_high-r.pooled_hedges_g]],
                fmt="D", color=[DARK, BLUE, ORANGE][i], capsize=3, ms=5)
ax.axvline(0, color=GREY, lw=.8); ax.set_yticks(yy, labels)
ax.set_xlabel("REML/HK pooled Hedges g (95% CI)"); ax.set_title("Compartment sensitivity")
fig.text(.20, .018, "Public animal/eye mapping was incomplete for several datasets; all cohorts were evaluated using the same\n"
         "deposited-sample eligibility criteria.", fontsize=8, color=GREY, va="bottom")
save(fig, "Fig2")


# Fig 3: donor-level human sensitivity analysis. No pooled regression line.
h = pd.read_csv(ROOT / "04_human_rebuild/human_primary_model_input.tsv", sep="\t")
h = h[h.analysis_universe == "HIGHER_CELL_COUNT_7"].copy()
mods = pd.read_csv(ROOT / "04_human_rebuild/human_primary_model_full_results.tsv", sep="\t")
mods = mods[(mods.analysis_universe == "HIGHER_CELL_COUNT_7") & (mods.term == "ADORA2A_z")]
lodo = pd.read_csv(ROOT / "04_human_rebuild/human_LODO.tsv", sep="\t")
lodo = lodo[(lodo.analysis_universe == "HIGHER_CELL_COUNT_7") & (lodo.model == "M0")]
fig, axs = plt.subplots(2, 2, figsize=(7.3, 7.25))
fig.subplots_adjust(left=.11, right=.97, top=.96, bottom=.10, hspace=.48, wspace=.42)
ax = axs[0, 0]; panel(ax, "A")
for ds, color, marker in [("GSE165784", BLUE, "o"), ("GSE245561", ORANGE, "s")]:
    d = h[h.dataset == ds]
    ax.scatter(d.ADORA2A_z, d.score_z, label=ds, color=color, marker=marker, s=34)
ax.set(xlabel="ADORA2A abundance (z)", ylabel="Six-gene donor score (sensitivity)",
       title="Higher-cell-count donor stratum (≥20 cells)")
ax.legend(frameon=False)

ax = axs[0, 1]; panel(ax, "B")
m = mods.set_index("model").loc[["M0", "M1", "M2"]].reset_index(); yy = np.arange(3)[::-1]
ax.errorbar(m.estimate, yy, xerr=[m.estimate-m.ci95_low, m.ci95_high-m.estimate],
            fmt="D", color=BLUE, ecolor=BLUE, capsize=3, ms=5)
ax.axvline(0, color=GREY, lw=.8)
ax.set_yticks(yy, ["M0", "M1: + dataset", "M2: + dataset + technical PC1"])
ax.set_xlabel("ADORA2A coefficient (95% CI)"); ax.set_title("Sensitivity-model estimates")

ax = axs[1, 0]; panel(ax, "C")
l = lodo.sort_values("estimate"); yy = np.arange(len(l))
ax.errorbar(l.estimate, yy, xerr=[l.estimate-l.ci95_low, l.ci95_high-l.estimate],
            fmt="o", color=BLUE, ecolor=BLUE, capsize=2, ms=4)
ax.axvline(0, color=GREY, lw=.8); ax.set_yticks(yy, l.left_out_donor_id)
ax.set_xlabel("M0 coefficient after donor omission"); ax.set_title("M0 leave-one-donor-out")

ax = axs[1, 1]; panel(ax, "D"); ax.set_axis_off()
diag = [[r.model, int(r.residual_df), f"{r.condition_number:.1f}", f"{r.vif:.1f}"] for _, r in m.iterrows()]
tbl = ax.table(cellText=diag, colLabels=["Model", "Residual\ndf", "Condition\nno.", "VIF"],
               loc="center", cellLoc="center", colColours=["#DDEBF7"]*4)
tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.55)
for cell in tbl.get_celld().values():
    cell.set_text_props(va="center")
for col in range(4):
    tbl[(0, col)].set_height(tbl[(0, col)].get_height() * 1.35)
ax.set_title("Identification diagnostics", pad=10)
ax.text(.5, .10, "M2 is a technical identification sensitivity;\nadjusted estimates are weakly identified.",
        ha="center", transform=ax.transAxes, color=GREY, fontsize=8)
save(fig, "Fig3")


# Fig 4: separate human contexts; coefficients are not cross-dataset comparable.
fig, axs = plt.subplots(3, 2, figsize=(7.3, 8.55))
fig.subplots_adjust(left=.18, right=.96, top=.96, bottom=.08, hspace=.68, wspace=.48)
axs = axs.ravel()
d = pd.read_csv(ROOT/"05_contextual_human/Context_GSE160306.tsv", sep="\t")
ax=axs[0]; panel(ax,"A"); yy=np.arange(len(d)); ax.errorbar(d.estimate,yy,xerr=[d.estimate-d.ci95_low,d.ci95_high-d.estimate],fmt="o",color=BLUE,capsize=3); ax.axvline(0,color=GREY,lw=.8); ax.set_yticks(yy,["Stage","Stage ×\nperipheral region"]); ax.set_xlabel("Coefficient (95% CI)"); ax.set_title("GSE160306 retina")
d=pd.read_csv(ROOT/"05_contextual_human/Context_GSE60436_samples.tsv",sep="\t"); ax=axs[1]; panel(ax,"B"); order=[x for x in ["normal_retina","inactive_FVM","active_FVM"] if x in set(d.group)]
for i,g in enumerate(order):
    vals=d.loc[d.group==g,"ADORA2A"].dropna(); ax.scatter(np.full(len(vals),i)+np.linspace(-.06,.06,len(vals)),vals,color=[GREY,BLUE,ORANGE][i],s=26)
ax.set_xticks(range(len(order)),[x.replace("_","\n") for x in order]); ax.set_ylabel("Deposited expression value"); ax.set_title("GSE60436 FVM")
d=pd.read_csv(ROOT/"05_contextual_human/Context_GSE276892.tsv",sep="\t"); ax=axs[2]; panel(ax,"C"); groups=list(dict.fromkeys(d.group))
for i,g in enumerate(groups):
    vals=d.loc[d.group==g,"ADORA2A"].dropna(); ax.scatter(np.full(len(vals),i)+np.linspace(-.08,.08,len(vals)),vals,color=[BLUE,ORANGE,GREEN,MAGENTA][i%4],s=25); ax.plot([i-.15,i+.15],[np.median(vals)]*2,color=DARK)
ax.set_xticks(range(len(groups)),groups); ax.set_ylabel("ADORA2A expression"); ax.set_title("GSE276892 vitreous cells")
d=pd.read_csv(ROOT/"05_contextual_human/Context_CD31_FVM_effects.tsv",sep="\t"); ax=axs[3]; panel(ax,"D"); yy=np.arange(len(d)); ax.errorbar(d.effect,yy,xerr=[d.effect-d.lo,d.hi-d.effect],fmt="o",color=BLUE,capsize=3); ax.axvline(0,color=GREY,lw=.8); ax.set_yticks(yy,["GSE94019\nPDR FVM vs control", "GSE307925\ntreated vs untreated", "GSE179568\nPDR vs macular pucker", "GSE179568\nPDR vs macular hole"]); ax.set_xlabel("Dataset-specific coefficient (95% CI)"); ax.set_title("CD31/FVM contexts")
ax.text(0,-.30,"Coefficient magnitudes are dataset-specific\nand are not directly comparable.",transform=ax.transAxes,fontsize=8,color=GREY,va="top")
d=pd.read_csv(ROOT/"05_contextual_human/Context_2026_MNV.tsv",sep="\t"); ax=axs[4]; panel(ax,"E"); d=d[d.file=="media-2.xlsx"].drop_duplicates("cell_type").sort_values("logFC_pseudobulk"); yy=np.arange(len(d)); ax.scatter(d.logFC_pseudobulk,yy,color=ORANGE,s=26); ax.axvline(0,color=GREY,lw=.8); ax.set_yticks(yy,d.cell_type); ax.set_xlabel("MNV vs control pseudobulk logFC"); ax.set_title("2026 MNV preprint v1")
d=pd.read_csv(ROOT/"05_contextual_human/Context_GSE234047.tsv",sep="\t"); ax=axs[5]; panel(ax,"F"); wide=d.pivot(index="section",columns="stratum",values="mean_count")
if {"NON_LESIONAL","LESION_OVERLAP"}.issubset(wide.columns):
    for _,r in wide.iterrows(): ax.plot([0,1],[r.NON_LESIONAL,r.LESION_OVERLAP],color=GREY,marker="o",ms=4)
ax.set_xticks([0,1],["Non-lesional","Lesion overlap"]); ax.set_xlim(-.15,1.15); ax.set_ylabel("Mean ADORA2A count per spot"); ax.set_title("GSE234047 spatial (one donor)")
save(fig,"Fig4")


# Fig 5: explicit supported / not supported / not tested boundaries.
fig, ax = plt.subplots(figsize=(7.3, 4.55)); ax.set_axis_off()
columns = [
    (.03, .64, .94, .25, "SUPPORTED", GREEN,
     "Positive average transcript direction across eligible mouse P17 OIR\nretinal transcriptomic contexts",
     "Bounded by heterogeneity, a prediction interval crossing zero,\nand biological-unit uncertainty."),
    (.03, .36, .94, .20, "NOT SUPPORTED", ORANGE,
     "Stable, reproducible positive human transcriptomic signal",
     "Donor-level sensitivity estimates were imprecise;\ncontext-specific results were not pooled."),
    (.03, .06, .94, .22, "NOT TESTED", GREY,
     "Receptor activation · causal importance\nTherapeutic efficacy",
     "All conclusions remain at transcript level."),
]
for x,y,w,h,label,color,statement,boundary in columns:
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.012",facecolor="white",edgecolor=color,linewidth=1.4))
    ax.add_patch(Rectangle((x,y),.19,h,facecolor=color,alpha=.14,edgecolor="none"))
    ax.text(x+.095,y+h/2,label,ha="center",va="center",fontweight="bold",fontsize=9,color=color)
    ax.text(x+.22,y+h*.64,statement,ha="left",va="center",fontsize=8.5,fontweight="bold",wrap=True)
    ax.text(x+.22,y+h*.28,boundary,ha="left",va="center",fontsize=8,color=GREY,wrap=True)
save(fig,"Fig5")

print(f"Figures written to {OUT}")
