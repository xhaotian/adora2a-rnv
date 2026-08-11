#!/usr/bin/env python3
"""Generate the supplementary quantitative figure collection."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

PROJECT=Path(__file__).resolve().parents[2]
ROOT=PROJECT/"review_remediation"; OUT=ROOT/"07_figures"; OUT.mkdir(parents=True,exist_ok=True)
BLUE,ORANGE,GREEN,GREY="#0072B2","#D55E00","#009E73","#777777"
plt.rcParams.update({"font.family":"Arimo","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})

with PdfPages(OUT/"Supplementary_Figures.pdf") as pdf:
    # S1
    scr=pd.read_csv(ROOT/"03_systematic_search/screening_log.tsv",sep="\t")
    counts=scr.decision.value_counts().sort_values()
    fig,ax=plt.subplots(figsize=(8,5));ax.barh(counts.index,counts.values,color=[GREEN if "INCLUDE" in x else GREY for x in counts.index]);ax.set_xlabel("GEO series");ax.set_title("S1 Fig. GEO screening decisions");
    for i,v in enumerate(counts.values): ax.text(v+.2,i,str(v),va="center")
    fig.tight_layout();pdf.savefig(fig);plt.close(fig)

    # S2
    d=pd.read_csv(ROOT/"02_mouse_unit_audit/MOUSE_FINAL_SAMPLE_EXPRESSION.tsv",sep="\t")
    selected=d.groupby("study").size().sort_values(ascending=False).head(8).index
    d=d[d.study.isin(selected)].copy();d["display_z"]=d.groupby("study").Adora2a.transform(lambda x:(x-x.mean())/(x.std(ddof=0) or 1))
    fig,axs=plt.subplots(2,4,figsize=(11,6.2));
    for ax,(study,z) in zip(axs.ravel(),d.groupby("study")):
        for i,g in enumerate(["Normoxia","OIR"]):
            vals=z.loc[z.group==g,"display_z"].dropna();ax.scatter(np.full(len(vals),i)+np.linspace(-.08,.08,len(vals)),vals,color=[BLUE,ORANGE][i],s=24)
            if len(vals): ax.plot([i-.14,i+.14],[np.median(vals)]*2,color="#222222")
        ax.set_xticks([0,1],["Normoxia","OIR"]);ax.set_title(study);ax.axhline(0,color="#BBBBBB",lw=.6)
    fig.suptitle("S2 Fig. Selected compatible mouse sample-level examples",fontweight="bold");fig.supxlabel("Deposited retinal sample/library units");fig.supylabel("Within-study standardized ADORA2A expression");fig.tight_layout(rect=[.02,.03,1,.95]);pdf.savefig(fig);plt.close(fig)

    # S3
    mod=pd.read_csv(ROOT/"04_human_rebuild/human_primary_model_full_results.tsv",sep="\t");mod=mod[mod.term=="ADORA2A_z"].copy()
    fig,ax=plt.subplots(figsize=(8,4.8)); labels=[]
    for j,(univ,color) in enumerate([("HIGHER_CELL_COUNT_7",BLUE),("EXPANDED_DONOR_9",ORANGE)]):
        z=mod[mod.analysis_universe==univ].set_index("model").loc[["M0","M1","M2"]]; y=np.arange(3)+(j-.5)*.18
        ax.errorbar(z.estimate,y,xerr=[z.estimate-z.ci95_low,z.ci95_high-z.estimate],fmt="o",color=color,capsize=3,label=univ.replace("_"," ").title())
    ax.axvline(0,color=GREY,lw=.8);ax.set_yticks(range(3),["M0","M1","M2"]);ax.set_xlabel("ADORA2A coefficient (95% CI)");ax.set_title("S3 Fig. Higher-cell-count and expanded donor sensitivity");ax.legend(["Higher-cell-count stratum (≥20 cells)","Expanded stratum (including <20 cells)"],frameon=False);fig.tight_layout();pdf.savefig(fig);plt.close(fig)

    # S4
    l=pd.read_csv(ROOT/"04_human_rebuild/human_LODO.tsv",sep="\t");l=l[l.analysis_universe=="HIGHER_CELL_COUNT_7"]
    fig,axs=plt.subplots(1,3,figsize=(11,5.2),sharey=True)
    donors=sorted(l.left_out_donor_id.unique());pos={x:i for i,x in enumerate(donors)}
    for ax,model in zip(axs,["M0","M1","M2"]):
        z=l[l.model==model].copy(); y=z.left_out_donor_id.map(pos);ax.errorbar(z.estimate,y,xerr=[z.estimate-z.ci95_low,z.ci95_high-z.estimate],fmt="o",color=BLUE,capsize=2);ax.axvline(0,color=GREY,lw=.8);ax.set_title(model);ax.set_xlabel("Coefficient (95% CI)")
    axs[0].set_yticks(range(len(donors)),donors);fig.suptitle("S4 Fig. Full leave-one-donor-out identification diagnostics",fontweight="bold");fig.tight_layout(rect=[0,0,1,.94]);pdf.savefig(fig);plt.close(fig)

    # S5
    sig=pd.read_csv(ROOT/"04_human_rebuild/human_signature_sensitivity_results.tsv",sep="\t")
    candidates=[c for c in ["signature","signature_name","model","estimate","ci95_low","ci95_high"] if c in sig.columns]
    if {"estimate","ci95_low","ci95_high"}.issubset(sig.columns):
        z=sig.copy(); namecol="signature" if "signature" in z.columns else "signature_name"; z=z[z.get("term",pd.Series(["ADORA2A_z"]*len(z))).astype(str).str.contains("ADORA2A",case=False,na=False)] if "term" in z.columns else z
        z=z.head(30);fig,ax=plt.subplots(figsize=(9,max(4,0.25*len(z)+1.5)));y=np.arange(len(z));ax.errorbar(z.estimate,y,xerr=[z.estimate-z.ci95_low,z.ci95_high-z.estimate],fmt="o",color=BLUE,capsize=2);ax.axvline(0,color=GREY,lw=.8);labels=z[namecol].astype(str)+((" · "+z.model.astype(str)) if "model" in z.columns else "");ax.set_yticks(y,labels);ax.set_xlabel("ADORA2A coefficient (95% CI)");ax.set_title("S5 Fig. All executed donor-level signature analyses");fig.tight_layout();pdf.savefig(fig);plt.close(fig)
    else:
        fig,ax=plt.subplots(figsize=(8,3));ax.axis("off");ax.text(.5,.5,"S5 Fig. All executed donor-level signature results are provided in Supplementary Table S2.",ha="center");pdf.savefig(fig);plt.close(fig)

    # S6: exploratory small-study-effect diagnostic
    e=pd.read_csv(ROOT/"02_mouse_unit_audit/MOUSE_FINAL_PRIMARY_STUDY_EFFECTS.tsv",sep="\t")
    test=pd.read_csv(ROOT/"02_mouse_compartment_audit/FUNNEL_ASYMMETRY_DIAGNOSTIC.tsv",sep="\t").iloc[0]
    e["se"]=np.sqrt(e.sampling_variance)
    fig,ax=plt.subplots(figsize=(7.2,5.2))
    ax.scatter(e.hedges_g,e.se,color=BLUE,s=32)
    ax.axvline(e.hedges_g.mean(),color=GREY,ls="--",lw=1)
    ax.invert_yaxis();ax.set_xlabel("Hedges g");ax.set_ylabel("Standard error")
    ax.set_title("S6 Fig. Exploratory funnel plot")
    ax.text(.98,.96,f"Pustejovsky–Rodgers test\nt = {test.test_statistic_t:.3f}, P = {test.p:.3f}, k = {int(test.k)}",
            transform=ax.transAxes,ha="right",va="top",fontsize=9)
    fig.tight_layout();pdf.savefig(fig);plt.close(fig)

print(OUT/"Supplementary_Figures.pdf")
